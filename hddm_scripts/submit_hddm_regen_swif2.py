#!/usr/bin/env python3
"""
submit_hddm_regen_swif2.py

Submit swif2 farm jobs to batch-convert RBHG txt → HDDM files in parallel.
All parameters (particle type, target, run number, vertex, file count, paths)
are read automatically from rbhg_config.json files found under --study-path.
No hardcoded dataset lists — works for CPP, GlueX, electrons, muons, etc.

Usage:
    python3 submit_hddm_regen_swif2.py [options]

Options:
    --study-path PATH       Root directory to scan for rbhg_config.json files
                            (default: output/RBHG/CPP_FFS1)
    --dry-run               Preview jobs without submitting anything
    --files-per-job N       Conversions per swif2 job  (default: 20)
    --workflow NAME         Override auto-generated workflow name
    --env-file PATH         GlueX environment XML      (default: version.xml)
    --skip-missing          Skip txt files that don't exist (default: abort)
    --force                 Queue files even if output HDDM already exists

Examples:
    python3 submit_hddm_regen_swif2.py --dry-run
    python3 submit_hddm_regen_swif2.py --force --files-per-job 20
    python3 submit_hddm_regen_swif2.py --study-path output/RBHG/CPP_FFS1_muons --dry-run
    python3 submit_hddm_regen_swif2.py --study-path output/RBHG/GlueX_study --force
"""

import subprocess
import os
import sys
import math
import json
import argparse
from datetime import datetime

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
FRAMEWORK          = "/work/halld/home/acschick/channels/swif2-generator-framework"
HDDM_SCRIPTS       = os.path.join(FRAMEWORK, "hddm_scripts")
ASCII2HDDM         = os.path.join(HDDM_SCRIPTS, "ascii2hddm.py")
WORKER_SCRIPT      = os.path.join(HDDM_SCRIPTS, "hddm_batch_worker.sh")
DEFAULT_STUDY_PATH  = os.path.join(FRAMEWORK, "output", "RBHG", "CPP_FFS1")
DEFAULT_ENV_FILE    = "/group/halld/www/halldweb/html/halld_versions/version.xml"
RUN_PERIODS_FILE    = os.path.join(FRAMEWORK, "RunPeriods.json")

# ---------------------------------------------------------------------------
# swif2 resource parameters for ascii2hddm conversion jobs.
# ascii2hddm.py is text→binary only, so 500 MB RAM is ample.
# ---------------------------------------------------------------------------
JOB_RESOURCES = {
    "-cores": "1",
    "-ram":   "500m",    # 0.5 GB -- text->binary conversion only
    "-disk":  "2g",
    "-time":  "60m",
    "-os":    "el9",
}

# ---------------------------------------------------------------------------
# Mappings read from rbhg_config.json
# ---------------------------------------------------------------------------
LEPTON_TO_PARTICLE = {
    "ee": "epem",
    "mu": "mupmum",
}

# Vertex string passed to ascii2hddm.py --vertex vx:vy:zmin:zmax (colon-separated, no spaces).
# "0:0:0:0" means "omit the flag entirely" (let hdgeant4 pull from CCDB).
EXPERIMENT_TO_VERTEX = {
    "CPP":   "0:0:1:1.02806",   # lead-foil target in Hall D
    "GlueX": "0:0:0:0",         # hydrogen / CCDB default
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def safe_job_name(s, maxlen=72):
    """Sanitize a string for use as a swif2 job/workflow name."""
    return s.replace("/", "_").replace(" ", "_").replace("-", "_")[:maxlen]


def workflow_exists(name):
    result = subprocess.run(
        ["swif2", "status", "-workflow", name],
        capture_output=True, text=True
    )
    return result.returncode == 0


def create_workflow(name, dry_run):
    """Create swif2 workflow; append _1, _2, ... if name already taken."""
    original = name
    n = 1
    while workflow_exists(name):
        print(f"  Workflow '{name}' already exists.")
        name = f"{original}_{n}"
        n += 1
        if n > 20:
            sys.exit("ERROR: Too many existing workflows with similar names.")
    if dry_run:
        print(f"  [DRY RUN] swif2 create -workflow {name}")
    else:
        result = subprocess.run(
            ["swif2", "create", "-workflow", name],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            print(f"WARNING: swif2 create may have failed:\n{result.stderr}")
    return name


def submit_job(workflow, job_name, logdir, worker_script, env_file, ascii2hddm,
               particle, target, run_number, vertex, txt_files, hddm_files, dry_run):
    """Build and submit (or print) one swif2 add-job command."""
    stdout_path = os.path.join(logdir, f"{job_name}.out")
    stderr_path = os.path.join(logdir, f"{job_name}.err")

    cmd = ["swif2", "add-job",
           "-workflow", workflow,
           "-name",     job_name,
           "-stdout",   stdout_path,
           "-stderr",   stderr_path]

    for flag, val in JOB_RESOURCES.items():
        cmd += [flag, val]

    # Worker: bash hddm_batch_worker.sh <envfile> <ascii2hddm> <particle> <target> <run>
    #         [--vertex vx vy zmin zmax]  txt_0 hddm_0  txt_1 hddm_1 ...
    cmd += ["bash", worker_script, env_file, ascii2hddm, particle, target, str(run_number)]

    if vertex and vertex.strip() != "0:0:0:0":
        cmd += ["--vertex", vertex]

    for txt, hddm in zip(txt_files, hddm_files):
        cmd += [txt, hddm]

    if dry_run:
        print(f"    [DRY RUN] -name {job_name}")
        print(f"              particle={particle} target={target} run={run_number} vertex='{vertex}'")
        print(f"              files: {os.path.basename(txt_files[0])}...{os.path.basename(txt_files[-1])} ({len(txt_files)})")
        return True

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR adding job {job_name}: {result.stderr.strip()}")
        return False
    return True


def extract(d, *keys, default=None):
    """Safely traverse a nested dict: extract(cfg, 'a', 'b', 'c')."""
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    return d


_run_periods_cache = None

def _load_run_periods():
    global _run_periods_cache
    if _run_periods_cache is None:
        try:
            with open(RUN_PERIODS_FILE) as f:
                _run_periods_cache = json.load(f)
        except Exception as e:
            print(f"  WARNING: Could not load RunPeriods.json: {e}")
            _run_periods_cache = {}
    return _run_periods_cache


def lookup_run_number(run_period_str, polarization_deg):
    """
    Fall back to RunPeriods.json when rbhg_config has run_number = null.
    run_period_str : e.g. "2205_AMO_FFS"  → period key "2205"
    polarization_deg: e.g. "135", "45", "AMO"  → key "135DEG", "45DEG", "AMO"
    Returns int run number or None.
    """
    rp = _load_run_periods()
    period_key = run_period_str.split("_")[0] if run_period_str else None
    if not period_key or period_key not in rp:
        return None
    pol_key = polarization_deg if polarization_deg == "AMO" else f"{polarization_deg}DEG"
    return extract(rp, period_key, "Polarizations", pol_key, "characteristic_run")


def load_dataset_info(config_path):
    """
    Read one rbhg_config.json and return a dict with all parameters needed
    to submit HDDM regeneration jobs.  Returns None on error.
    """
    try:
        with open(config_path) as f:
            raw = json.load(f)
    except Exception as e:
        print(f"  ERROR reading {config_path}: {e}")
        return None

    cfg = raw.get("rbhg_config", raw)   # handle both wrapped and flat forms

    experiment    = extract(cfg, "physics_settings", "experiment",      default="GlueX")
    lepton_type   = extract(cfg, "physics_settings", "lepton_type",     default="ee")
    target        = extract(cfg, "physics_settings", "target",          default=None)
    run_number    = extract(cfg, "physics_settings", "run_number",      default=None)
    run_period    = extract(cfg, "physics_settings", "run_period",      default=None)
    polarization  = extract(cfg, "physics_settings", "polarization_deg",default=None)
    total_jobs    = extract(cfg, "event_counts",     "total_jobs",      default=None)
    vectors_dir   = extract(cfg, "directory_paths", "rbhg_generation",
                            "vectors_hddm_directory", default=None)

    # run_number is often null in generated configs; fall back to RunPeriods.json
    if run_number is None and run_period and polarization:
        run_number = lookup_run_number(run_period, str(polarization))
        if run_number is not None:
            print(f"    [run_number] looked up from RunPeriods.json: {run_number}")

    missing = [name for name, val in [("target",      target),
                                       ("run_number",  run_number),
                                       ("total_jobs",  total_jobs),
                                       ("vectors_dir", vectors_dir)]
               if val is None]
    if missing:
        print(f"  WARNING: {config_path} is missing fields: {missing}")
        return None

    particle = LEPTON_TO_PARTICLE.get(lepton_type)
    if particle is None:
        print(f"  WARNING: unknown lepton_type '{lepton_type}' in {config_path}")
        return None

    return {
        "config_path": config_path,
        "vectors_dir": vectors_dir,
        "particle":    particle,
        "target":      target,
        "run_number":  int(run_number),
        "total_jobs":  int(total_jobs),
        "experiment":  experiment,
        "lepton_type": lepton_type,
        "vertex":      EXPERIMENT_TO_VERTEX.get(experiment, "0:0:0:0"),
    }


def find_configs(study_path):
    """
    Recursively find all rbhg_config.json files under study_path.
    Skips MCWrapper, logs, __pycache__, and FortranFiles subdirectories.
    """
    study_path = os.path.abspath(study_path)
    if os.path.isfile(study_path):
        return [study_path]
    if not os.path.isdir(study_path):
        print(f"ERROR: study path not found: '{study_path}'")
        sys.exit(1)

    SKIP_DIRS = {"MCWrapper", "logs", "__pycache__", "FortranFiles", "hists"}
    configs = []
    for root, dirs, files in os.walk(study_path):
        dirs[:] = sorted([d for d in dirs if d not in SKIP_DIRS])
        if "rbhg_config.json" in files:
            configs.append(os.path.join(root, "rbhg_config.json"))
    return sorted(configs)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Submit swif2 jobs to convert RBHG txt→HDDM (config-driven)"
    )
    parser.add_argument("--study-path",    default=DEFAULT_STUDY_PATH, metavar="PATH",
                        help="Root path to scan for rbhg_config.json files "
                             f"(default: output/RBHG/CPP_FFS1)")
    parser.add_argument("--dry-run",       action="store_true",
                        help="Preview only; nothing is submitted or modified")
    parser.add_argument("--files-per-job", type=int, default=20, metavar="N",
                        help="Number of txt→hddm conversions per swif2 job (default: 20)")
    parser.add_argument("--workflow",      default=None, metavar="NAME",
                        help="Override auto-generated workflow name")
    parser.add_argument("--env-file",      default=DEFAULT_ENV_FILE, metavar="PATH",
                        help=f"GlueX environment XML (default: version.xml)")
    parser.add_argument("--skip-missing",  action="store_true",
                        help="Skip missing txt files instead of aborting")
    parser.add_argument("--force",         action="store_true",
                        help="Queue files even if output HDDM already exists")
    args = parser.parse_args()

    timestamp   = datetime.now().strftime("%Y%m%d_%H%M")
    study_label = os.path.basename(os.path.normpath(args.study_path))
    workflow_name = args.workflow or f"hddm_regen_{study_label}_{timestamp}"

    logdir = os.path.join(FRAMEWORK, "output", "RBHG", "logs", "hddm_regen", timestamp)

    print("=" * 68)
    print("  submit_hddm_regen_swif2.py  (config-driven)")
    print(f"  Study path:    {args.study_path}")
    print(f"  Workflow:      {workflow_name}{' (DRY RUN)' if args.dry_run else ''}")
    print(f"  Files/job:     {args.files_per_job}   RAM: {JOB_RESOURCES['-ram']}")
    print(f"  Env file:      {args.env_file}")
    print(f"  Log directory: {logdir}")
    print("=" * 68)

    # Sanity checks
    for path, label in [(ASCII2HDDM, "ascii2hddm.py"), (WORKER_SCRIPT, "hddm_batch_worker.sh")]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    if not args.dry_run:
        os.makedirs(logdir, exist_ok=True)
    else:
        print(f"  [DRY RUN] Would create log directory: {logdir}")

    # Discover all configs under study path
    config_paths = find_configs(args.study_path)
    if not config_paths:
        print(f"ERROR: No rbhg_config.json found under {args.study_path}")
        sys.exit(1)
    print(f"\nFound {len(config_paths)} config(s):")
    for p in config_paths:
        print(f"  {os.path.relpath(p, FRAMEWORK)}")

    # Create workflow
    print(f"\nCreating swif2 workflow: {workflow_name}")
    workflow_name = create_workflow(workflow_name, args.dry_run)
    print(f"  -> {workflow_name}")

    # ---------------------------------------------------------------------------
    # Process each config / dataset
    # ---------------------------------------------------------------------------
    total_jobs_submitted = 0
    total_files_queued   = 0
    all_errors           = []

    for config_path in config_paths:
        rel = os.path.relpath(config_path, FRAMEWORK)
        print(f"\n--- {rel} ---")

        info = load_dataset_info(config_path)
        if info is None:
            all_errors.append(f"Failed to load: {config_path}")
            continue

        vectors_dir = info["vectors_dir"]
        n_files     = info["total_jobs"]

        print(f"    experiment={info['experiment']}  particle={info['particle']}  "
              f"target={info['target']}")
        print(f"    run={info['run_number']}  vertex='{info['vertex']}'  "
              f"total_jobs={n_files}")
        print(f"    vectors_dir: {vectors_dir}")

        if not os.path.isdir(vectors_dir):
            msg = f"  WARNING: vectors/ not found: {vectors_dir}"
            print(msg)
            all_errors.append(msg)
            continue

        # Collect (txt, hddm) pairs
        pairs = []
        abort = False
        for i in range(n_files):
            txt_file  = os.path.join(vectors_dir, f"vectors_{i}.txt")
            hddm_file = os.path.join(vectors_dir, f"vectors_{i}.hddm")

            if not os.path.exists(txt_file):
                msg = f"  {'WARNING' if args.skip_missing else 'ERROR'}: missing vectors_{i}.txt"
                print(msg)
                if args.skip_missing:
                    all_errors.append(msg)
                    continue
                else:
                    abort = True
                    all_errors.append(msg)
                    break

            if os.path.exists(hddm_file) and not args.force:
                continue   # already converted

            pairs.append((txt_file, hddm_file))

        if abort:
            print(f"  Aborting dataset. Use --skip-missing to continue past gaps.")
            continue

        if not pairs:
            print(f"  Nothing to convert (all HDDM already exist; use --force to overwrite).")
            continue

        n_batches  = math.ceil(len(pairs) / args.files_per_job)
        dir_label  = safe_job_name(
            os.path.relpath(os.path.dirname(config_path),
                            os.path.join(FRAMEWORK, "output", "RBHG")),
            maxlen=48
        )
        print(f"    Files: {len(pairs)}  ->  {n_batches} jobs x {args.files_per_job}/job")

        for batch_idx in range(n_batches):
            batch = pairs[batch_idx * args.files_per_job:
                          (batch_idx + 1) * args.files_per_job]
            job_name = safe_job_name(f"hddm_regen_{dir_label}_b{batch_idx:04d}")

            ok = submit_job(
                workflow      = workflow_name,
                job_name      = job_name,
                logdir        = logdir,
                worker_script = WORKER_SCRIPT,
                env_file      = args.env_file,
                ascii2hddm    = ASCII2HDDM,
                particle      = info["particle"],
                target        = info["target"],
                run_number    = info["run_number"],
                vertex        = info["vertex"],
                txt_files     = [p[0] for p in batch],
                hddm_files    = [p[1] for p in batch],
                dry_run       = args.dry_run,
            )
            if ok:
                total_jobs_submitted += 1
                total_files_queued   += len(batch)

        print(f"    check: {n_batches} jobs queued for this dataset")

    # ---------------------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------------------
    print("\n" + "=" * 68)
    print("  SUMMARY")
    print("=" * 68)
    print(f"  Workflow:       {workflow_name}")
    print(f"  Jobs submitted: {total_jobs_submitted}")
    print(f"  Files queued:   {total_files_queued}")

    if all_errors:
        print(f"\n  Warnings/errors ({len(all_errors)}):")
        for e in all_errors:
            print(f"    {e}")

    if not args.dry_run and total_jobs_submitted > 0:
        print(f"\n  To run the workflow:")
        print(f"    swif2 run -workflow {workflow_name}")
        print(f"\n  To monitor status:")
        print(f"    swif2 status -workflow {workflow_name}")
        print(f"    swif2 status -workflow {workflow_name} -jobs")
    elif args.dry_run:
        print("\n  DRY RUN complete. Re-run without --dry-run to submit.")

    print("=" * 68)


if __name__ == "__main__":
    main()
