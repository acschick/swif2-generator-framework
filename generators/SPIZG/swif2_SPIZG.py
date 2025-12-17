#!/usr/bin/env python3
"""
swif2_SPIZG.py - SPIZG Generator Workflow Submission

Submits SPIZG (pi0 Primakoff) generator jobs to the farm via swif2.
Uses RunPeriods.json for automatic configuration of cobrems distributions,
run periods, and event counts.
"""

import subprocess
import os
import sys
import json
import shutil
from datetime import datetime

# Parse command line arguments for interactive mode
INTERACTIVE_MODE = '--interactive' in sys.argv or '--local' in sys.argv
DRY_RUN = '--dry-run' in sys.argv
DEBUG_MODE = '--debug' in sys.argv
if INTERACTIVE_MODE:
    print("="*70)
    print("INTERACTIVE MODE: Jobs will run locally, not submitted to swif2")
    print("="*70)
if DRY_RUN:
    print("="*70)
    print("DRY RUN MODE: Will show what would be done without executing")
    print("="*70)
if DEBUG_MODE:
    print("="*70)
    print("DEBUG MODE: Running 2 jobs with 1000 events each for testing")
    print("="*70)


#############################################################################
####################### Configuration Loading ##############################
#############################################################################

def load_config(config_file="SPIZG.config"):
    """
    Load configuration from simple key-value format file.
    Format: KEY VALUE  (separated by whitespace)
    Lines starting with # are comments and ignored
    """
    config = {}
    
    # If config_file is an absolute path, use it directly
    # Otherwise, look for it relative to this script's directory
    if os.path.isabs(config_file):
        config_path = config_file
    else:
        config_path = os.path.join(os.path.dirname(__file__), config_file)
    
    if not os.path.exists(config_path):
        print(f"ERROR: Config file {config_path} not found!")
        sys.exit(1)
    
    with open(config_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue
            
            # Strip inline comments
            if '#' in line:
                line = line.split('#')[0].strip()
            
            parts = line.split()
            if len(parts) < 2:
                if len(parts) == 1:
                    config[parts[0]] = ""
                    continue
                print(f"Warning: Invalid config line {line_num}: {line}")
                continue
            
            key = parts[0]
            value = parts[1]
            
            # Convert boolean strings
            if value.lower() in ('true', 'false'):
                config[key] = value.lower() == 'true'
            # Convert numeric strings
            elif value.replace('.', '').replace('-', '').isdigit():
                if '.' in value:
                    config[key] = float(value)
                else:
                    config[key] = int(value)
            else:
                config[key] = value
    
    # Second pass: resolve variable substitutions like {USERNAME}
    for key, value in config.items():
        if isinstance(value, str) and '{' in value and '}' in value:
            resolved_value = value
            for var_key, var_value in config.items():
                if isinstance(var_value, str):
                    resolved_value = resolved_value.replace(f'{{{var_key}}}', str(var_value))
            config[key] = resolved_value
    
    print(f"Loaded {len(config)} configuration parameters from {config_file}")
    return config


def load_run_periods(json_path="RunPeriods.json"):
    """Load RunPeriods.json for automatic configuration."""
    script_dir = os.path.dirname(__file__)
    full_path = os.path.join(script_dir, json_path)
    
    if not os.path.exists(full_path):
        print(f"Warning: {json_path} not found. Automatic lookups disabled.")
        return {}
    
    with open(full_path, 'r') as f:
        return json.load(f)


def get_cobrems_file(run_period, polarization, run_periods_data):
    """Get coherent bremsstrahlung distribution file from RunPeriods.json."""
    try:
        return run_periods_data[run_period]['Polarizations'][polarization]['cobrems_distribution']
    except KeyError:
        return None


# Load configurations
config = load_config("SPIZG.config")
run_periods = load_run_periods()

def get_config(key, default=None):
    """Get a configuration value with optional default."""
    return config.get(key, default)


def replace_path_placeholders(path_string, generator_type, username, framework_home):
    """
    Replace placeholders in path strings with actual values.
    
    Args:
        path_string: Path string potentially containing {GENERATOR_TYPE}, {USERNAME}, {FRAMEWORK_HOME}
        generator_type: Generator type (e.g., 'RBHG', 'SPIZG')
        username: Username for paths
        framework_home: Framework home directory path
        
    Returns:
        Path string with placeholders replaced
    """
    if not isinstance(path_string, str):
        return path_string
    
    replacements = {
        '{GENERATOR_TYPE}': generator_type,
        '{USERNAME}': username,
        '{FRAMEWORK_HOME}': framework_home
    }
    
    result = path_string
    for placeholder, value in replacements.items():
        result = result.replace(placeholder, value)
    
    return result


#############################################################################
############################ SPIZG Parameters ###############################
#############################################################################

# Directory setup
FRAMEWORK_HOME = get_config("FRAMEWORK_HOME", os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
TEMPLATE_DIR = os.path.join(FRAMEWORK_HOME, "generators", "SPIZG", "template_SPIZG")
USERNAME = get_config("USERNAME", os.environ.get("USER", "user"))

# Generator parameters from config
NEVENTS = get_config("NEVENTS", 50000)
FFFILE = get_config("FFFILE", "")
POL_DEG = get_config("POL_DEG", 135)
POL_FRAC = get_config("POL_FRAC", 0.99)
PRIM = get_config("PRIM", True)
STRO = get_config("STRO", True)
INTER = get_config("INTER", True)
SPIZG_PHI = get_config("SPIZG_PHI", 0.0)
SPIZG_GAMMA = get_config("SPIZG_GAMMA", 7.806)
SPIZG_C = get_config("SPIZG_C", 500)
SPIZG_ISSSIN = get_config("SPIZG_ISSSIN", True)

# Run period handling
RUN_PERIOD = str(get_config("RUN_PERIOD", "1801"))
GENERATOR_TYPE = get_config("GENERATOR_TYPE", "SPIZG")

# Output directories (with placeholder replacement)
GENERATOR_OUTPUT_BASE = replace_path_placeholders(
    get_config("GENERATOR_OUTPUT_DIR_BASE", "{FRAMEWORK_HOME}/output/{GENERATOR_TYPE}/"),
    GENERATOR_TYPE, USERNAME, FRAMEWORK_HOME
)
LOG_DIR_BASE = replace_path_placeholders(
    get_config("LOG_DIR_BASE", "/farm_out/{USERNAME}/{GENERATOR_TYPE}/"),
    GENERATOR_TYPE, USERNAME, FRAMEWORK_HOME
)

# Environment and scripts
ENV_FILE = get_config("ENV_FILE", "/group/halld/www/halldweb/html/halld_versions/version.xml")
JOB_SCRIPT = get_config("JOB_SCRIPT", "swif2_SPIZG_script.sh")
# ascii2hddm.py is centralized in hddm_scripts/
ASCII2HDDM_SCRIPT = os.path.join(FrameworkHomeDirectory, "hddm_scripts", "ascii2hddm.py")

# SWIF2 job control
PROJECT = get_config("PROJECT", "halld")
TRACK = get_config("TRACK", "production")
DISK_USAGE = get_config("DISK_USAGE", "10GB")
RAM_USAGE = get_config("RAM_USAGE", "2GB")
TIME_LIMIT = get_config("TIME_LIMIT", "6hrs")


#############################################################################
############################ Helper Functions ###############################
#############################################################################

def replace_placeholders_in_file(template_file, output_file, replacements):
    """
    Replace placeholders in template file and write to output file.
    
    Args:
        template_file: Path to template file with placeholders
        output_file: Path to output file
        replacements: Dict of {placeholder: value} pairs
    """
    with open(template_file, 'r') as f:
        content = f.read()
    
    for placeholder, value in replacements.items():
        content = content.replace(placeholder, str(value))
    
    # Ensure parent directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    with open(output_file, 'w') as f:
        f.write(content)
    
    print(f"  Created {output_file} with {len(replacements)} replacements")


def compile_generator(source_file, output_exe):
    """
    Compile SPIZG C++ generator.
    
    Args:
        source_file: Path to C++ source file
        output_exe: Path to output executable
    
    Returns:
        bool: True if compilation successful
    """
    compile_script = os.path.join(TEMPLATE_DIR, "compile_spizg.sh")
    
    if not os.path.exists(compile_script):
        print(f"ERROR: Compilation script not found: {compile_script}")
        return False
    
    print(f"  Compiling {source_file}...")
    cmd = [compile_script, source_file, output_exe]
    
    if DRY_RUN:
        print(f"  [DRY RUN] Would execute: {' '.join(cmd)}")
        return True
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"  ✓ Compilation successful: {output_exe}")
            return True
        else:
            print(f"  ✗ Compilation failed!")
            print(result.stdout)
            print(result.stderr)
            return False
    except Exception as e:
        print(f"  ✗ Compilation error: {e}")
        return False


def create_output_directories(base_dir):
    """Create necessary output directories."""
    dirs = [
        base_dir,
        os.path.join(base_dir, "vectors"),
        os.path.join(base_dir, "hddm"),
        os.path.join(base_dir, "logs")
    ]
    
    for d in dirs:
        if not os.path.exists(d):
            if not DRY_RUN:
                os.makedirs(d, exist_ok=True)
                print(f"  Created directory: {d}")
            else:
                print(f"  Would create directory: {d}")


def create_workflow_name(run_period, polarization, tag=""):
    """Generate workflow name."""
    base = f"SPIZG_{run_period}_{polarization}"
    if tag:
        base += f"_{tag}"
    return base


def create_spizg_config(run_period, polarization_deg, n_jobs, events_per_job, 
                        output_dir, vectors_path, template_dir, jana_config_path):
    """
    Create SPIZG configuration JSON compatible with prepareSimulation.py.
    Mimics the structure of rbhg_config.json created by swif2_RBHG.py.
    
    Args:
        run_period: Run period (e.g., '2205')
        polarization_deg: Polarization angle (e.g., '135DEG')
        n_jobs: Number of generation jobs
        events_per_job: Events per job
        output_dir: Output directory path
        vectors_path: Path to vectors/HDDM directory
        template_dir: Template directory path
        jana_config_path: Path to JANA analysis config file
    
    Returns:
        dict: Configuration dictionary
    """
    total_events = n_jobs * events_per_job
    study_name = f"SPIZG_{run_period}_{polarization_deg}"
    
    config = {
        "rbhg_config": {
            "study_info": {
                "study_name": study_name,
                "study_type": "SPIZGen pi0 Primakoff",
                "generator": "SPIZGen",
                "description": f"SPIZGen pi0 Primakoff generator for {run_period} ({polarization_deg})",
                "creation_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            },
            "generation_info": {
                "generator_type": "SPIZG",
                "run_period": run_period,
                "polarization": polarization_deg,
                "total_events": total_events,
                "events_per_file": events_per_job,
                "total_jobs": n_jobs,
                "mode": "run_period"
            },
            "physics_settings": {
                "reaction": "pi0_primakoff",
                "target": "Pb208",
                "beam_type": "photon"
            },
            "output_settings": {
                "ascii_output": True,
                "hddm_output": True,
                "histograms": False,
                "verbose_output": True
            },
            "directory_paths": {
                "base_paths": {
                    "framework_home": FRAMEWORK_HOME,
                    "generator_type": GENERATOR_TYPE,
                    "output_base": os.path.dirname(output_dir)
                },
                "workflow_name": study_name,
                "spizg_generation": {
                    "output_directory": output_dir,
                    "vectors_hddm_directory": vectors_path,
                    "template_directory": template_dir,
                    "farm_logs": f"generation/{study_name}"
                },
                "swif2_workflow_names": {
                    "generation": f"gen_{study_name}",
                    "simulation": f"sim_{study_name}",
                    "dselector": f"dsel_{study_name}",
                    "analysis_hists": f"hist_{study_name}"
                },
                "templates": {
                    "spizg_template_directory": template_dir,
                    "mcwrapper_template_directory": os.path.join(os.path.dirname(template_dir), "template_MCWrapper")
                }
            },
            "downstream_configs": {
                "mc_simulation": {
                    "events_per_mc_file": events_per_job,
                    "target": "Pb208",
                    "input_hddm_pattern": "vectors_*.hddm",
                    "expected_input_files": n_jobs,
                    "jana_config": jana_config_path
                }
            },
            "mcwrapper_settings": {
                "jana_config": jana_config_path,
                "run_period": run_period
            }
        }
    }
    return config


def write_spizg_config(config, output_directory):
    """Write SPIZG configuration to JSON file."""
    config_file = os.path.join(output_directory, "spizg_config.json")
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=2)
    print(f"Configuration saved to: {config_file}")
    return config_file


def run_interactive(executable, events_per_job, output_dir):
    """
    Run generator locally in interactive mode for quick debugging.
    
    Args:
        executable: Path to compiled generator executable
        events_per_job: Number of events to generate
        output_dir: Output directory (OUTPUT/2205_135DEG)
    
    Returns:
        bool: True if successful
    """
    print("\n" + "="*70)
    print("INTERACTIVE MODE: Running generator locally")
    print("="*70)
    print(f"  Executable: {executable}")
    print(f"  Events: {events_per_job}")
    print(f"  Output: {output_dir}")
    print("="*70 + "\n")
    
    # Create vectors subdirectory if it doesn't exist
    vectors_path = os.path.join(output_dir, "vectors")
    os.makedirs(vectors_path, exist_ok=True)
    
    # Change to output directory (OUTPUT/2205_135DEG)
    # The executable will create files in vectors/ subdirectory
    original_dir = os.getcwd()
    os.chdir(output_dir)
    
    try:
        # Run generator with job number 1
        cmd = [executable, "1", str(events_per_job)]
        print(f"Running: {' '.join(cmd)}")
        print("-"*70)
        
        result = subprocess.run(cmd, capture_output=False, text=True)
        
        print("-"*70)
        
        if result.returncode != 0:
            print(f"\n✗ Generator failed with exit code {result.returncode}")
            return False
        
        # Check for output files in vectors/ subdirectory
        vectors_path = os.path.join(output_dir, "vectors")
        csv_file = os.path.join(vectors_path, "genPrim.csv")
        root_file = os.path.join(vectors_path, "monteInputHistImTest2.root")
        
        if os.path.exists(csv_file):
            file_size = os.path.getsize(csv_file)
            print(f"\n✓ Created CSV file: {csv_file} ({file_size} bytes)")
        else:
            print(f"\n✗ CSV file not found: {csv_file}")
            return False
        
        if os.path.exists(root_file):
            file_size = os.path.getsize(root_file)
            print(f"✓ Created ROOT file: {root_file} ({file_size} bytes)")
        else:
            print(f"⚠ ROOT file not found: {root_file}")
        
        # Count events in CSV
        with open(csv_file, 'r') as f:
            num_lines = sum(1 for line in f)
        print(f"✓ Generated {num_lines} events")
        
        print("\n" + "="*70)
        print("✓ INTERACTIVE RUN COMPLETE")
        print("="*70)
        print(f"Output directory: {vectors_path}")
        print(f"CSV file: genPrim.csv")
        print(f"ROOT file: monteInputHistImTest2.root")
        print("="*70 + "\n")
        
        return True
        
    except Exception as e:
        print(f"\n✗ Error running generator: {e}")
        return False
    finally:
        os.chdir(original_dir)


def submit_swif2_workflow(workflow_name, job_params, n_jobs):
    """
    Create and submit swif2 workflow.
    
    Args:
        workflow_name: Name of the workflow
        job_params: Dict containing job parameters
        n_jobs: Number of jobs to create
    
    Returns:
        bool: True if successful
    """
    print(f"\n{'='*70}")
    print(f"SWIF2 Workflow: {workflow_name}")
    print(f"{'='*70}")
    print(f"  Jobs: {n_jobs}")
    print(f"  Events per job: {job_params['events_per_job']}")
    print(f"  Total events: {n_jobs * job_params['events_per_job']}")
    print(f"  Disk: {job_params['-disk']}")
    print(f"  RAM: {job_params['-ram']}")
    print(f"  Time: {job_params['-time']}")
    
    if DRY_RUN or INTERACTIVE_MODE:
        print(f"\n  [{'DRY RUN' if DRY_RUN else 'INTERACTIVE'}] Would create workflow: {workflow_name}")
        return True
    
    # Create workflow - same as RBHG
    subprocess.run(["swif2", "create", "-workflow", workflow_name])
    print(f"  ✓ Created workflow: {workflow_name}")
    
    # Add jobs - same structure as RBHG
    for job_num in range(1, n_jobs + 1):
        job_name = f"{workflow_name}_{job_num}"
        
        # Build add_job_dict similar to RBHG
        add_job_dict = {
            "-workflow": workflow_name,
            "-name": job_name,
            "-stdout": os.path.join(job_params['log_dir'], f"stdout_{job_name}.out"),
            "-stderr": os.path.join(job_params['log_dir'], f"stderr_{job_name}.err")
        }
        
        # Add swif2 resource parameters
        for key in ['-cores', '-ram', '-disk', '-time', '-account', '-partition', '-os']:
            if key in job_params:
                add_job_dict[key] = job_params[key]
        
        # Build command list
        command_list = ["swif2", "add-job"]
        
        # Append options and values (same pattern as RBHG)
        for key, value in add_job_dict.items():
            if key.startswith("-"):
                command_list.extend([key, value])
            else:
                command_list.append(value)
        
        # Add bash script path
        command_list.append(job_params['script_path'])
        
        # Arguments for bash script (same pattern as RBHG)
        command_list.append(job_params['env_file'])
        # C++ command as single string (like RBHG's fortran_command)
        cpp_command = f"{job_params['executable']} {job_num} {job_params['events_per_job']}"
        command_list.append(cpp_command)
        command_list.append(job_params['vectors_path'])
        # HDDM command as single string (same format as RBHG)
        hddm_command = f"python3 {job_params['ascii2hddm_script']} vectors_{job_num}.txt vectors_{job_num}.hddm"
        command_list.append(hddm_command)
        
        # Execute the command
        subprocess.run(command_list)
        
        # Print progress every 100 jobs
        if job_num == 1 or job_num == n_jobs or job_num % 100 == 0:
            print(f"  Adding jobs... {job_num}/{n_jobs}", end='\r')
    
    print(f"\n  ✓ Added {n_jobs} jobs to workflow")
    
    # Run workflow - same as RBHG
    subprocess.run(["swif2", "run", "-workflow", workflow_name])
    print(f"  ✓ Workflow submitted: {workflow_name}")
    
    return True


#############################################################################
############################ Main Workflow ##################################
#############################################################################

def main():
    print("\n" + "="*70)
    print("SPIZG Generator Submission Script")
    print("="*70)
    print(f"Generator Type: {GENERATOR_TYPE}")
    print(f"Run Period: {RUN_PERIOD}")
    print(f"Polarization: {POL_DEG} degrees")
    print(f"Events: {NEVENTS}")
    print(f"Form Factor File: {FFFILE}")
    print("="*70 + "\n")
    
    # Get cobrems file from RunPeriods.json
    pol_key = f"{int(POL_DEG)}DEG" if POL_DEG != "AMO" else "AMO"
    cobrems_file = get_cobrems_file(RUN_PERIOD, pol_key, run_periods)
    
    if cobrems_file:
        print(f"✓ Found cobrems distribution: {cobrems_file}")
    else:
        print(f"⚠ Warning: No cobrems file found for {RUN_PERIOD} {pol_key}")
        print(f"  Continuing with config defaults...")
    
    # Create output directory
    output_dir = os.path.join(GENERATOR_OUTPUT_BASE, f"{RUN_PERIOD}_{pol_key}")
    create_output_directories(output_dir)
    
    # Prepare generator source
    print("\n" + "="*70)
    print("Preparing Generator")
    print("="*70)
    
    template_cpp = os.path.join(TEMPLATE_DIR, "temp_genNumTot.cpp")
    working_cpp = os.path.join(output_dir, f"genNumTot_{RUN_PERIOD}_{pol_key}.cpp")
    executable = os.path.join(output_dir, f"genNumTot_{RUN_PERIOD}_{pol_key}.exe")
    
    # Placeholder replacements
    replacements = {
        "SPIZG_NUMEVENTS = 1000000": f"SPIZG_NUMEVENTS = {NEVENTS}",
        "SPIZG_POL_FRAC = 1.0": f"SPIZG_POL_FRAC = {POL_FRAC}",
        "SPIZG_POL_DEG = 135.0": f"SPIZG_POL_DEG = {float(POL_DEG)}",
        'SPIZG_FFFILE = "/w/halld-scshelf2101/home/shannen/events/generators/FF_inter/LargeTableCoulBigBig2.csv"': f'SPIZG_FFFILE = "{FFFILE}"',
        "SPIZG_PRIM = true": f"SPIZG_PRIM = {'true' if PRIM else 'false'}",
        "SPIZG_STRO = true": f"SPIZG_STRO = {'true' if STRO else 'false'}",
        "SPIZG_INTER = true": f"SPIZG_INTER = {'true' if INTER else 'false'}",
        "SPIZG_PHI = 45": f"SPIZG_PHI = {SPIZG_PHI}",
        "SPIZG_GAMMA = 7.806": f"SPIZG_GAMMA = {SPIZG_GAMMA}",
        "SPIZG_C = 500": f"SPIZG_C = {SPIZG_C}",
        "SPIZG_ISSSIN = true": f"SPIZG_ISSSIN = {'true' if SPIZG_ISSSIN else 'false'}",
        'SPIZP_OUT = "./"': f'SPIZP_OUT = "{output_dir}/vectors/"'
    }
    
    replace_placeholders_in_file(template_cpp, working_cpp, replacements)
    
    # Compile
    if not compile_generator(working_cpp, executable):
        print("\n✗ Compilation failed. Aborting.")
        sys.exit(1)
    
    # Prepare job parameters
    print("\n" + "="*70)
    print("Preparing SWIF2 Submission")
    print("="*70)
    
    # Interactive mode: run 1 job with 500 events locally
    if INTERACTIVE_MODE:
        events_per_job = 500
        n_jobs = 1
        workflow_name = create_workflow_name(RUN_PERIOD, pol_key, tag="INTERACTIVE")
        print("\n⚠ INTERACTIVE MODE: 1 job × 500 events (local execution)")
    # Debug mode: override to 2 jobs with 1000 events each
    elif DEBUG_MODE:
        events_per_job = 1000
        n_jobs = 2
        workflow_name = create_workflow_name(RUN_PERIOD, pol_key, tag="DEBUG")
        print("\n⚠ DEBUG MODE: Overriding to 2 jobs × 1000 events = 2000 total events")
    else:
        events_per_job = int(get_config("NEVENTS_PERFILE_FIXED", 10000))
        n_jobs = max(1, int(NEVENTS / events_per_job))
        workflow_name = create_workflow_name(RUN_PERIOD, pol_key)
    
    job_params = {
        '-cores': '1',
        '-ram': RAM_USAGE,
        '-disk': DISK_USAGE,
        '-time': TIME_LIMIT,
        '-account': PROJECT,
        '-partition': TRACK,
        '-os': 'el9',
        'script_path': os.path.join(TEMPLATE_DIR, JOB_SCRIPT),
        'env_file': ENV_FILE,
        'executable': executable,
        'events_per_job': events_per_job,
        'vectors_path': os.path.join(output_dir, "vectors"),
        'log_dir': os.path.join(output_dir, "logs"),
        'ascii2hddm_script': ASCII2HDDM_SCRIPT
    }
    
    # Create configuration JSON for prepareSimulation.py
    print("\n" + "="*70)
    print("Creating configuration JSON for downstream pipeline...")
    print("="*70)
    
    # Get JANA config path - use jobs_submit if it exists, otherwise allow override from config
    # JANA config - look for it in generator directory, fall back to config value
    default_jana_config = os.path.join(FRAMEWORK_HOME, "generators", "SPIZG", "jobs_submit", "jana_analysis.config")
    jana_config_path = config.get("JANA_CONFIG", default_jana_config)
    
    if not os.path.exists(jana_config_path):
        print(f"WARNING: JANA config not found: {jana_config_path}")
        print(f"  Continuing anyway - prepareSimulation.py may need manual config")
    else:
        print(f"  Using JANA config: {jana_config_path}")
    
    spizg_config = create_spizg_config(
        run_period=RUN_PERIOD,
        polarization_deg=f"{POL_DEG}DEG",
        n_jobs=n_jobs,
        events_per_job=events_per_job,
        output_dir=output_dir,
        vectors_path=os.path.join(output_dir, "vectors"),
        template_dir=TEMPLATE_DIR,
        jana_config_path=jana_config_path
    )
    
    config_json_path = write_spizg_config(spizg_config, output_dir)
    
    # Handle interactive mode - run locally instead of submitting to farm
    if INTERACTIVE_MODE:
        print("\n" + "="*70)
        print("INTERACTIVE MODE - Running generator locally")
        print("="*70)
        print(f"  Events: {events_per_job}")
        print(f"  Output: {output_dir}")
        print("="*70 + "\n")
        
        success = run_interactive(
            executable=executable,
            events_per_job=events_per_job,
            output_dir=output_dir
        )
        
        if success:
            print("✓ Interactive run completed successfully")
            sys.exit(0)
        else:
            print("✗ Interactive run failed")
            sys.exit(1)
    
    # Submit to farm
    if submit_swif2_workflow(workflow_name, job_params, n_jobs):
        print("\n" + "="*70)
        print("✓ SUCCESS")
        print("="*70)
        print(f"Workflow: {workflow_name}")
        print(f"Jobs: {n_jobs}")
        print(f"Total events: {n_jobs * events_per_job}")
        print(f"Output: {output_dir}")
        print(f"Config JSON: {config_json_path}")
        print("="*70)
        print("\nNext steps:")
        print(f"  1. Monitor jobs: swif2 status {workflow_name}")
        print(f"  2. After completion, run prepareSimulation.py:")
        print(f"     python prepareSimulation.py {config_json_path}")
        print("="*70 + "\n")
    else:
        print("\n✗ Failed to submit workflow")
        sys.exit(1)


if __name__ == "__main__":
    main()
