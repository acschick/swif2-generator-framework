#!/usr/bin/env python3
"""
swif2_RBHG_DSelector.py

Submits DSelector analysis jobs for RBHG simulation output.
Reads configuration from JSON files and submits swif2 workflows to run
DSelector over MCWrapper ReactionFilter trees.

Usage:
    python swif2_RBHG_DSelector.py <config.json>
    python swif2_RBHG_DSelector.py <run_period_master_config.json>
    python swif2_RBHG_DSelector.py <ffs_master_config.json>
"""

import os
import sys
import glob
import json
import shutil
import subprocess
import argparse
from pathlib import Path

# Framework directory (for portable paths)
FRAMEWORK_DIR = Path(__file__).parent.resolve()

# Central DSelector location (relative to framework)
CENTRAL_DSELECTOR_PATH = str(FRAMEWORK_DIR / 'DSelector_Templates' / '2eMissingProton_Systematics' / 'DSelector_2eMissingProton_Systematics.C')

# Import utilities if available
try:
    from RBHG_utilities import (
        validate_config_file, extract_config_value,
        generate_workflow_name, SwifWorkflowManager
    )
except ImportError:
    # Fallback functions (RBHG_utilities module not required)
    def extract_config_value(config_data, keys, default=None):
        """Extract value from nested dictionary using list of possible keys"""
        if not isinstance(keys, list):
            keys = [keys]
        for key_path in keys:
            try:
                value = config_data
                for key in key_path.split('.'):
                    value = value[key]
                if value is not None:
                    return value
            except (KeyError, TypeError):
                continue
        return default
    
    # Stub functions for unused imports (suppress import warning)
    def validate_config_file(*args, **kwargs):
        """Placeholder - not used in this script"""
        pass
    
    def generate_workflow_name(*args, **kwargs):
        """Placeholder - not used in this script"""
        return "workflow"
    
    class SwifWorkflowManager:
        """Placeholder - not used in this script"""
        pass


def detect_tree_name(trees_directory):
    """
    Auto-detect tree name from first ROOT file in directory.
    Uses Python ROOT bindings after gxenv setup.
    
    Args:
        trees_directory: Path to directory containing tree ROOT files
        
    Returns:
        Tree name (str) or None if detection fails
    """
    try:
        # Find first ROOT file
        tree_files = list(Path(trees_directory).glob("tree_*.root"))
        if not tree_files:
            print(f"    WARNING: No tree files found for tree name detection")
            return None
        
        sample_file = str(tree_files[0])
        print(f"    Detecting tree name from: {os.path.basename(sample_file)}")
        
        # Use subprocess with ROOT Python - inherit full environment
        cmd = f"""python3 -c "import ROOT; f = ROOT.TFile('{sample_file}'); keys = [k for k in f.GetListOfKeys() if 'TTree' in k.GetClassName()]; print(keys[0].GetName() if keys else 'NOTREE')" """
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        tree_name = result.stdout.strip()
        if tree_name and tree_name != 'NOTREE':
            print(f"    Detected tree name: {tree_name}")
            return tree_name
        else:
            # Check stderr for errors
            if result.stderr:
                print(f"    WARNING: Tree detection error: {result.stderr.strip()}")
            else:
                print(f"    WARNING: Could not detect tree name from output: '{tree_name}'")
            return None
            
    except Exception as e:
        print(f"    WARNING: Tree name auto-detection failed: {e}")
        return None


def estimate_job_time_limit(trees_directory):
    """
    Estimate appropriate time limit based on average file size.
    
    Args:
        trees_directory: Path to directory containing tree ROOT files
        
    Returns:
        Time limit string (e.g., '1hrs', '30mins')
    """
    try:
        tree_files = list(Path(trees_directory).glob("tree_*.root"))
        if not tree_files:
            return '1hrs'  # Default
        
        # Sample first 10 files to estimate average size
        sample_files = tree_files[:min(10, len(tree_files))]
        total_size = sum(os.path.getsize(f) for f in sample_files)
        avg_size_mb = (total_size / len(sample_files)) / (1024 * 1024)
        
        # Estimate: ~1-2 minutes per 20MB, with overhead
        # 20MB -> 30 mins, 50MB -> 1 hr, 100MB -> 2 hrs
        if avg_size_mb < 10:
            time_limit = '30mins'
        elif avg_size_mb < 30:
            time_limit = '1hrs'
        elif avg_size_mb < 60:
            time_limit = '2hrs'
        else:
            time_limit = '4hrs'
        
        print(f"    Average file size: {avg_size_mb:.1f} MB -> Time limit: {time_limit}")
        return time_limit
        
    except Exception as e:
        print(f"    WARNING: Could not estimate time limit: {e}")
        return '1hrs'  # Safe default


def load_and_validate_config(config_path):
    """Load JSON config and determine what type it is (Level 1/2/3)"""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, 'r', encoding='utf-8') as f:
        config_data = json.load(f)
    
    # Level 3: FFS Master
    if 'ffs_master_config' in config_data:
        return config_data, 'ffs_master', 3
    
    # Level 2: Run Period Master
    run_period_indicators = [
        'run_period_master_config' in config_data,
        'individual_configs' in config_data,
        'run_periods' in config_data and 'study_info' in config_data
    ]
    if any(run_period_indicators):
        return config_data, 'run_period_master', 2
    
    # Level 1: Individual
    individual_indicators = [
        'rbhg_config' in config_data,
        'generation_info' in config_data.get('rbhg_config', {})
    ]
    if any(individual_indicators):
        return config_data, 'individual', 1
    
    raise ValueError(f"Cannot determine configuration type from {config_path}")


def process_individual_directory(config_path, dselector_path, dry_run=False, modify_output_names=False):
    """Process a single individual directory configuration"""
    config_data, config_type, config_level = load_and_validate_config(config_path)
    
    if config_level != 1:
        print(f"  ERROR: Expected individual config, got level {config_level}")
        return False
    
    # Extract paths from config
    input_trees_dir = extract_config_value(config_data,
        ['rbhg_config.directory_paths.downstream_analysis.mcwrapper_simulation.trees_directory',
         'rbhg_config.directory_paths.downstream_analysis.simulation.trees_directory',
         'downstream_analysis.mcwrapper_simulation.trees_directory',
         'downstream_analysis.simulation.trees_directory'], None)
    
    output_dir = extract_config_value(config_data,
        ['rbhg_config.directory_paths.downstream_analysis.dselector_analysis.output_directory',
         'downstream_analysis.dselector_analysis.output_directory'], None)
    
    workflow_name = extract_config_value(config_data,
        ['rbhg_config.directory_paths.swif2_workflow_names.dselector',
         'swif2_workflow_names.dselector'], None)
    
    # Extract info for constructing actual paths if config paths are wrong
    study_name = extract_config_value(config_data,
        ['rbhg_config.generation_info.study_name', 'study_name'], 'Study')
    nametag = extract_config_value(config_data,
        ['rbhg_config.generation_info.nametag', 'nametag'], '')
    form_factor = extract_config_value(config_data,
        ['rbhg_config.physics_settings.form_factor', 'form_factor'], 'FFN')
    
    # Extract actual run period from directory path (e.g., "1801" from "1801_0DEG_FFN_DBLRAD")
    # The config's run_period field contains things like "FULL2018_AMO_FFS" which is not the actual run period
    # The directory structure is: output/Study/Nametag/RunPeriod_Pol_FF_Rad/rbhg_config.json
    config_dir = os.path.dirname(config_path)
    dir_basename = os.path.basename(config_dir)
    # Directory name format: 1801_0DEG_FFN_DBLRAD, extract first part
    run_period = dir_basename.split('_')[0] if '_' in dir_basename else '1801'
    
    polarization_raw = extract_config_value(config_data,
        ['rbhg_config.physics_settings.polarization_deg', 'polarization'], '0')
    polarization = polarization_raw if polarization_raw.upper() == 'AMO' or polarization_raw.endswith('DEG') else polarization_raw + 'DEG'
    lepton = extract_config_value(config_data,
        ['rbhg_config.physics_settings.lepton_type', 'lepton_type'], 'ee')
    radiation_mode = extract_config_value(config_data,
        ['rbhg_config.physics_settings.radiation_mode', 'radiation_mode'], 'DBLRAD')
    
    # Construct path components (always needed for log dirs)
    base_output = f"/volatile/halld/home/acschick/RBHG/{study_name}"
    actual_dir_name = f"{run_period}_{polarization}_{form_factor}_{radiation_mode}"
    
    # If input_trees_dir doesn't exist, try to construct the actual path from MCWrapper output
    if not input_trees_dir or not os.path.exists(input_trees_dir):
        # Actual volatile structure: /volatile/.../FFS1/BCFFN/1801_0DEG_FFN_DBLRAD/Simulation/root/trees/
        # Try multiple naming patterns for the middle directory:
        # 1. nametag only (e.g., "SIM")
        # 2. nametag+form_factor concatenated (e.g., "SIMFF1", "BCFFN")
        # 3. nametag_form_factor with underscore (e.g., "SIM_FF1")
        # 4. form_factor only (if no nametag)
        
        possible_middle_dirs = []
        if nametag:
            possible_middle_dirs.append(nametag)  # Try just nametag first (e.g., "SIM")
            possible_middle_dirs.append(f"{nametag}{form_factor}")  # Then concatenated (e.g., "BCFFN")
            possible_middle_dirs.append(f"{nametag}_{form_factor}")  # Then with underscore
        else:
            possible_middle_dirs.append(form_factor)
        
        input_trees_dir_try = None
        nametag_ff_combined = None
        
        for middle_dir in possible_middle_dirs:
            test_path = os.path.join(base_output, middle_dir, actual_dir_name, "Simulation", "root", "trees")
            if os.path.exists(test_path):
                input_trees_dir_try = test_path
                nametag_ff_combined = middle_dir
                print(f"  Found actual trees at: {input_trees_dir_try}")
                break
        
        if input_trees_dir_try:
            input_trees_dir = input_trees_dir_try
            
            # Also update output_dir to match actual structure
            output_dir = os.path.join(base_output, nametag_ff_combined, actual_dir_name, "DSelector")
        else:
            print(f"  ERROR: Could not find input trees directory")
            print(f"    Config path: {input_trees_dir}")
            print(f"    Tried patterns with base: {base_output}")
            print(f"      Possible middle dirs: {possible_middle_dirs}")
            print(f"      Run dir: {actual_dir_name}")
            return False
    else:
        # Config path exists, derive nametag_ff_combined for log dirs
        nametag_ff_combined = f"{nametag}{form_factor}" if nametag else form_factor
    
    if not output_dir:
        print(f"  ERROR: Could not determine output directory")
        return False
    
    # Count tree files
    tree_files = list(Path(input_trees_dir).glob("tree_*.root"))
    if len(tree_files) == 0:
        print(f"  WARNING: No tree files found in {input_trees_dir}")
        return False
    
    print(f"  Found {len(tree_files)} tree files in {input_trees_dir}")
    
    # Construct workflow name if not in config
    if not workflow_name:
        workflow_name = f"dsel_{study_name}_{nametag}_{form_factor}_{run_period}_{polarization}_{lepton}"
    
    print(f"  Workflow: {workflow_name}")
    print(f"  Input: {input_trees_dir}")
    print(f"  Output: {output_dir}")
    
    # Auto-detect tree name from input files (do before dry-run check)
    tree_name = detect_tree_name(input_trees_dir)
    if not tree_name:
        print(f"  WARNING: Using default tree name 'epemmissprot__B2_Tree'")
        tree_name = "epemmissprot__B2_Tree"
    
    # Estimate appropriate time limit based on file sizes (do before dry-run check)
    time_limit = estimate_job_time_limit(input_trees_dir)
    
    if dry_run:
        print(f"  DRY RUN: Would submit DSelector jobs for {workflow_name}")
        print(f"    Tree name: {tree_name}")
        print(f"    Time limit: {time_limit}")
        return True
    
    # Create output directory for ROOT files
    os.makedirs(output_dir, exist_ok=True)
    
    # Create log directory on /farm_out (not /volatile)
    log_dir = f"/farm_out/acschick/DSelector_logs/{study_name}/{nametag_ff_combined}/{actual_dir_name}"
    os.makedirs(log_dir, exist_ok=True)
    print(f"  Log directory: {log_dir}")
    
    # Determine DSelector location and optionally modify output names
    if modify_output_names:
        # Copy DSelector and modify output filenames for uniqueness
        if not os.path.exists(dselector_path):
            print(f"  ERROR: DSelector not found at {dselector_path}")
            return False
        
        dselector_name = os.path.splitext(os.path.basename(dselector_path))[0]
        dselector_header = dselector_path.replace('.C', '.h')
        
        dest_selector = os.path.join(output_dir, f"{dselector_name}.C")
        dest_header = os.path.join(output_dir, f"{dselector_name}.h")
        
        shutil.copy(dselector_path, dest_selector)
        if os.path.exists(dselector_header):
            shutil.copy(dselector_header, dest_header)
        
        print(f"  Copied DSelector to {output_dir}")
        
        # Modify output names
        modify_dselector_output_names(dest_selector, study_name, nametag, form_factor, 
                                      run_period, polarization, lepton)
        
        selector_path_to_use = dest_selector
    else:
        # Use central DSelector location (no copying or modification)
        if not os.path.exists(CENTRAL_DSELECTOR_PATH):
            print(f"  ERROR: Central DSelector not found at {CENTRAL_DSELECTOR_PATH}")
            return False
        selector_path_to_use = CENTRAL_DSELECTOR_PATH
        print(f"  Using central DSelector (output names from DSelector itself)")
    
    # Create config file for launch script
    config_file = create_dselector_config(output_dir, log_dir, workflow_name, input_trees_dir, 
                                         selector_path_to_use, tree_name, time_limit)
    
    print(f"  Created config: {config_file}")
    
    # Submit to swif2
    success = submit_dselector_workflow(workflow_name, config_file, dry_run=False)
    
    return success


def modify_dselector_output_names(dselector_path, study, nametag, ff, period, pol, lepton):
    """Modify DSelector file to have unique output filenames"""
    with open(dselector_path, 'r') as f:
        content = f.read()
    
    # Construct unique filenames
    hist_file = f"Hist_{study}_{nametag}_{ff}_{period}_{pol}_{lepton}.root"
    flat_file = f"Flat_{study}_{nametag}_{ff}_{period}_{pol}_{lepton}.root"
    flat_tree = f"Flat{study}_{nametag}_{ff}_{period}_{pol}_{lepton}"
    
    # Replace output filename patterns
    import re
    content = re.sub(r'dOutputFileName\s*=\s*"[^"]*"', 
                     f'dOutputFileName = "{hist_file}"', content)
    content = re.sub(r'dFlatTreeFileName\s*=\s*"[^"]*"',
                     f'dFlatTreeFileName = "{flat_file}"', content)
    content = re.sub(r'dFlatTreeName\s*=\s*"[^"]*"',
                     f'dFlatTreeName = "{flat_tree}"', content)
    
    with open(dselector_path, 'w') as f:
        f.write(content)
    
    print(f"    Modified output names: {hist_file}, {flat_file}")


def create_dselector_config(output_dir, log_dir, workflow_name, input_dir, selector_path, tree_name, time_limit='1hrs'):
    """Create configuration file for DSelector launch script"""
    config_content = f"""# DSelector analysis configuration
# Generated by swif2_RBHG_DSelector.py

# WORKFLOW INFO
WORKFLOW           {workflow_name}
OUTDIR_LARGE       {output_dir}
OUTDIR_SMALL       {log_dir}
SELECTOR_NAME      {selector_path.replace('.C', '')}
INDATA_TOPDIR      {input_dir}

# SCICOMP JOB ACCOUNTING
PROJECT            halld
TRACK              production
OS                 el9

# JOB RESOURCES
NCORES             1
DISK               500MB
RAM                3GB
TIMELIMIT          {time_limit}

# JOB CONTROL
ENVFILE            /work/halld/home/acschick/channels/batch_submission/2026launch/launch_scripts/launch/version_7.5.0.xml
SCRIPTFILE         /work/halld/home/acschick/channels/batch_submission/2026launch/launch_scripts/root_analysis/script.sh

# ROOT CONFIG
ROOT_SCRIPT        /work/halld/home/acschick/channels/batch_submission/2026launch/launch_scripts/root_analysis/Run_Selector.C
TREE_NAME          {tree_name}
"""
    
    config_file = os.path.join(output_dir, "RBHG_dselector.config")
    with open(config_file, 'w') as f:
        f.write(config_content)
    
    return config_file


def submit_dselector_workflow(workflow_name, config_file, dry_run=False):
    """Submit DSelector workflow to swif2"""
    # Create workflow
    cmd_create = f"swif2 create {workflow_name}"
    print(f"  Creating workflow: {cmd_create}")
    
    if not dry_run:
        result = subprocess.run(cmd_create, shell=True, capture_output=True, text=True)
        if result.returncode != 0 and "already exists" not in result.stderr:
            print(f"    ERROR creating workflow: {result.stderr}")
            return False
    
    # Add jobs using launch script
    launch_script = "/work/halld/home/acschick/channels/batch_submission/2026launch/launch_scripts/launch/launch.py"
    cmd_add = f"{launch_script} {config_file} 40000 57000"
    print(f"  Adding jobs: {cmd_add}")
    print()  # Blank line before launch.py output
    
    if not dry_run:
        result = subprocess.run(cmd_add, shell=True, env=os.environ.copy())
        if result.returncode != 0:
            print(f"    ERROR: Job submission failed with exit code {result.returncode}")
            return False
        print()  # Blank line after launch.py output
    
    # Run workflow
    cmd_run = f"swif2 run {workflow_name}"
    print(f"  Running workflow: {cmd_run}")
    
    if not dry_run:
        result = subprocess.run(cmd_run, shell=True)
        if result.returncode != 0:
            print(f"    ERROR: Running workflow failed with exit code {result.returncode}")
            return False
    
    print(f"  ✓ Successfully submitted {workflow_name}")
    return True


def process_run_period_master(config_path, dselector_path, dry_run=False, modify_output_names=False):
    """Process a run period master configuration"""
    config_data, _, _ = load_and_validate_config(config_path)
    
    master_config = config_data.get('run_period_master_config', config_data)
    batch_processing = master_config.get('batch_processing', {})
    individual_configs = batch_processing.get('individual_configs', [])
    
    if not individual_configs:
        print(f"  ERROR: No individual_configs found in run period master")
        return False
    
    config_dir = os.path.dirname(os.path.abspath(config_path))
    success_count = 0
    
    # Expand glob patterns in individual_configs
    expanded_configs = []
    for config_rel_path in individual_configs:
        config_pattern = os.path.join(config_dir, config_rel_path)
        matches = glob.glob(config_pattern)
        if matches:
            expanded_configs.extend(matches)
        else:
            print(f"  WARNING: No files matched pattern: {config_rel_path}")
    
    total_count = len(expanded_configs)
    print(f"  Processing {total_count} individual configurations...")
    
    for config_full_path in expanded_configs:
        config_rel = os.path.relpath(config_full_path, config_dir)
        print(f"\n  Processing: {config_rel}")
        
        if process_individual_directory(config_full_path, dselector_path, dry_run, modify_output_names):
            success_count += 1
    
    print(f"\n  Processed {success_count}/{total_count} configurations successfully")
    return success_count > 0


def process_ffs_master(config_path, dselector_path, dry_run=False, modify_output_names=False):
    """Process an FFS master configuration"""
    config_data, _, _ = load_and_validate_config(config_path)
    
    ffs_config = config_data.get('ffs_master_config', config_data)
    datasets = ffs_config.get('datasets', {})
    
    if not datasets:
        print(f"  ERROR: No datasets found in FFS master")
        return False
    
    config_dir = os.path.dirname(os.path.abspath(config_path))
    success_count = 0
    total_count = len(datasets)
    
    print(f"  Processing {total_count} datasets...")
    
    for dataset_name, dataset_info in datasets.items():
        print(f"\n  Processing dataset: {dataset_name}")
        
        config_file = dataset_info.get('config_file') or dataset_info.get('run_period_master_config')
        if not config_file:
            print(f"    ERROR: No config_file found for dataset {dataset_name}")
            continue
        
        dataset_config_path = os.path.join(config_dir, config_file)
        
        if not os.path.exists(dataset_config_path):
            print(f"    ERROR: Config file not found: {dataset_config_path}")
            continue
        
        # Check if it's Level 1 or Level 2
        dataset_config_data, _, config_level = load_and_validate_config(dataset_config_path)
        
        if config_level == 1:
            if process_individual_directory(dataset_config_path, dselector_path, dry_run, modify_output_names):
                success_count += 1
        elif config_level == 2:
            if process_run_period_master(dataset_config_path, dselector_path, dry_run, modify_output_names):
                success_count += 1
    
    print(f"\n  Processed {success_count}/{total_count} datasets successfully")
    return success_count > 0


def main():
    parser = argparse.ArgumentParser(
        description='Submit DSelector analysis jobs for RBHG simulation output',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process single directory
  python swif2_RBHG_DSelector.py output/FFS1/qDATAq_FFN/1801_0DEG_FFN_DBLRAD/rbhg_config.json
  
  # Process run period master
  python swif2_RBHG_DSelector.py output/FFS1/qDATAq_FFN/run_period_master_config.json
  
  # Process entire FFS study
  python swif2_RBHG_DSelector.py output/FFS1/ffs_master_config.json
  
  # Dry run to see what would be submitted
  python swif2_RBHG_DSelector.py output/FFS1/ffs_master_config.json --dry-run
        """
    )
    
    parser.add_argument('config', help='Path to configuration JSON file')
    parser.add_argument('--dselector', 
                       default=str(FRAMEWORK_DIR / 'DSelector_Templates' / '2eMissingProton_Systematics' / 'DSelector_2eMissingProton_Systematics.C'),
                       help='Path to DSelector .C file (default: %(default)s)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be done without actually submitting jobs')
    parser.add_argument('--modify-output-names', action='store_true',
                       help='Copy and modify DSelector to set unique output filenames')
    parser.add_argument('--verbose', action='store_true',
                       help='Print detailed information')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.config):
        print(f"ERROR: Config file not found: {args.config}")
        sys.exit(1)
    
    if not os.path.exists(args.dselector):
        print(f"ERROR: DSelector file not found: {args.dselector}")
        sys.exit(1)
    
    print(f"RBHG DSelector Submission")
    print(f"{'='*60}")
    print(f"Config: {args.config}")
    print(f"DSelector: {args.dselector}")
    if args.dry_run:
        print(f"DRY RUN MODE - No jobs will be submitted")
    print(f"{'='*60}\n")
    
    # Load and determine config type
    try:
        config_data, config_type, config_level = load_and_validate_config(args.config)
        print(f"Configuration type: {config_type} (Level {config_level})\n")
    except Exception as e:
        print(f"ERROR loading config: {e}")
        sys.exit(1)
    
    # Process based on level
    success = False
    if config_level == 1:
        success = process_individual_directory(args.config, args.dselector, args.dry_run, args.modify_output_names)
    elif config_level == 2:
        success = process_run_period_master(args.config, args.dselector, args.dry_run, args.modify_output_names)
    elif config_level == 3:
        success = process_ffs_master(args.config, args.dselector, args.dry_run, args.modify_output_names)
    
    if success:
        print(f"\n{'='*60}")
        print(f"SUCCESS: DSelector submission completed")
        print(f"{'='*60}")
        sys.exit(0)
    else:
        print(f"\n{'='*60}")
        print(f"FAILED: DSelector submission encountered errors")
        print(f"{'='*60}")
        sys.exit(1)


if __name__ == "__main__":
    main()
