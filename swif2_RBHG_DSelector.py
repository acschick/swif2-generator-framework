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

# Import utilities if available
try:
    from RBHG_utilities import (
        validate_config_file, extract_config_value,
        generate_workflow_name, SwifWorkflowManager
    )
except ImportError:
    print("Warning: Could not import RBHG_utilities, using fallback functions")
    
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


def process_individual_directory(config_path, dselector_path, dry_run=False):
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
    
    # If input_trees_dir doesn't exist, try to construct the actual path from MCWrapper output
    if not input_trees_dir or not os.path.exists(input_trees_dir):
        # Try actual MCWrapper output structure: /volatile/.../FFS1/<nametag>/<run_period>_<pol>_<ff>_<rad>/Simulation/root/trees/
        base_output = f"/volatile/halld/home/acschick/RBHG/{study_name}"
        actual_dir_name = f"{run_period}_{polarization}_{form_factor}_{radiation_mode}"
        input_trees_dir_try = os.path.join(base_output, nametag, actual_dir_name, "Simulation", "root", "trees")
        
        if os.path.exists(input_trees_dir_try):
            print(f"  Found actual trees at: {input_trees_dir_try}")
            input_trees_dir = input_trees_dir_try
            
            # Also update output_dir to match actual structure
            output_dir = os.path.join(base_output, nametag, actual_dir_name, "DSelector")
        else:
            print(f"  ERROR: Could not find input trees directory")
            print(f"    Config path: {input_trees_dir}")
            print(f"    Tried: {input_trees_dir_try}")
            return False
    
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
    
    if dry_run:
        print(f"  DRY RUN: Would submit DSelector jobs for {workflow_name}")
        return True
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Copy DSelector files to output directory
    dselector_dir = os.path.dirname(dselector_path)
    dselector_name = os.path.splitext(os.path.basename(dselector_path))[0]
    dselector_header = dselector_path.replace('.C', '.h')
    
    dest_selector = os.path.join(output_dir, f"{dselector_name}.C")
    dest_header = os.path.join(output_dir, f"{dselector_name}.h")
    
    shutil.copy(dselector_path, dest_selector)
    if os.path.exists(dselector_header):
        shutil.copy(dselector_header, dest_header)
    
    print(f"  Copied DSelector to {output_dir}")
    
    # Modify DSelector output filenames to be unique
    modify_dselector_output_names(dest_selector, study_name, nametag, form_factor, 
                                   run_period, polarization, lepton)
    
    # Create config file for launch script
    config_file = create_dselector_config(output_dir, workflow_name, input_trees_dir, 
                                         dest_selector, study_name)
    
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


def create_dselector_config(output_dir, workflow_name, input_dir, selector_path, study_name):
    """Create configuration file for DSelector launch script"""
    config_content = f"""# DSelector analysis configuration
# Generated by swif2_RBHG_DSelector.py

# WORKFLOW INFO
WORKFLOW           {workflow_name}
OUTDIR_LARGE       {output_dir}
OUTDIR_SMALL       {output_dir}
SELECTOR_NAME      {selector_path.replace('.C', '')}
INDATA_TOPDIR      {input_dir}

# SCICOMP JOB ACCOUNTING
PROJECT            halld
TRACK              production
OS                 el9

# JOB RESOURCES
NCORES             1
DISK               80MB
RAM                2GB
TIMELIMIT          10hrs

# JOB CONTROL
ENVFILE            /group/halld/www/halldweb/html/halld_versions/version.xml
SCRIPTFILE         /work/halld/home/acschick/channels/batch_submission/2025launch/root_analysis/script.sh

# ROOT CONFIG
ROOT_SCRIPT        /work/halld/home/acschick/channels/batch_submission/2025launch/root_analysis/Run_Selector.C
TREE_NAME          epemmissprot__B2_Tree
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
    launch_script = "/work/halld/home/acschick/channels/batch_submission/2024launch/launch/launch2.py"
    cmd_add = f"{launch_script} {config_file} 40000 57000"
    print(f"  Adding jobs: {cmd_add}")
    
    if not dry_run:
        result = subprocess.run(cmd_add, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"    ERROR adding jobs: {result.stderr}")
            return False
    
    # Run workflow
    cmd_run = f"swif2 run {workflow_name}"
    print(f"  Running workflow: {cmd_run}")
    
    if not dry_run:
        result = subprocess.run(cmd_run, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"    ERROR running workflow: {result.stderr}")
            return False
    
    print(f"  ✓ Successfully submitted {workflow_name}")
    return True


def process_run_period_master(config_path, dselector_path, dry_run=False):
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
        
        if process_individual_directory(config_full_path, dselector_path, dry_run):
            success_count += 1
    
    print(f"\n  Processed {success_count}/{total_count} configurations successfully")
    return success_count > 0


def process_ffs_master(config_path, dselector_path, dry_run=False):
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
            if process_individual_directory(dataset_config_path, dselector_path, dry_run):
                success_count += 1
        elif config_level == 2:
            if process_run_period_master(dataset_config_path, dselector_path, dry_run):
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
                       default='/work/halld/home/acschick/tmva-hpo/DSelector/DSelector_2eMissingProton.C',
                       help='Path to DSelector .C file (default: %(default)s)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be done without actually submitting jobs')
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
        success = process_individual_directory(args.config, args.dselector, args.dry_run)
    elif config_level == 2:
        success = process_run_period_master(args.config, args.dselector, args.dry_run)
    elif config_level == 3:
        success = process_ffs_master(args.config, args.dselector, args.dry_run)
    
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
