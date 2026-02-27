#!/usr/bin/env python3
"""
Merge All Histograms Script
============================
Merge histogram txt files into ROOT files for all directories in a study.

Usage:
  # Interactive (runs locally on all directories sequentially):
  ./merge_all_histograms.py output/RBHG/CPP_FFS1/ --interactive
  
  # Farm submission (one job per directory):
  ./merge_all_histograms.py output/RBHG/CPP_FFS1/ --farm
  
  # Use master JSON to auto-discover directories:
  ./merge_all_histograms.py output/RBHG/CPP_FFS1/BC_FF1/run_period_master_config.json --farm
  ./merge_all_histograms.py output/RBHG/CPP_FFS1/ffs_master_config.json --interactive
  
  # Specify cleanup mode (move/delete/none):
  ./merge_all_histograms.py output/RBHG/CPP_FFS1/ --interactive --cleanup move

This script:
- Recursively finds all directories with hists/ subdirectories
- For each found directory, merges histogram txt files into merged_histograms.root
- Can run interactively (sequential) or submit farm jobs (parallel)
- Uses mergeHists.py from generators/RBHG/histogram_scripts/
"""

import os
import sys
import json
import subprocess
import glob
from pathlib import Path


def find_histogram_directories(base_path):
    """
    Recursively find all directories containing a 'hists/' subdirectory.
    
    Args:
        base_path: Root directory to search from
        
    Returns:
        List of paths to directories containing hists/
    """
    hist_dirs = []
    
    for root, dirs, files in os.walk(base_path):
        if 'hists' in dirs:
            hists_path = os.path.join(root, 'hists')
            # Check if there are actually histogram files
            txt_files = glob.glob(os.path.join(hists_path, 'array_*.txt'))
            if txt_files:
                hist_dirs.append(root)
    
    return sorted(hist_dirs)


def load_master_config_directories(json_path):
    """
    Load directories from a master config JSON file.
    Handles both Level 2 (run_period_master_config.json) and Level 3 (ffs_master_config.json).
    
    Args:
        json_path: Path to master config JSON
        
    Returns:
        List of directory paths from the config
    """
    with open(json_path, 'r') as f:
        config = json.load(f)
    
    directories = []
    
    # Level 3: ffs_master_config.json has run_period_configs at top level
    if 'run_period_configs' in config:
        for rp_config_path in config['run_period_configs']:
            if os.path.exists(rp_config_path):
                rp_dirs = load_master_config_directories(rp_config_path)
                directories.extend(rp_dirs)
    
    # Level 2: run_period_master_config.json has individual_configs
    elif 'individual_configs' in config:
        for ind_config_path in config['individual_configs']:
            if os.path.exists(ind_config_path):
                # Extract directory from config path (remove /rbhg_config.json)
                directory = os.path.dirname(ind_config_path)
                # Check if it has histograms
                hists_path = os.path.join(directory, 'hists')
                if os.path.exists(hists_path):
                    txt_files = glob.glob(os.path.join(hists_path, 'array_*.txt'))
                    if txt_files:
                        directories.append(directory)
    
    return directories


def generate_descriptive_filename(directory):
    """
    Generate a descriptive ROOT filename from directory structure.
    
    Args:
        directory: Path to directory (e.g., output/RBHG/CPP_FFS1/BC_FFN/2205_135DEG_FFN_IRADOFF)
        
    Returns:
        Descriptive filename (e.g., CPP_FFS1_BC_FFN_2205_135DEG_IRADOFF_PRESIM.root)
    """
    # Get the path components
    parts = Path(directory).parts
    
    # Find where 'output' and 'RBHG' are, then take everything after
    try:
        rbhg_index = parts.index('RBHG')
        relevant_parts = list(parts[rbhg_index + 1:])  # Everything after RBHG/
    except (ValueError, IndexError):
        # Fallback: just use last few directory names
        relevant_parts = list(parts[-4:]) if len(parts) >= 4 else list(parts)
    
    if len(relevant_parts) >= 3:
        # Expected structure: [study_name, dataset_formfactor, individual_runperiod_pol_ff_extras]
        study_name = relevant_parts[0]
        dataset_ff = relevant_parts[1]  # e.g., BC_FFN or BC_FF1
        individual = relevant_parts[2]  # e.g., 2205_135DEG_FFN_IRADOFF
        
        # Extract form factor from dataset level (FFN or FF1)
        form_factor = None
        if dataset_ff.endswith('_FFN'):
            form_factor = 'FFN'
        elif dataset_ff.endswith('_FF1'):
            form_factor = 'FF1'
        
        # Remove redundant form factor from individual level
        if form_factor and f'_{form_factor}_' in individual:
            individual = individual.replace(f'_{form_factor}_', '_')
        elif form_factor and individual.endswith(f'_{form_factor}'):
            individual = individual[:-len(f'_{form_factor}')]
        
        # Construct filename with PRESIM label
        filename = f"{study_name}_{dataset_ff}_{individual}_PRESIM.root"
    else:
        # Fallback: just join all parts
        filename = "_".join(relevant_parts) + "_PRESIM.root"
    
    return filename


def merge_interactive(directory, cleanup_mode='none'):
    """
    Merge histograms in a directory interactively (locally).
    
    Args:
        directory: Path to directory containing hists/ subdirectory
        cleanup_mode: 'none', 'move', or 'delete'
    """
    # Import mergeHists dynamically
    framework_home = os.path.dirname(os.path.abspath(__file__))
    merge_script = os.path.join(framework_home, 'generators/RBHG/histogram_scripts/mergeHists.py')
    
    if not os.path.exists(merge_script):
        print(f"ERROR: Cannot find mergeHists.py at {merge_script}")
        return False
    
    # Import as module
    import importlib.util
    spec = importlib.util.spec_from_file_location("mergeHists", merge_script)
    mergeHists = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mergeHists)
    
    # Generate descriptive filename
    output_filename = generate_descriptive_filename(directory)
    
    print(f"Merging: {directory}")
    print(f"  Output: {output_filename}")
    try:
        mergeHists.merge_histograms(directory, 
                                    cleanup_mode=cleanup_mode if cleanup_mode != 'none' else None,
                                    output_filename=output_filename)
        print(f"  ✓ Complete")
        return True
    except Exception as e:
        print(f"  ✗ ERROR: {e}")
        return False


def submit_farm_job(directory, cleanup_mode='none'):
    """
    Submit a SLURM farm job to merge histograms for one directory.
    
    Args:
        directory: Path to directory containing hists/ subdirectory
        cleanup_mode: 'none', 'move', or 'delete'
    """
    framework_home = os.path.dirname(os.path.abspath(__file__))
    merge_script = os.path.join(framework_home, 'generators/RBHG/histogram_scripts/mergeHists.py')
    
    # Generate descriptive filename
    output_filename = generate_descriptive_filename(directory)
    
    # Determine cleanup argument
    cleanup_arg = f"cleanup_mode='{cleanup_mode}'" if cleanup_mode != 'none' else "cleanup_mode=None"
    
    # Create SLURM script
    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name=mergeHists
#SBATCH --output={directory}/mergeHists_%j.out
#SBATCH --error={directory}/mergeHists_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=500mb
#SBATCH --cpus-per-task=1
#SBATCH --partition=production

# Source the gxenv command
source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh

# Set up the environment that has ROOT built with python
gxenv 

# Execute the Python script
cd {framework_home}
python3 -c "
import sys
sys.path.insert(0, '{os.path.dirname(merge_script)}')
from mergeHists import merge_histograms
merge_histograms('{directory}', {cleanup_arg}, output_filename='{output_filename}')
"
"""
    
    # Write temporary sbatch script
    temp_script = os.path.join(directory, 'temp_merge_sbatch.sh')
    with open(temp_script, 'w') as f:
        f.write(sbatch_script)
    
    # Submit job
    try:
        result = subprocess.run(['sbatch', temp_script], capture_output=True, text=True, check=True)
        job_id = result.stdout.strip().split()[-1]
        print(f"Submitted: {directory} (Job ID: {job_id})")
        print(f"  Output: {output_filename}")
        os.remove(temp_script)
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR submitting {directory}: {e}")
        if os.path.exists(temp_script):
            os.remove(temp_script)
        return False


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Merge histogram txt files into ROOT files for all directories in a study',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive merge for all directories under CPP_FFS1/:
  ./merge_all_histograms.py output/RBHG/CPP_FFS1/ --interactive
  
  # Farm submission using master config:
  ./merge_all_histograms.py output/RBHG/CPP_FFS1/BC_FF1/run_period_master_config.json --farm
  
  # With cleanup (move txt files to archived_histograms/):
  ./merge_all_histograms.py output/RBHG/CPP_FFS1/ --interactive --cleanup move
        """
    )
    
    parser.add_argument('path', help='Directory path or master config JSON file')
    
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--interactive', action='store_true',
                           help='Run merging interactively (sequential, local)')
    mode_group.add_argument('--farm', action='store_true',
                           help='Submit farm jobs (parallel, SLURM)')
    
    parser.add_argument('--cleanup', choices=['none', 'move', 'delete'], default='none',
                       help='Cleanup mode for original txt files (default: none)')
    
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be merged without actually doing it')
    
    args = parser.parse_args()
    
    # Determine if input is JSON or directory
    if args.path.endswith('.json'):
        if not os.path.exists(args.path):
            print(f"ERROR: Config file not found: {args.path}")
            sys.exit(1)
        print(f"Loading directories from master config: {args.path}")
        directories = load_master_config_directories(args.path)
        print(f"Found {len(directories)} directories with histograms from config")
    else:
        if not os.path.isdir(args.path):
            print(f"ERROR: Directory not found: {args.path}")
            sys.exit(1)
        print(f"Searching for histogram directories under: {args.path}")
        directories = find_histogram_directories(args.path)
        print(f"Found {len(directories)} directories with histograms")
    
    if not directories:
        print("No directories with histogram files found.")
        sys.exit(0)
    
    print("\nDirectories to process:")
    for d in directories:
        print(f"  - {d}")
    
    if args.dry_run:
        print("\n[DRY RUN] No actual merging performed.")
        sys.exit(0)
    
    print(f"\nMode: {'Interactive' if args.interactive else 'Farm submission'}")
    print(f"Cleanup: {args.cleanup}")
    print()
    
    # Confirm
    if not args.dry_run:
        response = input(f"Merge histograms in each of {len(directories)} directories? [y/N]: ")
        if response.lower() != 'y':
            print("Cancelled.")
            sys.exit(0)
    
    # Process directories
    success_count = 0
    fail_count = 0
    
    print()
    for i, directory in enumerate(directories, 1):
        print(f"[{i}/{len(directories)}] ", end="", flush=True)
        if args.interactive:
            success = merge_interactive(directory, args.cleanup)
        else:  # farm
            success = submit_farm_job(directory, args.cleanup)
        
        if success:
            success_count += 1
        else:
            fail_count += 1
    
    # Summary
    print("\n" + "="*70)
    print("Merge Summary")
    print("="*70)
    print(f"Total directories: {len(directories)}")
    print(f"Successful: {success_count}")
    print(f"Failed: {fail_count}")
    
    if args.farm:
        print("\nFarm jobs submitted. Check progress with: squeue -u $USER")
        print("Job outputs will be in each directory: mergeHists_<jobid>.out")


if __name__ == "__main__":
    main()
