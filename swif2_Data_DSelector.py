#!/usr/bin/env python3
"""
swif2_Data_DSelector.py

Submits DSelector analysis jobs for real experimental data.
Reads configuration from workflows_root_analysis.config and RunPeriods.json
to automatically create SWIF2 workflows organized by run period and polarization.

Usage:
    python swif2_Data_DSelector.py [--config workflows_root_analysis.config]
    python swif2_Data_DSelector.py --dry-run
    python swif2_Data_DSelector.py --run-period 1801 --polarizations 0DEG,45DEG
    python swif2_Data_DSelector.py --run-period 2205 --polarizations 45DEG,135DEG --targets FULL,EMPTY
"""

import os
import sys
import json
import subprocess
import argparse
from pathlib import Path

# Framework paths
FRAMEWORK_DIR = Path(__file__).parent.resolve()
DEFAULT_CONFIG = FRAMEWORK_DIR / "workflows_root_analysis.config"
RUN_PERIODS_JSON = FRAMEWORK_DIR / "RunPeriods.json"
BUNDLED_LAUNCH_SCRIPT = FRAMEWORK_DIR / "DSelector_Templates" / "launch.py"
BUNDLED_SCRIPT_SH = FRAMEWORK_DIR / "DSelector_Templates" / "script.sh"
BUNDLED_RUN_SELECTOR = FRAMEWORK_DIR / "DSelector_Templates" / "Run_Selector.C"
DEFAULT_ENVFILE = "/group/halld/www/halldweb/html/halld_versions/version.xml"


def parse_config(config_path):
    """Parse workflows_root_analysis.config file"""
    config = {}
    
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                
                # Handle boolean values
                if value.upper() in ['TRUE', 'FALSE']:
                    value = value.upper() == 'TRUE'
                
                config[key] = value
    
    return config


def load_run_periods_json():
    """Load RunPeriods.json for RCDB queries and metadata"""
    if not RUN_PERIODS_JSON.exists():
        raise FileNotFoundError(f"RunPeriods.json not found at {RUN_PERIODS_JSON}")
    
    with open(RUN_PERIODS_JSON, 'r') as f:
        return json.load(f)


def expand_list(value_str, available_options):
    """Expand comma-separated list or 'ALL' keyword"""
    if value_str.upper() == 'ALL':
        return available_options
    return [v.strip() for v in value_str.split(',')]


def get_rcdb_query(run_period, polarization, run_periods_data, target_type=None):
    """
    Get RCDB query string for run period and polarization.
    launch.py will use this query to filter runs.
    
    Args:
        run_period: Run period identifier (e.g., '1801')
        polarization: Polarization string (e.g., '0DEG')
        run_periods_data: Data from RunPeriods.json
        target_type: Optional target type (e.g., 'FULL' or 'EMPTY')
        
    Returns:
        RCDB query string, or None if not available
    """
    try:
        period_data = run_periods_data[run_period]
        
        # Check if polarization data exists
        pol_data = period_data.get('Polarizations', {}).get(polarization, {})
        if not pol_data:
            print(f"    WARNING: No data for {run_period} {polarization}")
            return None
        
        # Get RCDB query template from run period level
        rcdb_query_template = period_data.get('RCDB_query', '')
        
        if not rcdb_query_template:
            print(f"    WARNING: No RCDB query template for {run_period}")
            return None
        
        # Substitute polarization angle. AMO is represented in RCDB with
        # polarization_angle==-1.0, not as a 0-degree diamond orientation.
        pol_map = {
            '0DEG': '0',
            '45DEG': '45',
            '90DEG': '90',
            '135DEG': '135',
            'AMO': '-1.0'
        }
        
        pol_deg = pol_map.get(polarization, '0')
        rcdb_query = rcdb_query_template.replace('{polDeg}', pol_deg)
        if target_type:
            rcdb_query = rcdb_query.replace('{targetType}', target_type)
        
        print(f"    RCDB query: {rcdb_query}")
        return rcdb_query
        
    except KeyError as e:
        print(f"    ERROR: Missing data in RunPeriods.json: {e}")
        return None


def get_run_range(run_period, run_periods_data):
    """
    Get run range for a run period from RunPeriods.json.
    
    Args:
        run_period: Run period identifier (e.g., '1801')
        run_periods_data: Data from RunPeriods.json
        
    Returns:
        Tuple of (min_run, max_run) or None
    """
    try:
        period_data = run_periods_data[run_period]
        run_range_str = period_data.get('run_range', '')
        
        if not run_range_str:
            print(f"    WARNING: No run_range for {run_period}")
            return None
        
        # Parse "40856 - 42559" format
        parts = run_range_str.replace('-', ' ').split()
        if len(parts) >= 2:
            return (int(parts[0]), int(parts[-1]))
        
        return None
        
    except (KeyError, ValueError) as e:
        print(f"    ERROR parsing run_range: {e}")
        return None


def get_analysis_version(config, run_period, polarization, run_periods_data):
    """Choose the analysis tree version/path key for this run period and polarization."""
    pol_key = f'ANALYSIS_VERSION_{run_period}_{polarization}'
    if pol_key in config:
        return config[pol_key]
    
    period_key = f'ANALYSIS_VERSION_{run_period}'
    if period_key in config:
        analysis_version = config[period_key]
        if '{polarization}' in analysis_version:
            return analysis_version.replace('{polarization}', polarization)
        if '{Polarization}' in analysis_version:
            return analysis_version.replace('{Polarization}', polarization)
        return analysis_version
    
    data_paths = run_periods_data.get(run_period, {}).get('tape_analysisTrees_paths', {}).get('DATA', {})
    matching_versions = [key for key in data_paths if polarization in key]
    if len(matching_versions) == 1:
        return matching_versions[0]
    
    return 'ver01'


def get_mss_path(run_period, tree_type, analysis_version, run_periods_data):
    """
    Get MSS path to merged tree files from RunPeriods.json.
    
    Args:
        run_period: Run period (e.g., '1801')
        tree_type: Tree type (e.g., 'epemmissprot__B2_U1' or 'tree_epemmissprot__B2_U1')
        analysis_version: Analysis version (e.g., 'ver17')
        run_periods_data: Data from RunPeriods.json
        
    Returns:
        Path string or None if not found
    """
    try:
        period_data = run_periods_data[run_period]
        tape_paths = period_data.get('tape_analysisTrees_paths', {})
        
        # Look under DATA category for real experimental data
        data_paths = tape_paths.get('DATA', {})
        
        if analysis_version in data_paths:
            version_trees = data_paths[analysis_version]
            
            # Try with tree_ prefix first (JSON format), then without
            tree_key = f"tree_{tree_type}" if not tree_type.startswith('tree_') else tree_type
            tree_key_alt = tree_type.replace('tree_', '', 1) if tree_type.startswith('tree_') else f"tree_{tree_type}"
            
            if tree_key in version_trees:
                mss_path = version_trees[tree_key]
            elif tree_key_alt in version_trees:
                mss_path = version_trees[tree_key_alt]
            else:
                mss_path = None
            
            if mss_path:
                # Check if path exists
                if not os.path.exists(mss_path):
                    print(f"    WARNING: MSS path in RunPeriods.json does not exist: {mss_path}")
                    return None
                
                print(f"    MSS path: {mss_path}")
                return mss_path
        
        # Fallback: construct path manually if not in RunPeriods.json
        print(f"    WARNING: No MSS path in RunPeriods.json for DATA/{analysis_version}/{tree_type}, constructing manually")
        period_map = {
            '2017': '2017-01',
            '1801': '2018-01',
            '1808': '2018-08',
            '2205': '2022-05'
        }
        
        if run_period not in period_map:
            print(f"    ERROR: Unknown run period {run_period}")
            return None
        
        period_str = period_map[run_period]
        mss_path = f"/mss/halld/RunPeriod-{period_str}/analysis/{analysis_version}/tree_{tree_type}/merged"
        
        if not os.path.exists(mss_path):
            print(f"    ERROR: Constructed MSS path does not exist: {mss_path}")
            return None
        
        print(f"    MSS path (constructed): {mss_path}")
        return mss_path
        
    except KeyError as e:
        print(f"    ERROR: Missing data in RunPeriods.json: {e}")
        return None


def get_tree_name_from_type(tree_type):
    """
    Extract the ROOT tree name from the tree type.
    
    Tree type format: "epemmissprot__B2_U1" (or with "tree_" prefix)
    Tree name format: "epemmissprot__B2_U1_Tree"
    
    Simply strip "tree_" prefix if present and append "_Tree".
    
    Args:
        tree_type: Tree type string (e.g., 'epemmissprot__B2_U1' or 'tree_epemmissprot__B2_U1')
        
    Returns:
        Tree name for the ROOT TTree object
    """
    # Strip "tree_" prefix if present
    base_type = tree_type.replace('tree_', '', 1) if tree_type.startswith('tree_') else tree_type
    
    # Construct tree name by appending "_Tree"
    return f"{base_type}_Tree"


def get_targets(config, args, run_period, run_periods_data):
    """Return target types to process, or [None] when target splitting is disabled."""
    if args.targets:
        return args.targets.split(',')
    
    if 'TARGETS' in config:
        return expand_list(config['TARGETS'], run_periods_data.get(run_period, {}).get('TargetTypes', []))
    
    target_types = run_periods_data.get(run_period, {}).get('TargetTypes', [])
    return target_types if target_types else [None]


def get_polarizations_for_period(config, args, run_period, run_periods_data):
    """Return requested polarizations, using the run-period JSON order for ALL."""
    period_pols = list(run_periods_data.get(run_period, {}).get('Polarizations', {}).keys())
    fallback_pols = ['0DEG', '45DEG', '90DEG', '135DEG', 'AMO']
    available_pols = period_pols if period_pols else fallback_pols
    
    if args.polarizations:
        return [pol.strip() for pol in args.polarizations.split(',')]
    
    pol_config = config.get('POLARIZATIONS', 'ALL')
    return expand_list(pol_config, available_pols)


def create_job_config(output_dir, log_dir, workflow_name, input_dir, 
                     selector_path, tree_name, rcdb_query, config):
    """Create jobs_root_analysis.config for launch.py"""
    
    time_limit = config.get('TIMELIMIT', '2hrs')
    cores = config.get('CORES', '1')
    ram = config.get('RAM', '4GB')
    disk = config.get('DISK', '10GB')
    account = config.get('SWIF2_ACCOUNT', 'halld')
    partition = config.get('SWIF2_PARTITION', 'production')
    
    # Convert selector_path to string and remove .C extension
    selector_name = str(selector_path).replace('.C', '')
    
    config_content = f"""# DSelector analysis configuration for real data
# Generated by swif2_Data_DSelector.py

# WORKFLOW INFO
WORKFLOW           {workflow_name}
OUTDIR_LARGE       {output_dir}
OUTDIR_SMALL       {log_dir}
SELECTOR_NAME      {selector_name}
INDATA_TOPDIR      {input_dir}

# RCDB QUERY (used by launch.py to filter runs)
RCDB_QUERY         "{rcdb_query}"

# SCICOMP JOB ACCOUNTING
PROJECT            {account}
TRACK              {partition}
OS                 el9

# JOB RESOURCES
NCORES             {cores}
DISK               {disk}
RAM                {ram}
TIMELIMIT          {time_limit}

# JOB CONTROL
ENVFILE            {DEFAULT_ENVFILE}
SCRIPTFILE         {BUNDLED_SCRIPT_SH}

# ROOT CONFIG
ROOT_SCRIPT        {BUNDLED_RUN_SELECTOR}
TREE_NAME          {tree_name}
"""
    
    config_file = output_dir / "jobs_root_analysis.config"
    config_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(config_file, 'w') as f:
        f.write(config_content)
    
    return config_file


def submit_workflow(workflow_name, config_file, run_range, dry_run=False):
    """Submit workflow using launch.py and swif2"""
    
    if dry_run:
        print(f"    DRY RUN: Would submit workflow {workflow_name}")
        print(f"      Config: {config_file}")
        print(f"      Run range: {run_range}")
        return True
    
    # Create workflow
    cmd_create = f"swif2 create {workflow_name}"
    print(f"    Creating workflow: {workflow_name}")
    
    result = subprocess.run(cmd_create, shell=True, capture_output=True, text=True)
    if result.returncode != 0 and "already exists" not in result.stderr:
        print(f"      ERROR creating workflow: {result.stderr}")
        return False
    
    # Add jobs using bundled launch.py
    cmd_add = f"python3 {BUNDLED_LAUNCH_SCRIPT} {config_file} {run_range[0]} {run_range[1]}"
    print(f"    Adding jobs: {run_range[0]}-{run_range[1]}")
    
    result = subprocess.run(cmd_add, shell=True, env=os.environ.copy())
    if result.returncode != 0:
        print(f"      ERROR: Job submission failed")
        return False
    
    # Run workflow
    cmd_run = f"swif2 run {workflow_name}"
    print(f"    Running workflow")
    
    result = subprocess.run(cmd_run, shell=True)
    if result.returncode != 0:
        print(f"      ERROR: Failed to run workflow")
        return False
    
    print(f"    ✓ Successfully submitted {workflow_name}")
    return True


def process_data_selection(config, run_periods_data, args):
    """
    Main processing logic: iterate over run periods and polarizations,
    create workflows for each combination.
    """
    
    # Determine which run periods and polarizations to process
    available_periods = list(run_periods_data.keys())
    if args.run_period:
        run_periods = [args.run_period]
    else:
        run_periods = expand_list(config.get('RUN_PERIODS', 'ALL'), available_periods)
    
    # Get other config parameters
    tree_type = config.get('TREE_TYPE', 'epemmissprot__B2_U1')
    dselector_path = Path(config.get('DSELECTOR_PATH', 
                          FRAMEWORK_DIR / 'DSelector_Templates' / '2eMissingProton_Systematics' / 'DSelector_2eMissingProton_Systematics.C'))
    output_base = Path(config.get('OUTPUT_BASE_DIR', '/volatile/halld/home/acschick/RealDataAnalysis'))
    log_base = Path(config.get('LOG_BASE_DIR', '/farm_out/acschick/DSelector_logs/RealData'))
    workflow_prefix = config.get('WORKFLOW_PREFIX', 'DSELECTOR_DATA')
    dry_run = config.get('DRY_RUN', False) or args.dry_run
    
    if not dselector_path.exists():
        print(f"ERROR: DSelector not found: {dselector_path}")
        return False
    
    print(f"\n{'='*70}")
    print(f"Processing Configuration:")
    print(f"  Run Periods: {', '.join(run_periods)}")
    if args.polarizations:
        print(f"  Polarizations: {args.polarizations}")
    else:
        print(f"  Polarizations: {config.get('POLARIZATIONS', 'ALL')} (per run period)")
    print(f"  Tree Type: {tree_type}")
    print(f"  DSelector: {dselector_path.name}")
    print(f"  Output Base: {output_base}")
    if dry_run:
        print(f"  DRY RUN MODE")
    print(f"{'='*70}\n")
    
    success_count = 0
    total_workflows = 0
    
    # Process each run period + polarization combination
    for run_period in run_periods:
        print(f"\n{'='*70}")
        print(f"Run Period: {run_period}")
        print(f"{'='*70}")
        
        targets = get_targets(config, args, run_period, run_periods_data)
        if targets != [None]:
            print(f"  Targets: {', '.join(targets)}")
        
        polarizations = get_polarizations_for_period(config, args, run_period, run_periods_data)
        
        for polarization in polarizations:
            print(f"\n  Polarization: {polarization}")
            print(f"  {'-'*66}")
            
            analysis_version = get_analysis_version(config, run_period, polarization, run_periods_data)
            print(f"    Analysis Version: {analysis_version}")
            
            for target_type in targets:
                if target_type:
                    print(f"\n    Target: {target_type}")
                    print(f"    {'-'*62}")
                
                total_workflows += 1
                
                # Get RCDB query (launch.py will use this to filter runs)
                rcdb_query = get_rcdb_query(run_period, polarization, run_periods_data, target_type)
                
                if not rcdb_query:
                    print(f"    Skipping {run_period} {polarization} (no RCDB query)")
                    continue
                
                # Get run range from RunPeriods.json
                run_range = get_run_range(run_period, run_periods_data)
                
                if not run_range:
                    print(f"    Skipping {run_period} {polarization} (no run range)")
                    continue
                
                print(f"    Run range: {run_range[0]}-{run_range[1]}")
                
                # Get MSS input path from RunPeriods.json
                mss_input_dir = get_mss_path(run_period, tree_type, analysis_version, run_periods_data)
                
                if not mss_input_dir:
                    print(f"    Skipping {run_period} {polarization} (MSS path not found)")
                    continue
                
                # Create workflow name
                workflow_parts = [workflow_prefix, run_period, polarization]
                if target_type:
                    workflow_parts.append(target_type)
                workflow_name = "_".join(workflow_parts)
                
                # Create output directories (handle {RunPeriod} placeholder)
                output_dir = Path(str(output_base).replace('{RunPeriod}', run_period)) / polarization
                log_dir = Path(str(log_base).replace('{RunPeriod}', run_period)) / polarization
                if target_type:
                    output_dir = output_dir / target_type
                    log_dir = log_dir / target_type
                
                output_dir.mkdir(parents=True, exist_ok=True)
                log_dir.mkdir(parents=True, exist_ok=True)
                
                # Get tree name from tree_type or config override
                tree_name = config.get('TREE_NAME')
                if not tree_name:
                    tree_name = get_tree_name_from_type(tree_type)
                
                print(f"    Tree name: {tree_name}")
                print(f"    Workflow: {workflow_name}")
                print(f"    Output: {output_dir}")
                
                # Create job config file
                job_config = create_job_config(
                    output_dir, log_dir, workflow_name,
                    mss_input_dir, dselector_path, tree_name, rcdb_query, config
                )
                
                # Submit workflow
                if submit_workflow(workflow_name, job_config, run_range, dry_run):
                    success_count += 1
    
    # Summary
    print(f"\n{'='*70}")
    print(f"Submission Summary:")
    print(f"  Successful: {success_count}/{total_workflows} workflows")
    print(f"{'='*70}\n")
    
    return success_count > 0


def main():
    parser = argparse.ArgumentParser(
        description='Submit DSelector analysis workflows for real experimental data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all run periods and polarizations from config
  python swif2_Data_DSelector.py
  
  # Process specific run period
  python swif2_Data_DSelector.py --run-period 1801
  
  # Process specific polarizations
  python swif2_Data_DSelector.py --polarizations 0DEG,45DEG

  # Process CPP full and empty target data separately
  python swif2_Data_DSelector.py --run-period 2205 --polarizations 45DEG,135DEG --targets FULL,EMPTY
  
  # Dry run to preview
  python swif2_Data_DSelector.py --dry-run
  
  # Use custom config
  python swif2_Data_DSelector.py --config my_workflows.config
        """
    )
    
    parser.add_argument('--config', type=Path, default=DEFAULT_CONFIG,
                       help=f'Path to workflow config file (default: {DEFAULT_CONFIG.name})')
    parser.add_argument('--run-period', 
                       help='Process specific run period (overrides config)')
    parser.add_argument('--polarizations',
                       help='Comma-separated polarizations to process (overrides config)')
    parser.add_argument('--targets',
                       help='Comma-separated target types to process (e.g. FULL,EMPTY; overrides config)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be done without submitting')
    parser.add_argument('--verbose', action='store_true',
                       help='Print detailed information')
    
    args = parser.parse_args()
    
    # Load configuration
    try:
        print(f"Loading configuration from: {args.config}")
        config = parse_config(args.config)
    except Exception as e:
        print(f"ERROR loading config: {e}")
        sys.exit(1)
    
    # Load RunPeriods.json
    try:
        run_periods_data = load_run_periods_json()
    except Exception as e:
        print(f"ERROR loading RunPeriods.json: {e}")
        sys.exit(1)
    
    # Verify bundled scripts exist
    if not BUNDLED_LAUNCH_SCRIPT.exists():
        print(f"ERROR: Bundled launch.py not found: {BUNDLED_LAUNCH_SCRIPT}")
        sys.exit(1)
    
    # Process data
    success = process_data_selection(config, run_periods_data, args)
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
