#!/usr/bin/env python3

"""
New JSON-driven prepareSimulation.py
Handles any level of JSON configuration (Level 1, 2, or 3)
Includes event counting validation and completion status tracking
"""

import os
import sys
import json
import argparse
from pathlib import Path
import math
# Add framework directories to Python path for imports
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)  # For swif2_common.py
sys.path.insert(0, os.path.join(script_dir, 'hddm_scripts'))  # For hddm_event_utils.py

# Import our utility modules
try:
    from swif2_common import (
        create_swif2_workflow, add_swif2_job, submit_workflow,
        generate_workflow_name, SwifWorkflowManager
    )
    from hddm_event_utils import (
        merge_and_validate_hddm, count_directory_hddm_events,
        save_event_count_report
    )
    print("Successfully imported utility modules")
except ImportError as e:
    print(f"ERROR: Could not import utility modules: {e}")
    print("Make sure swif2_common.py and hddm_event_utils.py are available")
    print("Run 'gxenv' before running this script if HDDM tools are not available")
    sys.exit(1)

def extract_config_value(config_data, field_paths, default=None):
    """
    Flexibly extract values from JSON config using multiple possible paths
    
    Args:
        config_data (dict): Configuration dictionary
        field_paths (list): List of possible field paths to try (supports dot notation)
        default: Default value if none of the paths work
        
    Returns:
        The first successfully extracted value, or default
    """
    for path in field_paths:
        try:
            # Handle dot notation like 'directory_paths.rbhg_generation.vectors_hddm_directory'
            if '.' in path:
                value = config_data
                for part in path.split('.'):
                    value = value[part]
                if value is not None:
                    return value
            else:
                # Simple field access
                value = config_data.get(path)
                if value is not None:
                    return value
        except (KeyError, TypeError):
            continue
    
    return default

def load_and_validate_config(config_path):
    """
    Load JSON config and determine what type it is (Level 1/2/3)
    
    Args:
        config_path (str): Path to JSON configuration file
        
    Returns:
        tuple: (config_data, config_type, config_level)
            config_type: 'individual', 'run_period_master', 'ffs_master'
            config_level: 1, 2, or 3
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config_data = json.load(f)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in {config_path}: {e}")
    
    # Flexible configuration type detection
    # Check in order from most specific to least specific to avoid false positives
    
    # Level 3: FFS Master - Look for ffs_master_config key (most specific)
    if 'ffs_master_config' in config_data:
        return config_data, 'ffs_master', 3
    
    # Level 2: Run Period Master - Look for run_period_master_config key or coordination patterns
    run_period_indicators = [
        'run_period_master_config' in config_data,
        'individual_configs' in config_data,
        'run_periods' in config_data and 'study_info' in config_data,
        'coordination_summary' in config_data
    ]
    if any(run_period_indicators):
        return config_data, 'run_period_master', 2
    
    # Level 1: Individual - Look for individual directory patterns  
    individual_indicators = [
        'rbhg_config' in config_data,
        'generation_info' in config_data.get('rbhg_config', {}),
        'physics_settings' in config_data.get('rbhg_config', {}),
        'directory_paths' in config_data.get('rbhg_config', {}),
        'downstream_analysis' in config_data.get('rbhg_config', {}).get('directory_paths', {}),
        any(key in config_data for key in ['directory_name', 'vectors_hddm_directory', 'output_directory'])
    ]
    if any(individual_indicators):
        return config_data, 'individual', 1
    
    # If we can't determine the type, provide helpful error
    available_keys = list(config_data.keys())
    raise ValueError(f"Cannot determine configuration type from {config_path}\n"
                   f"Available top-level keys: {available_keys}\n"
                   f"Expected patterns for Level 1: rbhg_config (with generation_info, physics_settings, directory_paths)\n"
                   f"Expected patterns for Level 2: run_period_master_config, individual_configs\n"
                   f"Expected patterns for Level 3: ffs_master_config, datasets")

def update_completion_status(config_path, step, status, details=None):
    """
    Update JSON configuration with completion status
    
    Args:
        config_path (str): Path to JSON configuration file
        step (str): Pipeline step name ('generation', 'preparation', 'simulation', etc.)
        status (str): Status ('pending', 'running', 'completed', 'failed')
        details (dict): Optional additional details
    """
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config_data = json.load(f)
        
        # Initialize completion_status if it doesn't exist
        if 'completion_status' not in config_data:
            config_data['completion_status'] = {}
        
        # Update the status
        status_entry = {
            'status': status,
            'timestamp': os.popen('date').read().strip()
        }
        
        if details:
            status_entry.update(details)
        
        config_data['completion_status'][step] = status_entry
        
        # Write back to file
        with open(config_path, 'w', encoding='utf-8') as f:
            json.dump(config_data, f, indent=2)
            
        print(f"Updated {step} status to '{status}' in {config_path}")
        
    except Exception as e:
        print(f"WARNING: Could not update completion status in {config_path}: {e}")

def check_and_create_backfill_spec(config_data, expected_total, actual_total, events_per_file):
    """
    Check if backfill jobs are needed and create specification
    
    Args:
        config_data (dict): Configuration data
        expected_total (int): Expected total events
        actual_total (int): Actual merged events
        events_per_file (int): Events per job file
        
    Returns:
        dict: Backfill specification, or None if no backfill needed
    """
    if expected_total is None:
        print("  NOTE: No expected_total_events specified, skipping backfill check")
        return None
    
    if actual_total >= expected_total:
        print(f"  Event count OK: {actual_total} >= {expected_total} (expected)")
        return None
    
    missing_events = expected_total - actual_total
    backfill_jobs = (missing_events + events_per_file - 1) // events_per_file
    
    # Calculate seed offset to ensure different random sequence
    # Use a large offset (100000) to avoid overlap with original job seeds
    original_jobs = extract_config_value(config_data,
        ['rbhg_config.event_counts.total_jobs', 'total_jobs', 
         'event_counts.total_jobs'], 0)
    
    seed_offset = max(100000, original_jobs + 10000)
    
    backfill_spec = {
        'missing_events': missing_events,
        'backfill_jobs': backfill_jobs,
        'events_per_job': events_per_file,
        'seed_offset': seed_offset,
        'workflow_suffix': '_backfill',
        'reason': f'Original generation completed with {actual_total}/{expected_total} events',
        'recommended_action': f'Submit {backfill_jobs} additional jobs with seed offset {seed_offset}'
    }
    
    return backfill_spec

def process_individual_directory(config_data, config_path):
    """
    Process a single directory (Level 1 config)
    
    Args:
        config_data (dict): Individual directory configuration
        config_path (str): Path to the configuration file
        
    Returns:
        dict: Processing results
    """
    print(f"\nProcessing individual directory: {extract_config_value(config_data, ['rbhg_config.generation_info.study_name', 'rbhg_config.generation_info.nametag', 'directory_name'], 'Unknown')}")
    
    results = {
        'success': False,
        'merged_hddm_file': None,
        'event_count': 0,
        'mcwrapper_created': False
    }
    
    # Update status to running
    update_completion_status(config_path, 'preparation', 'running')
    
    try:
        # Extract paths from config - flexible for different JSON structures
        vectors_dir = extract_config_value(config_data, 
            ['rbhg_config.directory_paths.rbhg_generation.vectors_hddm_directory',
             'vectors_hddm_directory', 'directory_paths.rbhg_generation.vectors_hddm_directory', 
             'vectorspath', 'vectors_directory'])
        
        directory_name = extract_config_value(config_data,
            ['rbhg_config.generation_info.study_name', 'rbhg_config.generation_info.nametag',
             'directory_name', 'generation_info.study_name', 'study_name', 
             'workflow_name', 'name'])
        
        if not vectors_dir:
            raise ValueError("Could not find vectors directory path in configuration")
        if not directory_name:
            raise ValueError("Could not find directory/study name in configuration")
        
        # Expected events per file - flexible field names
        expected_events_per_file = extract_config_value(config_data,
            ['rbhg_config.event_counts.events_per_file', 'rbhg_config.downstream_configs.mc_simulation.events_per_mc_file',
             'nevents_perfile', 'event_counts.events_per_file', 'events_per_file'], 
            default=50000)
        
        print(f"  Vectors directory: {vectors_dir}")
        print(f"  Expected events per file: {expected_events_per_file}")
        
        # Step 1: Merge HDDM files and validate event counts
        # Use the generation directory name as the base for the merged HDDM filename
        # e.g. 1801_0DEG_DBLRAD_GlueX_ee_Berlin_v113.hddm for clarity
        gen_dir_name = os.path.basename(os.path.dirname(vectors_dir))
        merged_hddm_file = os.path.join(vectors_dir, f"{gen_dir_name}.hddm")

        print(f"  Merging HDDM files to: {merged_hddm_file}")
        
        merge_results = merge_and_validate_hddm(
            input_dir=vectors_dir,
            output_file=merged_hddm_file,
            expected_events_per_file=expected_events_per_file,
            pattern="vectors_*.hddm",
            validate_before=True,
            validate_after=True
        )
        
        if not merge_results['merge_successful']:
            raise RuntimeError("HDDM merge failed")
        
        if not merge_results['validation_passed']:
            print("  WARNING: Event count validation failed, but continuing...")
        
        results['merged_hddm_file'] = merged_hddm_file
        results['event_count'] = merge_results['event_counts'].get('output_actual', 0)
        
        # Step 2: Check if backfill jobs are needed
        expected_total_events = extract_config_value(config_data,
            ['rbhg_config.event_counts.total_events', 'total_events',
             'event_counts.total_events'], None)
        
        backfill_info = check_and_create_backfill_spec(
            config_data, 
            expected_total_events, 
            results['event_count'],
            expected_events_per_file
        )
        
        if backfill_info:
            results['backfill_needed'] = True
            results['backfill_spec'] = backfill_info
            print(f"  BACKFILL REQUIRED: {backfill_info['missing_events']} events missing")
            print(f"    Will need {backfill_info['backfill_jobs']} additional jobs")
        else:
            results['backfill_needed'] = False
        
        # Step 3: Create MCWrapper configuration directory
        # This goes on /work/ next to vectors/, not on /volatile/
        mcwrapper_dir = os.path.join(os.path.dirname(vectors_dir), 'MCWrapper')
        
        if not os.path.exists(mcwrapper_dir):
            os.makedirs(mcwrapper_dir, exist_ok=True)
            print(f"  Created MCWrapper directory: {mcwrapper_dir}")
        
        # Create MCWrapper configuration files (reuse logic from original prepSim.py)
        # Pass the actual merged event count to use in gluex_MC_sub scripts
        mcwrapper_success = create_mcwrapper_config(config_data, mcwrapper_dir, merged_hddm_file, results['event_count'])
        results['mcwrapper_created'] = mcwrapper_success
        
        if mcwrapper_success:
            results['success'] = True
            
            # Update completion status
            completion_details = {
                'merged_hddm_file': merged_hddm_file,
                'event_count': results['event_count'],
                'expected_events': expected_total_events,
                'mcwrapper_directory': mcwrapper_dir,
                'validation_passed': merge_results['validation_passed']
            }
            
            # Add backfill info if needed
            if results.get('backfill_needed'):
                completion_details['backfill_required'] = True
                completion_details['backfill_spec'] = results['backfill_spec']
            
            update_completion_status(config_path, 'preparation', 'completed', completion_details)
            print(f"  SUCCESS: Individual directory processing completed")
        else:
            raise RuntimeError("MCWrapper configuration creation failed")
            
    except Exception as e:
        print(f"  ERROR: Individual directory processing failed: {e}")
        update_completion_status(config_path, 'preparation', 'failed', {'error': str(e)})
    
    return results

def create_mcwrapper_config(config_data, mcwrapper_dir, merged_hddm_file, actual_event_count):
    """
    Create MCWrapper configuration files using templates and preset configs
    
    Args:
        config_data (dict): Configuration data
        mcwrapper_dir (str): MCWrapper directory path
        merged_hddm_file (str): Path to merged HDDM file
        actual_event_count (int): Actual number of events in merged HDDM file
        
    Returns:
        bool: True if successful
    """
    import subprocess
    import platform
    import shutil
    import glob as glob_module
    import json
    
    try:
        template_dir = extract_config_value(config_data,
            ['rbhg_config.directory_paths.templates.mcwrapper_template_directory',
             'template_directory'], 
            os.path.join(os.path.dirname(__file__), 'template_MCWrapper'))
        
        # Extract information from config - flexible field access  
        run_period = extract_config_value(config_data,
            ['rbhg_config.physics_settings.run_period', 'run_period', 
             'physics_settings.run_period'], '1801')
        
        run_number = extract_config_value(config_data,
            ['rbhg_config.physics_settings.run_number', 'run_number',
             'physics_settings.run_number'], '42559')
        
        polarization_raw = extract_config_value(config_data,
            ['rbhg_config.physics_settings.polarization_deg', 'polarization', 
             'physics_settings.polarization_deg', 'polarization_deg'], '45')
        # Only append 'DEG' if not already present (e.g., AMO should stay as AMO, not AMODEG)
        polarization = polarization_raw if polarization_raw.upper() == 'AMO' or polarization_raw.endswith('DEG') else polarization_raw + 'DEG'
        
        particle_type = extract_config_value(config_data,
            ['rbhg_config.physics_settings.lepton_type', 'particle_type', 
             'physics_settings.lepton_type', 'lepton_type'], 'ee')
        
        experiment_type = extract_config_value(config_data,
            ['rbhg_config.physics_settings.experiment', 'experiment_type', 
             'physics_settings.experiment', 'experiment'], 'GlueX')
        
        target = extract_config_value(config_data,
            ['rbhg_config.physics_settings.target', 'rbhg_config.downstream_configs.mc_simulation.target',
             'target'], 'p')
        
        username = extract_config_value(config_data,
            ['rbhg_config.user_settings.username', 'username', 'user_settings.username'],
            os.environ.get('USER', 'acschick'))
        
        # Use the actual event count from the merged HDDM file (passed as parameter)
        # This is the true count after filtering empty files, not the expected count
        nevents = actual_event_count
        print(f"    Using actual merged event count: {nevents}")
        
        # Extract MCWrapper settings from JSON (these were added manually or by swif2_RBHG.py)
        mcw_settings = extract_config_value(config_data,
            ['rbhg_config.mcwrapper_settings', 'mcwrapper_settings'], {})
        
        # Get settings from JSON or use defaults
        env_xml = mcw_settings.get('recon_env', '/group/halld/www/halldweb/html/halld_versions/version.xml')
        ana_env_xml = mcw_settings.get('analysis_env', '')
        # Fallback: if environment XML not provided, use RunPeriods.json entry
        if not env_xml:
            try:
                framework_home = extract_config_value(config_data,
                    ['rbhg_config.directory_paths.base_paths.framework_home', 'framework_home'], os.path.dirname(__file__))
                runperiods_path = os.path.join(framework_home, 'RunPeriods.json')
                with open(runperiods_path, 'r') as rp_f:
                    rp_data = json.load(rp_f)
                rp_entry = rp_data.get(str(run_period)) or rp_data.get('default') or {}
                env_xml = rp_entry.get('recon_env', env_xml)
                ana_env_xml = rp_entry.get('analysis_env', ana_env_xml)
                if env_xml:
                    print(f"    Fallback: using recon_env from RunPeriods.json: {env_xml}")
            except Exception as e:
                print(f"    WARNING: Could not read recon_env from RunPeriods.json: {e}")
        bkg = mcw_settings.get('background_version', 'NONE')
        rcdb_query = mcw_settings.get('rcdb_query', '')
        
        # Substitute {polDeg} placeholder in RCDB query with actual polarization value
        # For AMO (amorphous), use '-1.0'; for degree values like '0DEG', extract the number
        if '{polDeg}' in rcdb_query:
            if polarization.upper() == 'AMO':
                pol_value = '-1.0'
            else:
                # Extract numeric value from polarization (e.g., '0DEG' -> '0', '45DEG' -> '45')
                pol_value = polarization.replace('DEG', '').replace('deg', '')
            rcdb_query = rcdb_query.replace('{polDeg}', pol_value)
        
        variation = mcw_settings.get('variation', 'mc')
        batch_system = mcw_settings.get('batch_system', 'swif2')
        run_range = mcw_settings.get('run_range', f'{run_number} {run_number}')

        # Honor USE_CHARACTERISTIC_RUN if requested. This can be provided in
        # the per-generation mcwrapper_settings as 'USE_CHARACTERISTIC_RUN' (or
        # lower-case 'use_characteristic_run'), or as a top-level flag in the
        # generation JSON under rbhg_config (rbhg_config.mcwrapper_use_characteristic_run).
        use_characteristic = mcw_settings.get('USE_CHARACTERISTIC_RUN', mcw_settings.get('use_characteristic_run', False))
        if not use_characteristic:
            # Check top-level config fields as a fallback
            use_characteristic = extract_config_value(config_data,
                ['rbhg_config.mcwrapper_use_characteristic_run', 'rbhg_config.use_characteristic_run'], False)

        if use_characteristic:
            # Attempt to read RunPeriods.json from framework home
            framework_home = extract_config_value(config_data,
                ['rbhg_config.directory_paths.base_paths.framework_home', 'framework_home'], os.path.dirname(__file__))
            runperiods_path = os.path.join(framework_home, 'RunPeriods.json')
            try:
                with open(runperiods_path, 'r') as rp_f:
                    rp_data = json.load(rp_f)

                rp_entry = rp_data.get(str(run_period)) or rp_data.get('default') or {}
                pol_entry = rp_entry.get('Polarizations', {}).get(polarization, {})
                char_run = pol_entry.get('characteristic_run') or rp_entry.get('characteristic_run')

                if char_run:
                    # Use single characteristic run as the run range (single value)
                    run_range = str(char_run)
                    print(f"    USE_CHARACTERISTIC_RUN enabled: setting run_range to characteristic run {run_range}")
                else:
                    print(f"    WARNING: USE_CHARACTERISTIC_RUN enabled but no characteristic_run found for run_period={run_period}, polarization={polarization}; falling back to configured run_range={run_range}")
            except Exception as e:
                print(f"    WARNING: failed to read RunPeriods.json ({runperiods_path}): {e}; using configured run_range={run_range}")
        
        # Pipeline control flags
        geant = mcw_settings.get('geant', 1)
        mcsmear = mcw_settings.get('mcsmear', 1)
        recon = mcw_settings.get('recon', 1)
        clean_geant = mcw_settings.get('clean_geant', 1)
        clean_mcsmear = mcw_settings.get('clean_mcsmear', 1)
        clean_recon = mcw_settings.get('clean_recon', 0)
        batch_mode = mcw_settings.get('batch_mode', 2)
        
        # Check for 10x flag - if set, use 100k events per file instead of 500k
        # This prevents jobs from running too long
        tenx_flag = extract_config_value(config_data,
            ['rbhg_config.event_counts.tenx_flag', 'event_counts.tenx_flag', 'tenx_flag'], False)
        default_per_file = 100000 if tenx_flag else 500000
        events_per_file = mcw_settings.get('events_per_file', default_per_file)
        
        jana_config = mcw_settings.get('jana_config', '')
        # Fallback: if jana_config not provided in mcwrapper_settings, try to
        # auto-select from RunPeriods.json using run_period, polarization, and particle_type
        if not jana_config:
            try:
                framework_home = extract_config_value(config_data,
                    ['rbhg_config.directory_paths.base_paths.framework_home', 'framework_home'], os.path.dirname(__file__))
                runperiods_path = os.path.join(framework_home, 'RunPeriods.json')
                with open(runperiods_path, 'r') as rp_f:
                    rp_data = json.load(rp_f)

                rp_entry = rp_data.get(str(run_period)) or rp_data.get('default') or {}

                # Build prioritized candidate keys based on particle type
                p = (particle_type or '').lower()
                candidates = []
                if p in ['ee', 'epem', 'epemissp', 'epimu', 'epimu', 'ep'] or 'ep' in p:
                    candidates = ['jana_epem_ReaFil_config', 'jana_epem_ml_skim_config', 'jana_epem_config']
                elif p in ['mupmum', 'mumu', 'mu'] or 'mu' in p:
                    candidates = ['jana_mupmum_ReaFil_config', 'jana_mupmum_config']
                else:
                    candidates = ['jana_epem_config', 'jana_mupmum_config', 'jana_epimu_config']

                for key in candidates:
                    val = rp_entry.get(key, '') or rp_entry.get(key.replace('_ReaFil', ''), '')
                    if val:
                        jana_config = val
                        print(f"    Fallback: selected JANA config from RunPeriods.json key '{key}': {jana_config}")
                        break

                # If selected jana_config is relative (e.g., 'standard.config'), try to resolve
                if jana_config and not os.path.isabs(jana_config):
                    # Look in the mcwrapper template jana_configs directory if available
                    try:
                        possible = os.path.join(template_dir, 'jana_configs', jana_config)
                        if os.path.exists(possible):
                            jana_config = possible
                            print(f"    Resolved relative jana_config to: {jana_config}")
                    except Exception:
                        pass

            except Exception as e:
                print(f"    WARNING: Could not determine jana_config from RunPeriods.json: {e}")
        
        # Get the generation workflow name and directory structure to maintain consistency
        # The directory name from generation looks like: 1801_0DEG_DBLRAD_GlueX_ee_Berlin_v113
        gen_output_dir = extract_config_value(config_data,
            ['rbhg_config.directory_paths.rbhg_generation.output_directory',
             'directory_paths.rbhg_generation.output_directory'], '')
        
        # Get the last directory component (e.g., 1801_0DEG_DBLRAD_GlueX_ee_Berlin_v113)
        gen_dir_name = os.path.basename(gen_output_dir) if gen_output_dir else ''
        
        # Get generation workflow name and replace 'gen_' with 'SIM_'
        gen_workflow_name = extract_config_value(config_data,
            ['rbhg_config.directory_paths.swif2_workflow_names.generation',
             'swif2_workflow_names.generation',
             'directory_paths.swif2_workflow_names.generation'], '')
        
        if gen_workflow_name and gen_workflow_name.startswith('gen_'):
            workflow_name = 'SIM_' + gen_workflow_name[4:]  # Replace 'gen_' with 'SIM_'
        else:
            # Fallback construction if workflow name not found
            workflow_name = f"SIM_{gen_dir_name}" if gen_dir_name else "SIM_Unknown"
        
        # Get study name and nametag for directory structure
        study_name = extract_config_value(config_data,
            ['rbhg_config.generation_info.study_name', 'study_name',
             'generation_info.study_name'], 'Study')
        
        nametag = extract_config_value(config_data,
            ['rbhg_config.generation_info.nametag', 'nametag',
             'generation_info.nametag'], '')
        
        # Construct output directory path
        # Format: /volatile/halld/home/{user}/RBHG/{study}/{nametag}/{gen_dir_name}/Simulation
        # Example: /volatile/halld/home/acschick/RBHG/CobremTest/A3_FFN/1801_0DEG_DBLRAD_GlueX_ee_Berlin_v113/Simulation
        if nametag:
            output_base_dir = f"/volatile/halld/home/{username}/RBHG/{study_name}/{nametag}/{gen_dir_name}/Simulation"
        else:
            output_base_dir = f"/volatile/halld/home/{username}/RBHG/{study_name}/{gen_dir_name}/Simulation"
        
        # Get HDDM file size for disk request
        hddm_size_bytes = os.path.getsize(merged_hddm_file)
        # Multiply by 2.2 and round up to nearest GB
        hddm_size_gb = math.ceil((hddm_size_bytes / (1024**3)) * 2.2)
        disk_request_gb = hddm_size_gb
        disk_request = f"{disk_request_gb}GB"
        
        print(f"    Creating MCWrapper config for {run_period} {polarization} {particle_type}")
        print(f"    Using template directory: {template_dir}")
        print(f"    HDDM file: {merged_hddm_file}")
        print(f"    Target: {target}, Events: {nevents}")
        print(f"    Disk request: {disk_request}")
        
        # Helper function to run sed commands
        def run_sed_replacement(file_path, replacements_dict):
            """Run sed replacements on a file"""
            # Detect OS for sed syntax
            os_type = platform.system()
            sed_option = '-i' if os_type == "Linux" else '-i ""'
            
            for placeholder, replacement in replacements_dict.items():
                # Escape special characters for sed
                escaped_replacement = replacement.replace('/', '\\/')
                cmd = f"sed {sed_option} 's|{placeholder}|{escaped_replacement}|g' {file_path}"
                try:
                    subprocess.run(cmd, shell=True, check=True, capture_output=True)
                except subprocess.CalledProcessError as e:
                    print(f"      WARNING: sed replacement failed for {placeholder}: {e}")
            
        success_count = 0
        
        # Determine what to use for ENVIRONMENT_FILE in MC.config
        # For swif2cont (containers): Must use XML directly, can't use custom .csh
        # For CPP: Need custom .csh for geometry settings
        if batch_system == 'swif2cont' and experiment_type != 'CPP':
            # Use XML directly for containerized old software (1801, 1808)
            environment_file_path = env_xml
            print(f"    Using XML environment directly (swif2cont): {env_xml}")
        else:
            # Create and use MCEnvironment.csh for CPP or native swif2
            # 1. Create MCEnvironment.csh from template
            print(f"    Creating MCEnvironment.csh...")
            env_csh_path = os.path.join(mcwrapper_dir, 'MCEnvironment.csh')
            template_env = os.path.join(template_dir, 'temp_MCEnvironment.csh')
            
            try:
                shutil.copy(template_env, env_csh_path)
                
                env_replacements = {
                    'RBHGCUEUSERNAME': username,
                    'RBHGPATHTOXML': env_xml
                }
                
                # CPP-specific settings
                if experiment_type == 'CPP':
                    env_replacements['RBHGCPPGEOMETRY'] = 'setenv JANA_GEOMETRY_URL "ccdb:///GEOMETRY/cpp_HDDS.xml"'
                    env_replacements['RBHGCPPMCVARIATION'] = 'setenv JANA_CALIB_CONTEXT "variation=mc_cpp"'
                else:
                    env_replacements['RBHGCPPGEOMETRY'] = ''
                    env_replacements['RBHGCPPMCVARIATION'] = ''
                
                run_sed_replacement(env_csh_path, env_replacements)
                print(f"      Created: {env_csh_path}")
                environment_file_path = env_csh_path
                success_count += 1
            except Exception as e:
                print(f"      ERROR creating MCEnvironment.csh: {e}")
                environment_file_path = env_xml  # Fallback to XML
        
        # 2. Copy JANA config from mcwrapper_settings
        print(f"    Creating jana config...")
        jana_config_dest = os.path.join(mcwrapper_dir, os.path.basename(jana_config))
        
        try:
            if jana_config and os.path.exists(jana_config):
                shutil.copy(jana_config, jana_config_dest)
                print(f"      Copied: {jana_config} -> {jana_config_dest}")
                success_count += 1
            else:
                print(f"      ERROR: JANA config not found: {jana_config}")
        except Exception as e:
            print(f"      ERROR copying JANA config: {e}")
        
        # 3. Create MC.config from template
        print(f"    Creating MC.config...")
        mc_config_path = os.path.join(mcwrapper_dir, 'MC.config')
        template_mc = os.path.join(template_dir, 'temp_MC.config')
        
        try:
            shutil.copy(template_mc, mc_config_path)
            
            # Determine container OS settings based on batch system
            if batch_system == 'swif2cont':
                container_os_lines = """GENERATOR_OS=LOCAL
POSTGEN_OS=LOCAL
SIMULATION_OS=CENTOS7
MCSMEAR_OS=CENTOS7
RECON_OS=CENTOS7
ANA_OS=CENTOS7"""
            else:
                container_os_lines = ""
            
            mc_replacements = {
                'RBHGWORKFLOWNAME': workflow_name,
                'RBHGOUTPUTBASEDIR': output_base_dir,
                'RBHGSPATHTOHDDMFILE': merged_hddm_file,
                'RBHGINCLUDEBACKGROUND': bkg,
                'RBHG_RECON_ENVIRONMENT': environment_file_path,
                'RBHG_ANALYSIS_ENVIRONMENT': f'ANA_ENVIRONMENT_FILE={ana_env_xml}' if ana_env_xml else '',
                'RBHG_RCDBQUERY': rcdb_query,
                'RBHG_JANACONFIGFILE': jana_config_dest,
                'RBHG_EXPERIMENT': experiment_type,
                'RBHGDISKREQUEST': disk_request,
                'RBHGVARIATIONMC': variation,
                'RBHGBATCHSYSTEM': batch_system,
                'RBHGGENERATOROSSETTING': container_os_lines.split('\n')[0] if container_os_lines else '',
                'RBHGPOSTGENOSSETTING': container_os_lines.split('\n')[1] if container_os_lines else '',
                'RBHGSIMULATIONOSSETTING': container_os_lines.split('\n')[2] if container_os_lines else '',
                'RBHGMCSMEAROSSETTING': container_os_lines.split('\n')[3] if container_os_lines else '',
                'RBHGRECONOSSETTING': container_os_lines.split('\n')[4] if container_os_lines else '',
                'RBHGANAOSSETTING': container_os_lines.split('\n')[5] if container_os_lines else ''
            }
            
            run_sed_replacement(mc_config_path, mc_replacements)
            print(f"      Created: {mc_config_path}")
            success_count += 1
        except Exception as e:
            print(f"      ERROR creating MC.config: {e}")
        
        # 4. Create gluex_MC_sub.csh from template  
        print(f"    Creating gluex_MC_sub.csh...")
        sub_csh_path = os.path.join(mcwrapper_dir, 'gluex_MC_sub.csh')
        template_sub = os.path.join(template_dir, 'temp_gluex_MC_sub.csh')
        
        try:
            shutil.copy(template_sub, sub_csh_path)
            
            # Format run range (handle various formats: "40856 42559", "40856 - 42559", "40856-42559", or "42559")
            # Clean up any existing dashes and spaces, then format consistently
            run_range_cleaned = run_range.replace('-', ' ').strip()
            run_range_parts = run_range_cleaned.split()
            if len(run_range_parts) == 2:
                run_range_formatted = f"{run_range_parts[0]}-{run_range_parts[1]}"
            elif len(run_range_parts) == 1:
                run_range_formatted = run_range_parts[0]
            else:
                # Fallback: just use as-is
                run_range_formatted = run_range.replace(' ', '-')

            # Determine farm_out base for logs
            farm_out_base = extract_config_value(config_data,
                ['rbhg_config.directory_paths.base_paths.farm_out_base',
                 'rbhg_config.directory_paths.base_paths.farm_out',
                 'rbhg_config.directory_paths.farm_out_base',
                 'rbhg_config.directory_paths.base_paths.farm_out_base'],
                None)
            if not farm_out_base:
                # Fallback to constructed path using generator_type from JSON
                generator_type = extract_config_value(config_data, ['rbhg_config.directory_paths.base_paths.generator_type'], 'RBHG')
                farm_out_base = f"/farm_out/{username}/{generator_type}/{study_name}"

            # Log dir should point to simulation logs under farm_out
            logdir = os.path.join(farm_out_base, 'simulation')
            # Ensure trailing slash for readability in configs
            if not logdir.endswith('/'):
                logdir = logdir + '/'
            
            # Simplified replacements - minimal working format for containerized MCWrapper
            sub_replacements = {
                'RBHGPATHTOMCCONFIG': mc_config_path,
                'RBHGRUNRANGE': run_range_formatted,
                'RBHGNEVENTS': str(nevents),
                'RBHGPERFILE': str(events_per_file),
                'RBHGBATCHMODE': str(batch_mode),
                'RBHGLOGDIR': logdir
            }
            
            run_sed_replacement(sub_csh_path, sub_replacements)
            
            # Make script executable
            os.chmod(sub_csh_path, 0o755)
            
            print(f"      Created: {sub_csh_path}")
            success_count += 1
        except Exception as e:
            print(f"      ERROR creating submission script: {e}")
        
        # Determine expected file count
        # 3 files if using XML directly (no MCEnvironment.csh): MC.config, jana_config, gluex_MC_sub.csh
        # 4 files if creating MCEnvironment.csh: MC.config, MCEnvironment.csh, jana_config, gluex_MC_sub.csh
        expected_count = 3 if (batch_system == 'swif2cont' and experiment_type != 'CPP') else 4
        
        print(f"    Created {success_count}/{expected_count} MCWrapper configuration files")
        return success_count == expected_count
        
    except Exception as e:
        print(f"    ERROR creating MCWrapper config: {e}")
        import traceback
        traceback.print_exc()
        return False

def process_run_period_master(config_data, config_path):
    """
    Process run period master (Level 2 config)
    
    Args:
        config_data (dict): Run period master configuration
        config_path (str): Path to the configuration file
        
    Returns:
        dict: Processing results
    """
    print(f"\nProcessing run period master configuration")
    print(f"  Dataset: {config_data.get('dataset_name', 'Unknown')}")
    print(f"  Run period: {config_data.get('run_period', 'Unknown')}")
    
    results = {
        'success': True,
        'processed_configs': [],
        'failed_configs': []
    }
    
    # Update master status
    update_completion_status(config_path, 'preparation', 'running')
    
    # Process each individual configuration
    individual_config_patterns = extract_config_value(config_data, 
        ['run_period_master_config.batch_processing.individual_configs',
         'batch_processing.individual_configs', 'individual_configs'], [])
    
    if not individual_config_patterns:
        print("  WARNING: No individual_configs found in run period master")
        return results
    
    print(f"  Processing {len(individual_config_patterns)} individual configuration patterns")
    
    # Convert patterns to actual file paths
    actual_config_paths = []
    config_base_dir = os.path.dirname(config_path)
    
    for pattern in individual_config_patterns:
        print(f"    Searching for pattern: {pattern}")
        
        # Handle relative patterns by making them relative to the config directory
        if not os.path.isabs(pattern):
            search_pattern = os.path.join(config_base_dir, pattern)
        else:
            search_pattern = pattern
        
        # Use glob to find matching files
        import glob
        matching_files = glob.glob(search_pattern)
        
        if matching_files:
            actual_config_paths.extend(matching_files)
            print(f"      Found {len(matching_files)} matching files")
        else:
            print(f"      WARNING: No files found for pattern {pattern}")
    
    if not actual_config_paths:
        print("  ERROR: No individual configuration files found")
        results['success'] = False
        return results
    
    print(f"  Total individual configs to process: {len(actual_config_paths)}")
    
    for i, individual_config_path in enumerate(actual_config_paths, 1):
        print(f"\n  [{i}/{len(actual_config_paths)}] Processing: {os.path.basename(individual_config_path)}")
        
        try:
            # Load and process individual config
            individual_data, _, _ = load_and_validate_config(individual_config_path)
            individual_results = process_individual_directory(individual_data, individual_config_path)
            
            if individual_results['success']:
                # Extract directory name more safely
                directory_name = extract_config_value(individual_data,
                    ['rbhg_config.generation_info.study_name', 'rbhg_config.generation_info.nametag',
                     'directory_name', 'generation_info.study_name', 'study_name']) or os.path.basename(individual_config_path)
                
                results['processed_configs'].append({
                    'config_path': individual_config_path,
                    'directory': directory_name,
                    'event_count': individual_results['event_count'],
                    'merged_hddm': individual_results['merged_hddm_file']
                })
            else:
                results['failed_configs'].append(individual_config_path)
                results['success'] = False
                
        except Exception as e:
            print(f"    ERROR processing {individual_config_path}: {e}")
            results['failed_configs'].append(individual_config_path)
            results['success'] = False
    
    # Update master completion status
    completion_details = {
        'total_configs': len(actual_config_paths),
        'successful_configs': len(results['processed_configs']),
        'failed_configs': len(results['failed_configs']),
        'processed_directories': [r['directory'] for r in results['processed_configs']],
        'total_events': sum(r['event_count'] for r in results['processed_configs'])
    }
    
    if results['success']:
        update_completion_status(config_path, 'preparation', 'completed', completion_details)
        print(f"\n  SUCCESS: Run period master processing completed")
        print(f"    Processed: {len(results['processed_configs'])} configurations")
        print(f"    Total events: {completion_details['total_events']}")
    else:
        update_completion_status(config_path, 'preparation', 'failed', completion_details)
        print(f"\n  PARTIAL SUCCESS: {len(results['processed_configs'])}/{len(actual_config_paths)} completed")
    
    return results

def process_ffs_master(config_data, config_path):
    """
    Process FFS master (Level 3 config)
    
    Args:
        config_data (dict): FFS master configuration
        config_path (str): Path to the configuration file
        
    Returns:
        dict: Processing results
    """
    print(f"\nProcessing Form Factor Study (FFS) master configuration")
    
    # Handle nested structure: datasets might be in ffs_master_config or at top level
    if 'ffs_master_config' in config_data:
        ffs_config = config_data['ffs_master_config']
    else:
        ffs_config = config_data
    
    study_name = ffs_config.get('study_info', {}).get('study_name', 'Unknown')
    print(f"  Study: {study_name}")
    
    results = {
        'success': True,
        'processed_datasets': [],
        'failed_datasets': []
    }
    
    # Update FFS master status
    update_completion_status(config_path, 'preparation', 'running')
    
    # Process each dataset
    datasets = ffs_config.get('datasets', {})
    print(f"  Processing {len(datasets)} datasets")
    
    for dataset_name, dataset_info in datasets.items():
        print(f"\n  Processing dataset: {dataset_name}")
        
        try:
            # Get the config path - could be 'config_file' or 'run_period_master_config'
            config_path_rel = dataset_info.get('config_file') or dataset_info.get('run_period_master_config')
            if not config_path_rel:
                raise ValueError(f"No config_file or run_period_master_config found for dataset {dataset_name}")
            
            # Make path absolute relative to FFS master config directory
            ffs_dir = os.path.dirname(os.path.abspath(config_path))
            dataset_config_path = os.path.join(ffs_dir, config_path_rel)
            
            if not os.path.exists(dataset_config_path):
                raise FileNotFoundError(f"Config file not found: {dataset_config_path}")
            
            # Load and determine what type of config this is
            dataset_config_data, config_type, config_level = load_and_validate_config(dataset_config_path)
            
            # Process based on config level
            if config_level == 1:
                # Individual config - process directly
                dataset_results = process_individual_directory(dataset_config_data, dataset_config_path)
                if dataset_results['success']:
                    results['processed_datasets'].append({
                        'dataset_name': dataset_name,
                        'config_path': dataset_config_path,
                        'configurations': 1,
                        'total_events': dataset_results.get('event_count', 0)
                    })
                else:
                    results['failed_datasets'].append(dataset_name)
                    results['success'] = False
            elif config_level == 2:
                # Run period master - process as master
                dataset_results = process_run_period_master(dataset_config_data, dataset_config_path)
                if dataset_results['success']:
                    results['processed_datasets'].append({
                        'dataset_name': dataset_name,
                        'config_path': dataset_config_path,
                        'configurations': len(dataset_results['processed_configs']),
                        'total_events': sum(r['event_count'] for r in dataset_results['processed_configs'])
                    })
                else:
                    results['failed_datasets'].append(dataset_name)
                    results['success'] = False
            else:
                raise ValueError(f"Unexpected config level {config_level} for dataset {dataset_name}")
                
        except Exception as e:
            print(f"    ERROR processing dataset {dataset_name}: {e}")
            results['failed_datasets'].append(dataset_name)
            results['success'] = False
    
    # Update FFS master completion status
    completion_details = {
        'total_datasets': len(datasets),
        'successful_datasets': len(results['processed_datasets']),
        'failed_datasets': len(results['failed_datasets']),
        'dataset_summary': results['processed_datasets']
    }
    
    if results['success']:
        update_completion_status(config_path, 'preparation', 'completed', completion_details)
        print(f"\n  SUCCESS: FFS master processing completed")
        print(f"    Datasets processed: {len(results['processed_datasets'])}")
        
        total_events = sum(d['total_events'] for d in results['processed_datasets'])
        total_configs = sum(d['configurations'] for d in results['processed_datasets'])
        print(f"    Total configurations: {total_configs}")
        print(f"    Total events: {total_events}")
    else:
        update_completion_status(config_path, 'preparation', 'failed', completion_details)
        print(f"\n  PARTIAL SUCCESS: {len(results['processed_datasets'])}/{len(datasets)} datasets completed")
    
    return results

def main():
    """Main function to handle command line arguments and route processing"""
    
    parser = argparse.ArgumentParser(
        description="JSON-driven simulation preparation tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process individual directory (Level 1)
  python prepareSimulation.py path/to/rbhg_config.json
  
  # Process run period master (Level 2)  
  python prepareSimulation.py path/to/run_period_master_config.json
  
  # Process FFS master (Level 3)
  python prepareSimulation.py path/to/ffs_master_config.json
  
  # Dry run (validate configs without processing)
  python prepareSimulation.py path/to/config.json --dry-run
        """
    )
    
    parser.add_argument('config',
                       help='Path to JSON configuration file (Level 1, 2, or 3)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Validate configuration without processing')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Load and validate configuration
    try:
        config_data, config_type, config_level = load_and_validate_config(args.config)
        
        print(f"Loaded configuration: {args.config}")
        print(f"Configuration type: {config_type} (Level {config_level})")
        
        if args.dry_run:
            print("DRY RUN MODE - Configuration validated successfully")
            return 0
        
        # Route to appropriate processor
        if config_level == 1:
            results = process_individual_directory(config_data, args.config)
        elif config_level == 2:
            results = process_run_period_master(config_data, args.config)
        elif config_level == 3:
            results = process_ffs_master(config_data, args.config)
        else:
            raise ValueError(f"Unknown configuration level: {config_level}")
        
        # Report final results
        if results['success']:
            print(f"\nSUCCESS: Processing completed for {config_type} configuration")
            return 0
        else:
            print(f"\nFAILED: Processing failed for {config_type} configuration")
            return 1
            
    except Exception as e:
        print(f"ERROR: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
