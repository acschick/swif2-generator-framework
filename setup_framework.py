#!/usr/bin/env python3
"""
Framework Setup Script
======================
One-time configuration script that auto-detects user environment and updates
all configuration files with correct paths, username, and shell preferences.

Run this once after cloning the repository:
    ./setup_framework.py

This will:
1. Auto-detect your username
2. Auto-detect framework installation directory
3. Detect your shell type (bash/csh)
4. Generate user-local real-data DSelector defaults
5. Update all .config files in GeneratorConfigExamples/
6. Save settings for downstream scripts (prepareSimulation.py, etc.)
"""

import os
import sys
import json
import glob
from pathlib import Path


def detect_username():
    """Detect current username from environment."""
    username = os.getenv('USER') or os.getenv('USERNAME')
    if not username:
        print("WARNING: Could not auto-detect username")
        username = input("Enter your username: ").strip()
    return username


def detect_framework_home():
    """Detect framework home directory from script location."""
    # Script is in framework root, so just get the directory containing this script
    framework_home = os.path.dirname(os.path.abspath(__file__))
    return framework_home


def detect_shell():
    """Detect user's shell type (bash or csh/tcsh)."""
    shell_path = os.getenv('SHELL', '')
    
    # Check if it's a csh variant
    if 'csh' in shell_path.lower():
        shell_type = 'csh'
    elif 'bash' in shell_path.lower():
        shell_type = 'bash'
    else:
        # Default or unknown - ask user
        print(f"\nDetected shell: {shell_path}")
        print("Could not auto-detect shell type.")
        while True:
            choice = input("Which shell do you use? (bash/csh): ").strip().lower()
            if choice in ['bash', 'csh', 'tcsh']:
                shell_type = 'csh' if choice in ['csh', 'tcsh'] else 'bash'
                break
            print("Please enter 'bash' or 'csh'")
    
    return shell_type


def update_config_file(config_path, username, framework_home):
    """
    Update a single config file with detected username and framework_home.
    Only updates the values, preserves all comments and structure.
    """
    with open(config_path, 'r') as f:
        lines = f.readlines()
    
    updated_lines = []
    for line in lines:
        # Skip comments and empty lines
        if line.strip().startswith('#') or not line.strip():
            updated_lines.append(line)
            continue
        
        # Parse line to get key and value
        parts = line.split()
        if len(parts) >= 1:
            key = parts[0]
            
            # Update USERNAME
            if key == 'USERNAME':
                # Preserve any inline comments
                comment = ''
                if '#' in line:
                    comment = ' ' + line[line.index('#'):]
                else:
                    comment = '\n'
                updated_lines.append(f"USERNAME                     {username}{comment}")
            
            # Update FRAMEWORK_HOME (leave it blank to enable auto-detection)
            elif key == 'FRAMEWORK_HOME':
                # Intentionally left blank - enables portable auto-detection
                # Generators will auto-detect framework root from script location
                updated_lines.append(f"FRAMEWORK_HOME               # (blank = auto-detect) Detected as: {framework_home}\n")
            
            else:
                # Keep line as-is
                updated_lines.append(line)
        else:
            updated_lines.append(line)
    
    # Write updated config
    with open(config_path, 'w') as f:
        f.writelines(updated_lines)
    
    print(f"  ✓ Updated {os.path.basename(config_path)}")


def build_framework_settings(username, framework_home, shell_type):
    """Build the setup settings shared by downstream scripts."""
    volatile_base = f"/volatile/halld/home/{username}"
    farm_out_base = f"/farm_out/{username}"
    return {
        'username': username,
        'framework_home': framework_home,
        'shell_type': shell_type,
        'volatile_base': volatile_base,
        'farm_out_base': farm_out_base,
        'real_data_output_base': f"{volatile_base}/RealDataDSelector/{{RunPeriod}}",
        'real_data_log_base': f"{farm_out_base}/DSelector_logs/RealDataDSelector/{{RunPeriod}}",
        'rbhg_output_base': f"{volatile_base}/RBHG",
        'spizg_output_base': f"{volatile_base}/SPIZG",
        'root_analysis_launch_script': os.path.join(framework_home, 'DSelector_Templates', 'launch.py'),
        'root_analysis_scriptfile': os.path.join(framework_home, 'DSelector_Templates', 'script.sh'),
        'root_analysis_rootscript': os.path.join(framework_home, 'DSelector_Templates', 'Run_Selector.C'),
        'root_analysis_envfile': '/group/halld/www/halldweb/html/halld_versions/version.xml',
        'setup_complete': True
    }


def save_framework_settings(settings):
    """
    Save framework settings to a JSON file for other scripts to use.
    This allows prepareSimulation.py and other downstream scripts to
    automatically use the correct shell type, etc.
    """
    settings_path = os.path.join(settings['framework_home'], 'framework_settings.json')
    with open(settings_path, 'w') as f:
        json.dump(settings, indent=2, fp=f)
    
    return settings_path


def make_directory(path):
    """Create a directory if possible; warn rather than failing setup."""
    try:
        os.makedirs(path, exist_ok=True)
        print(f"  ✓ {path}")
        return True
    except OSError as e:
        print(f"  WARNING: Could not create {path}: {e}")
        return False


def write_default_data_dselector_config(settings, overwrite=False):
    """Create the default real-data DSelector config used by swif2_Data_DSelector.py."""
    config_path = os.path.join(settings['framework_home'], 'data_dselector.config')
    legacy_config_path = os.path.join(settings['framework_home'], 'workflows_root_analysis.config')
    if os.path.exists(config_path) and not overwrite:
        print(f"  ✓ Existing {os.path.basename(config_path)} left unchanged")
        return config_path
    if os.path.exists(legacy_config_path) and not overwrite:
        print("  NOTE: Found legacy workflows_root_analysis.config")
        print("        Writing a fresh data_dselector.config instead of copying old absolute paths.")
    
    content = f"""# Real-data DSelector workflow configuration
# Generated by setup_framework.py
#
# This file is read by swif2_Data_DSelector.py. Values use KEY=VALUE syntax.

RUN_PERIODS=1801,1808
POLARIZATIONS=ALL

TREE_TYPE=epemmissprot__B2_U1
TREE_NAME=epemmissprot__B2_U1_Tree
ANALYSIS_VERSION_1801=ver17
ANALYSIS_VERSION_1808=ver08
ANALYSIS_VERSION_2205=ver01

DSELECTOR_PATH={{FRAMEWORK_HOME}}/DSelector_Templates/2eMissingProton_Systematics/DSelector_2eMissingProton_Systematics.C

OUTPUT_BASE_DIR={settings['real_data_output_base']}
LOG_BASE_DIR={settings['real_data_log_base']}
WORKFLOW_PREFIX=DSELECTOR_DATA
EXISTING_OUTPUT_MODE=fail

TIMELIMIT=8hrs
CORES=1
RAM=3GB
DISK=10GB
SWIF2_ACCOUNT=halld
SWIF2_PARTITION=production
DRY_RUN=False
"""
    with open(config_path, 'w') as f:
        f.write(content)
    print(f"  ✓ Wrote {config_path}")
    return config_path


def update_runperiods_cobrems(framework_home):
    """
    Update RunPeriods.json to point cobrems_distribution entries to the
    framework-local generators/CobremsDistributions/ directory via a
    portable placeholder.
    
    Converts absolute paths like:
      /work/.../CobremsDistributions/cbrem_30327_0deg_rcdb.dat
    To:
      {FRAMEWORK_HOME}/generators/CobremsDistributions/cbrem_30327_0deg_rcdb.dat
    """
    runperiods_path = os.path.join(framework_home, 'RunPeriods.json')
    
    if not os.path.exists(runperiods_path):
        print(f"  WARNING: RunPeriods.json not found at {runperiods_path}")
        return
    
    # Read current RunPeriods.json
    with open(runperiods_path, 'r') as f:
        runperiods = json.load(f)
    
    updated_count = 0
    cobrems_dir = '{FRAMEWORK_HOME}/generators/CobremsDistributions'
    
    # Iterate through all run periods and polarizations
    for period_key, period_data in runperiods.items():
        if not isinstance(period_data, dict):
            continue
        
        polarizations = period_data.get('Polarizations', {})
        for pol_key, pol_data in polarizations.items():
            if 'cobrems_distribution' in pol_data:
                old_path = pol_data['cobrems_distribution']
                
                # Extract just the filename from the old path
                filename = os.path.basename(old_path)
                
                # Construct new framework-local path
                new_path = f"{cobrems_dir}/{filename}"
                
                # Update the entry
                pol_data['cobrems_distribution'] = new_path
                updated_count += 1
    
    # Write updated RunPeriods.json
    with open(runperiods_path, 'w') as f:
        json.dump(runperiods, f, indent=2)
    
    print(f"  ✓ Updated {updated_count} cobrems_distribution paths in RunPeriods.json")
    print(f"    New path format: {cobrems_dir}/<filename>")


def main():
    print("="*70)
    print("SWIF2 Generator Framework - One-Time Setup")
    print("="*70)
    print()
    
    # Step 1: Detect username
    print("[1/5] Detecting username...")
    username = detect_username()
    print(f"  ✓ Username: {username}")
    print()
    
    # Step 2: Detect framework home
    print("[2/5] Detecting framework installation directory...")
    framework_home = detect_framework_home()
    print(f"  ✓ Framework home: {framework_home}")
    print()
    
    settings = build_framework_settings(username, framework_home, 'unknown')
    
    # Step 3: Create output directories
    print("[3/6] Creating output directories...")
    output_dirs = [
        os.path.join(framework_home, 'output', 'RBHG'),
        os.path.join(framework_home, 'output', 'SPIZG'),
        settings['rbhg_output_base'],
        settings['spizg_output_base'],
        os.path.dirname(settings['real_data_output_base'].replace('{RunPeriod}', 'placeholder')),
        os.path.dirname(settings['real_data_log_base'].replace('{RunPeriod}', 'placeholder')),
    ]
    for output_dir in output_dirs:
        make_directory(output_dir)
    print()
    
    # Update RunPeriods.json cobrems paths
    print("Updating RunPeriods.json cobrems distribution paths...")
    update_runperiods_cobrems(framework_home)
    print()
    
    # Step 4: Detect shell type
    print("[4/6] Detecting shell type...")
    shell_type = detect_shell()
    print(f"  ✓ Shell type: {shell_type}")
    settings = build_framework_settings(username, framework_home, shell_type)
    print()
    
    # Step 5: Generate workflow defaults
    print("[5/6] Generating real-data DSelector defaults...")
    write_default_data_dselector_config(settings)
    print()
    
    # Step 6: Update config files
    print("[6/6] Updating configuration files...")
    config_dir = os.path.join(framework_home, 'GeneratorConfigExamples')
    config_files = glob.glob(os.path.join(config_dir, '*.config'))
    
    if not config_files:
        print(f"  WARNING: No .config files found in {config_dir}")
    else:
        for config_file in sorted(config_files):
            update_config_file(config_file, username, framework_home)
    print()
    
    # Save settings for downstream scripts
    print("Saving framework settings...")
    settings_path = save_framework_settings(settings)
    print(f"  ✓ Settings saved to: {settings_path}")
    print()
    
    # Summary
    print("="*70)
    print("Setup Complete!")
    print("="*70)
    print(f"""
Configuration Summary:
  Username:       {username}
  Framework Home: {framework_home}
  Shell Type:     {shell_type}

Your configuration files in GeneratorConfigExamples/ have been updated.
Downstream scripts will automatically use these settings.

You can now run the framework:
  ./swif2_gen.py GeneratorConfigExamples/RBHG.config
  ./swif2_gen.py GeneratorConfigExamples/SPIZG.config

To reconfigure, simply run this script again:
  ./setup_framework.py
""")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nSetup cancelled by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR during setup: {e}")
        sys.exit(1)
