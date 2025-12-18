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
4. Update all .config files in GeneratorConfigExamples/
5. Save settings for downstream scripts (prepareSimulation.py, etc.)
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
            
            # Update FRAMEWORK_HOME (leave it blank but add comment showing auto-detected value)
            elif key == 'FRAMEWORK_HOME':
                # Keep it blank (auto-detects), but update comment to show auto-detected path
                updated_lines.append(f"FRAMEWORK_HOME               # Auto-detected: {framework_home}\n")
            
            else:
                # Keep line as-is
                updated_lines.append(line)
        else:
            updated_lines.append(line)
    
    # Write updated config
    with open(config_path, 'w') as f:
        f.writelines(updated_lines)
    
    print(f"  ✓ Updated {os.path.basename(config_path)}")


def save_framework_settings(username, framework_home, shell_type):
    """
    Save framework settings to a JSON file for other scripts to use.
    This allows prepareSimulation.py and other downstream scripts to
    automatically use the correct shell type, etc.
    """
    settings = {
        'username': username,
        'framework_home': framework_home,
        'shell_type': shell_type,
        'setup_complete': True
    }
    
    settings_path = os.path.join(framework_home, 'framework_settings.json')
    with open(settings_path, 'w') as f:
        json.dump(settings, indent=2, fp=f)
    
    return settings_path


def main():
    print("="*70)
    print("SWIF2 Generator Framework - One-Time Setup")
    print("="*70)
    print()
    
    # Step 1: Detect username
    print("[1/4] Detecting username...")
    username = detect_username()
    print(f"  ✓ Username: {username}")
    print()
    
    # Step 2: Detect framework home
    print("[2/4] Detecting framework installation directory...")
    framework_home = detect_framework_home()
    print(f"  ✓ Framework home: {framework_home}")
    print()
    
    # Step 3: Detect shell type
    print("[3/4] Detecting shell type...")
    shell_type = detect_shell()
    print(f"  ✓ Shell type: {shell_type}")
    print()
    
    # Step 4: Update config files
    print("[4/4] Updating configuration files...")
    config_dir = os.path.join(framework_home, 'GeneratorConfigExamples')
    config_files = glob.glob(os.path.join(config_dir, '*.config'))
    
    if not config_files:
        print(f"  WARNING: No .config files found in {config_dir}")
    else:
        for config_file in sorted(config_files):
            update_config_file(config_file, username, framework_home)
    print()
    
    # Step 5: Save settings for downstream scripts
    print("Saving framework settings...")
    settings_path = save_framework_settings(username, framework_home, shell_type)
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
