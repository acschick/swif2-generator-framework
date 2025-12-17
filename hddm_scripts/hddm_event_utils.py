#!/usr/bin/env python3

"""
Enhanced event counting and validation utilities for HDDM files
Ensures event integrity throughout the analysis pipeline
"""

import os
import sys
import glob
import json
import subprocess
from pathlib import Path
import hddm_s

def count_hddm_events(hddm_filename, verbose=False):
    """
    Count events in a single HDDM file
    
    Args:
        hddm_filename (str): Path to HDDM file
        verbose (bool): Print progress messages
        
    Returns:
        int: Number of events in the file
    """
    if not os.path.exists(hddm_filename):
        print(f"ERROR: File not found: {hddm_filename}")
        return 0
    
    try:
        count = 0
        for i, _ in enumerate(hddm_s.istream(hddm_filename), 1):
            if verbose and i % 10000 == 0:
                print(f"...processed {i} events")
            count = i
        
        if verbose:
            print(f"Total events in {os.path.basename(hddm_filename)}: {count}")
        
        return count
        
    except Exception as e:
        print(f"ERROR reading {hddm_filename}: {e}")
        return 0

def count_directory_hddm_events(directory, pattern="vectors_*.hddm", verbose=False):
    """
    Count events in all HDDM files matching pattern in a directory
    
    Args:
        directory (str): Directory path
        pattern (str): File pattern to match
        verbose (bool): Print detailed information
        
    Returns:
        dict: File-by-file event counts and total
    """
    results = {
        'files': {},
        'total_events': 0,
        'total_files': 0,
        'directory': directory,
        'pattern': pattern
    }
    
    search_path = os.path.join(directory, pattern)
    hddm_files = glob.glob(search_path)
    
    if not hddm_files:
        print(f"WARNING: No files matching {pattern} found in {directory}")
        return results
    
    if verbose:
        print(f"Found {len(hddm_files)} HDDM files in {directory}")
    
    for hddm_file in sorted(hddm_files):
        basename = os.path.basename(hddm_file)
        event_count = count_hddm_events(hddm_file, verbose=False)
        results['files'][basename] = event_count
        results['total_events'] += event_count
        
        if verbose:
            print(f"  {basename}: {event_count} events")
    
    results['total_files'] = len(hddm_files)
    
    if verbose:
        print(f"Total: {results['total_events']} events in {results['total_files']} files")
    
    return results

def validate_expected_events(directory, expected_events_per_file, pattern="vectors_*.hddm", tolerance=0):
    """
    Validate that HDDM files contain the expected number of events
    
    Args:
        directory (str): Directory containing HDDM files
        expected_events_per_file (int): Expected events per file (e.g., 50000)
        pattern (str): File pattern to match
        tolerance (int): Acceptable deviation from expected count
        
    Returns:
        dict: Validation results
    """
    results = count_directory_hddm_events(directory, pattern, verbose=False)
    
    validation = {
        'valid': True,
        'total_expected': expected_events_per_file * results['total_files'],
        'total_actual': results['total_events'],
        'discrepancy': results['total_events'] - (expected_events_per_file * results['total_files']),
        'files_with_issues': [],
        'summary': {}
    }
    
    # Check each file
    for filename, actual_count in results['files'].items():
        expected = expected_events_per_file
        difference = actual_count - expected
        
        if abs(difference) > tolerance:
            validation['valid'] = False
            validation['files_with_issues'].append({
                'filename': filename,
                'expected': expected,
                'actual': actual_count,
                'difference': difference
            })
    
    # Overall validation
    if abs(validation['discrepancy']) > tolerance * results['total_files']:
        validation['valid'] = False
    
    validation['summary'] = {
        'total_files': results['total_files'],
        'expected_per_file': expected_events_per_file,
        'total_expected': validation['total_expected'],
        'total_actual': validation['total_actual'],
        'discrepancy': validation['discrepancy'],
        'files_with_issues': len(validation['files_with_issues']),
        'validation_passed': validation['valid']
    }
    
    return validation

def save_event_count_report(directory, output_file=None, pattern="vectors_*.hddm"):
    """
    Generate and save a detailed event count report
    
    Args:
        directory (str): Directory to analyze
        output_file (str): Output JSON file path (optional)
        pattern (str): File pattern to match
        
    Returns:
        str: Path to the saved report file
    """
    results = count_directory_hddm_events(directory, pattern, verbose=True)
    
    # Add metadata
    results['timestamp'] = str(Path().cwd())
    results['analysis_timestamp'] = os.popen('date').read().strip()
    
    # Determine output file
    if output_file is None:
        output_file = os.path.join(directory, "event_count_report.json")
    
    # Save report
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"Event count report saved to: {output_file}")
    return output_file

def merge_and_validate_hddm(input_dir, output_file, expected_events_per_file=50000, 
                          pattern="vectors_*.hddm", validate_before=True, validate_after=True):
    """
    Merge HDDM files with event count validation
    
    Args:
        input_dir (str): Directory containing input HDDM files
        output_file (str): Output merged HDDM file
        expected_events_per_file (int): Expected events per input file
        pattern (str): Pattern for input files
        validate_before (bool): Validate input files before merging
        validate_after (bool): Validate output file after merging
        
    Returns:
        dict: Merge and validation results
    """
    results = {
        'input_dir': input_dir,
        'output_file': output_file,
        'expected_events_per_file': expected_events_per_file,
        'merge_successful': False,
        'validation_passed': False,
        'input_validation': None,
        'output_validation': None,
        'event_counts': {}
    }
    
    # Pre-merge validation
    if validate_before:
        print("Validating input files...")
        input_validation = validate_expected_events(input_dir, expected_events_per_file, pattern)
        results['input_validation'] = input_validation
        
        if not input_validation['valid']:
            print("WARNING: Input validation failed!")
            print(f"Expected {input_validation['summary']['total_expected']} events, "
                  f"found {input_validation['summary']['total_actual']} events")
            
            for issue in input_validation['files_with_issues']:
                print(f"  {issue['filename']}: expected {issue['expected']}, "
                      f"got {issue['actual']} (diff: {issue['difference']})")
    
    # Get input file list and total expected events
    input_files = glob.glob(os.path.join(input_dir, pattern))
    if not input_files:
        print(f"ERROR: No files matching {pattern} found in {input_dir}")
        return results
    
    # Filter out empty files (0 events) before merging
    print(f"Checking {len(input_files)} files for empty files...")
    non_empty_files = []
    empty_files = []
    for f in input_files:
        event_count = count_hddm_events(f, verbose=False)
        if event_count > 0:
            non_empty_files.append(f)
        else:
            empty_files.append(os.path.basename(f))
    
    if empty_files:
        print(f"WARNING: Skipping {len(empty_files)} empty files: {', '.join(empty_files)}")
        missing_events = len(empty_files) * expected_events_per_file
        print(f"         Missing approximately {missing_events} events from failed jobs")
        print(f"         Consider re-running generator for these jobs or check SWIF2 logs for failures")
    
    if not non_empty_files:
        print(f"ERROR: No non-empty HDDM files found!")
        return results
    
    print(f"Using {len(non_empty_files)} non-empty files for merge")
    
    total_expected_events = len(input_files) * expected_events_per_file
    results['event_counts']['expected_total'] = total_expected_events
    results['event_counts']['input_files'] = len(input_files)
    results['event_counts']['non_empty_files'] = len(non_empty_files)
    results['event_counts']['empty_files'] = len(empty_files)
    
    # Perform merge using hddm_merge_files
    try:
        import subprocess
        # Note: hddm_merge_files uses -oOutputfile (no space!)
        cmd = ["hddm_merge_files", f"-o{output_file}"] + sorted(non_empty_files)
        print(f"Merging {len(non_empty_files)} HDDM files...")
        print(f"Command: hddm_merge_files -o{os.path.basename(output_file)} [... {len(non_empty_files)} files]")
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=input_dir)
        
        if result.returncode == 0:
            results['merge_successful'] = True
            print(f"Merge completed successfully: {output_file}")
        else:
            print(f"Merge failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return results
            
    except Exception as e:
        print(f"ERROR during merge: {e}")
        return results
    
    # Post-merge validation
    if validate_after and os.path.exists(output_file):
        print("Validating merged output file...")
        actual_events = count_hddm_events(output_file, verbose=True)
        results['event_counts']['output_actual'] = actual_events
        
        event_difference = actual_events - total_expected_events
        results['event_counts']['difference'] = event_difference
        
        if event_difference == 0:
            results['validation_passed'] = True
            print(f"SUCCESS: Event count verified ({actual_events} events)")
        else:
            print(f"WARNING: Event count mismatch!")
            print(f"Expected: {total_expected_events}, Actual: {actual_events}, "
                  f"Difference: {event_difference}")
    
    return results

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage:")
        print("  count_events.py <hddm_file>                    # Count single file")
        print("  count_events.py --dir <directory>              # Count all vectors_*.hddm in directory")
        print("  count_events.py --validate <directory> <expected_events_per_file>")
        print("  count_events.py --merge <input_dir> <output_file> [expected_events_per_file]")
        sys.exit(1)
    
    if sys.argv[1] == "--dir":
        # Count directory
        directory = sys.argv[2]
        results = count_directory_hddm_events(directory, verbose=True)
        save_event_count_report(directory)
        
    elif sys.argv[1] == "--validate":
        # Validate directory
        directory = sys.argv[2]
        expected = int(sys.argv[3]) if len(sys.argv) > 3 else 50000
        validation = validate_expected_events(directory, expected)
        
        print(f"Validation Results:")
        print(f"  Files checked: {validation['summary']['total_files']}")
        print(f"  Expected per file: {validation['summary']['expected_per_file']}")
        print(f"  Total expected: {validation['summary']['total_expected']}")
        print(f"  Total actual: {validation['summary']['total_actual']}")
        print(f"  Discrepancy: {validation['summary']['discrepancy']}")
        print(f"  Validation: {'PASSED' if validation['valid'] else 'FAILED'}")
        
        if validation['files_with_issues']:
            print(f"  Files with issues: {len(validation['files_with_issues'])}")
            for issue in validation['files_with_issues'][:5]:  # Show first 5
                print(f"    {issue['filename']}: {issue['actual']} (expected {issue['expected']})")
    
    elif sys.argv[1] == "--merge":
        # Merge and validate
        input_dir = sys.argv[2]
        output_file = sys.argv[3]
        expected = int(sys.argv[4]) if len(sys.argv) > 4 else 50000
        
        results = merge_and_validate_hddm(input_dir, output_file, expected)
        
        print("\nMerge Results:")
        print(f"  Merge successful: {results['merge_successful']}")
        print(f"  Validation passed: {results['validation_passed']}")
        if 'event_counts' in results:
            counts = results['event_counts']
            print(f"  Expected events: {counts.get('expected_total', 'unknown')}")
            print(f"  Actual events: {counts.get('output_actual', 'unknown')}")
            print(f"  Difference: {counts.get('difference', 'unknown')}")
    
    else:
        # Count single file
        hddm_file = sys.argv[1]
        count = count_hddm_events(hddm_file, verbose=True)
        print(f"Total events: {count}")