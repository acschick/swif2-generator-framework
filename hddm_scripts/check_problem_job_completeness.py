#!/usr/bin/env python3
"""
Check if "problem" jobs actually completed successfully
Compares event counts in vectors.txt vs .hddm files
"""

import os
import sys
import glob
from hddm_event_utils import count_hddm_events


def count_txt_events(txt_file):
    """Count lines (events) in a vectors txt file"""
    if not os.path.exists(txt_file):
        return 0
    
    with open(txt_file, 'r') as f:
        return sum(1 for line in f if line.strip())


def check_job_completion(vectors_dir, job_numbers):
    """Check if jobs are actually complete by comparing txt and hddm event counts
    
    Args:
        vectors_dir: Path to vectors directory
        job_numbers: List of job numbers to check
        
    Returns:
        dict: Results for each job with completion status
    """
    results = []
    
    for job_num in job_numbers:
        txt_file = os.path.join(vectors_dir, f"vectors_{job_num:04d}.txt")
        hddm_pattern = os.path.join(vectors_dir, f"*_{job_num:04d}.hddm")
        hddm_files = glob.glob(hddm_pattern)
        
        result = {
            'job_num': job_num,
            'txt_exists': os.path.exists(txt_file),
            'hddm_exists': len(hddm_files) > 0,
            'txt_events': 0,
            'hddm_events': 0,
            'hddm_file': None,
            'complete': False
        }
        
        # Count txt events
        if result['txt_exists']:
            result['txt_events'] = count_txt_events(txt_file)
        
        # Count hddm events
        if result['hddm_exists']:
            hddm_file = hddm_files[0]  # Should only be one
            result['hddm_file'] = os.path.basename(hddm_file)
            result['hddm_events'] = count_hddm_events(hddm_file, verbose=False)
        
        # Check if complete (both files exist and have matching event counts)
        if result['txt_exists'] and result['hddm_exists']:
            if result['txt_events'] == result['hddm_events'] and result['txt_events'] > 0:
                result['complete'] = True
        
        results.append(result)
    
    return results


def print_results(results):
    """Print formatted results table"""
    
    print("\n" + "="*90)
    print("JOB COMPLETION ANALYSIS")
    print("="*90)
    print(f"{'Job#':<8} {'TXT':<12} {'HDDM':<12} {'Match':<8} {'Status':<15} {'HDDM File'}")
    print("-"*90)
    
    complete_count = 0
    incomplete_count = 0
    
    for r in results:
        txt_str = f"{r['txt_events']} events" if r['txt_exists'] else "MISSING"
        hddm_str = f"{r['hddm_events']} events" if r['hddm_exists'] else "MISSING"
        
        if r['complete']:
            match = "✓ YES"
            status = "COMPLETE"
            complete_count += 1
        elif r['txt_exists'] and r['hddm_exists']:
            match = "✗ NO"
            status = "MISMATCH"
            incomplete_count += 1
        elif r['txt_exists'] and not r['hddm_exists']:
            match = "-"
            status = "TXT ONLY"
            incomplete_count += 1
        elif r['hddm_exists'] and not r['txt_exists']:
            match = "-"
            status = "HDDM ONLY"
            incomplete_count += 1
        else:
            match = "-"
            status = "NO OUTPUT"
            incomplete_count += 1
        
        hddm_name = r['hddm_file'] if r['hddm_file'] else "-"
        
        print(f"{r['job_num']:<8} {txt_str:<12} {hddm_str:<12} {match:<8} {status:<15} {hddm_name}")
    
    print("="*90)
    print(f"\nSummary:")
    print(f"  Complete jobs (txt + hddm match):     {complete_count}")
    print(f"  Incomplete/mismatched jobs:           {incomplete_count}")
    print(f"  Total jobs analyzed:                  {len(results)}")
    
    if complete_count > 0:
        print(f"\n  NOTE: {complete_count} job(s) appear to have completed successfully despite")
        print(f"        SLURM_TIMEOUT. They may not need to be retried.")
        print(f"        Consider keeping their output and only retrying incomplete jobs.")
    
    print("="*90)
    
    return complete_count, incomplete_count


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Check if problem jobs actually completed',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  # Check specific job numbers:
  ./hddm_scripts/check_problem_job_completeness.py output/RBHG/CPP_FFS1/BC_FFN/2205_135DEG_FFN_IRADOFF/vectors \\
    1809 2225 2253 2610 2803 3766 3844 6462
        """
    )
    
    parser.add_argument('vectors_dir', help='Path to vectors directory')
    parser.add_argument('job_numbers', nargs='+', type=int, help='Job numbers to check')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.vectors_dir):
        print(f"ERROR: Directory not found: {args.vectors_dir}")
        return 1
    
    print(f"Analyzing jobs in: {args.vectors_dir}")
    print(f"Checking {len(args.job_numbers)} job(s)...")
    
    results = check_job_completion(args.vectors_dir, args.job_numbers)
    complete_count, incomplete_count = print_results(results)
    
    return 0 if incomplete_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
