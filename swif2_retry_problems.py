#!/usr/bin/env python3
"""
SWIF2 Problem Job Retry Script
================================
Automatically retry problem jobs in a SWIF2 workflow after cleaning up partial output.

Usage:
  # Single workflow - dry run:
  ./swif2_retry_problems.py gen_CPP_FFS1_BCFFN_2205_135DEG_IRADOFF_ee
  
  # Single workflow - execute:
  ./swif2_retry_problems.py gen_CPP_FFS1_BCFFN_2205_135DEG_IRADOFF_ee --execute
  
  # All workflows matching pattern:
  ./swif2_retry_problems.py gen_CPP --pattern --execute
  ./swif2_retry_problems.py gen_CPP_FFS1 --pattern --execute
  
  # With vectors cleanup (single workflow only):
  ./swif2_retry_problems.py gen_CPP_FFS1_BCFFN_2205_135DEG_IRADOFF_ee \
    --vectors output/RBHG/CPP_FFS1/BC_FFN/2205_135DEG_FFN_IRADOFF/vectors --execute
"""

import os
import subprocess
import argparse
import glob

# Import completion checking functions
try:
    from hddm_scripts.hddm_event_utils import count_hddm_events
    HDDM_AVAILABLE = True
except ImportError:
    count_hddm_events = None
    HDDM_AVAILABLE = False


def run_command(cmd, capture=True):
    """Run shell command"""
    try:
        if capture:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
            return result.stdout.strip()
        else:
            subprocess.run(cmd, shell=True, check=True)
            return True  # Return True on success when not capturing output
    except subprocess.CalledProcessError as e:
        print(f"ERROR: {cmd}")
        if e.stderr:
            print(f"{e.stderr}")
        return None


def get_problem_attempts(workflow_name):
    """Get problem attempt IDs, job numbers, and problem types using swif2 status -problems
    
    Returns:
        list of dicts with keys: 'attempt_id', 'job_number', 'problem_type'
    """
    cmd = f"swif2 status {workflow_name} -problems"
    output = run_command(cmd)
    
    if not output:
        return []
    
    # Parse job IDs, job names, and problem types from output
    # Format: 
    #   job_id                        = 54033459
    #   job_name                      = gen_CPP_FFS1_BCFFN_2205_135DEG_IRADOFF_ee_1809
    #   job_attempt_problem           = SLURM_TIMEOUT
    
    jobs = []
    current_job = {}
    
    for line in output.split('\n'):
        line = line.strip()
        
        if line.startswith('job_id'):
            # Start a new job entry
            if current_job:
                jobs.append(current_job)
            current_job = {}
            
            parts = line.split('=')
            if len(parts) == 2:
                job_id = parts[1].strip()
                if job_id.isdigit():
                    current_job['attempt_id'] = job_id
        
        elif line.startswith('job_name'):
            parts = line.split('=')
            if len(parts) == 2:
                job_name = parts[1].strip()
                job_num_str = job_name.split('_')[-1]
                if job_num_str.isdigit():
                    current_job['job_number'] = int(job_num_str)
        
        elif line.startswith('job_attempt_problem'):
            parts = line.split('=')
            if len(parts) == 2:
                problem_type = parts[1].strip()
                current_job['problem_type'] = problem_type
    
    # Don't forget the last job
    if current_job:
        jobs.append(current_job)
    
    return jobs


def derive_vectors_path(workflow_name):
    """Derive vectors directory path from workflow name.
    
    Workflow format: gen_{FFS}_{DATASET}_{RUNPERIOD}_{POLARIZATION}_{EXTRAS}_{LEPTON}
    Example: gen_CPP_FFS1_BCFFN_2205_135DEG_IRADOFF_ee
    Path: output/RBHG/CPP_FFS1/BC_FFN/2205_135DEG_FFN_IRADOFF/vectors/
    """
    parts = workflow_name.split('_')
    
    if len(parts) < 7 or not parts[0] == 'gen':
        return None
    
    # Parse workflow components
    ffs = '_'.join(parts[1:3])  # CPP_FFS1
    dataset_raw = parts[3]  # BCFFN, BCFF1, SIMFF1, qDATAqFFN
    
    # Map dataset to directory name with underscore
    dataset_map = {
        'BCFFN': 'BC_FFN',
        'BCFF1': 'BC_FF1',
        'SIMFF1': 'SIM_FF1',
        'qDATAqFFN': 'qDATAq_FFN'
    }
    dataset = dataset_map.get(dataset_raw, dataset_raw)
    
    # Extract form factor from dataset
    if 'FFN' in dataset:
        ff = 'FFN'
    elif 'FF1' in dataset:
        ff = 'FF1'
    else:
        ff = 'FFN'  # default
    
    # Reconstruct run directory name: {RUNPERIOD}_{POLARIZATION}_{FF}_{EXTRAS}
    # From gen_CPP_FFS1_BCFFN_2205_135DEG_IRADOFF_ee
    # We want: 2205_135DEG_FFN_IRADOFF
    run_config_parts = parts[4:-1]  # ['2205', '135DEG', 'IRADOFF']
    run_dir = '_'.join([run_config_parts[0], run_config_parts[1], ff] + run_config_parts[2:])
    
    vectors_path = os.path.join('output', 'RBHG', ffs, dataset, run_dir, 'vectors')
    return vectors_path


def get_matching_workflows(pattern):
    """Get list of workflows matching a pattern from swif2 list"""
    print(f"Finding workflows matching: {pattern}*")
    
    # Use swif2 list to get all workflows for current user
    cmd = "swif2 list"
    output = run_command(cmd)
    
    if not output:
        return []
    
    matching = []
    for line in output.split('\n'):
        # Parse workflow name from first column (before first whitespace or |)
        parts = line.split()
        if parts:
            workflow_name = parts[0].strip()
            if workflow_name.startswith(pattern) and not workflow_name.startswith('#'):
                matching.append(workflow_name)
    
    return matching


def cleanup_vectors_files(vectors_dir, job_nums, dry_run=True):
    """Clean up partial vectors.txt and .hddm files for failed jobs"""
    if not vectors_dir:
        return []
        
    if not os.path.exists(vectors_dir):
        print(f"  Warning: Vectors directory not found: {vectors_dir}")
        return []
    
    print(f"  Looking for partial outputs for job numbers: {', '.join(str(n) for n in job_nums)}")
    
    files_deleted = []
    files_to_delete = []
    
    # First collect all files that would be deleted
    for job_num in job_nums:
        # Look for vectors_####.txt pattern
        pattern_txt = os.path.join(vectors_dir, f"vectors_{job_num:04d}.txt")
        pattern_hddm = os.path.join(vectors_dir, f"*_{job_num:04d}.hddm")
        
        if os.path.exists(pattern_txt):
            files_to_delete.append(pattern_txt)
        
        files_to_delete.extend(glob.glob(pattern_hddm))
    
    if not files_to_delete:
        print(f"  No partial output files found in: {vectors_dir}")
        return []
    
    # Now delete or show what would be deleted
    for filepath in files_to_delete:
        if dry_run:
            print(f"  Would delete: {filepath}")
        else:
            try:
                os.remove(filepath)
                print(f"  Deleted: {filepath}")
                files_deleted.append(filepath)
            except OSError as e:
                print(f"  ERROR deleting {filepath}: {e}")
    
    return files_deleted


def count_txt_events(txt_file):
    """Count lines (events) in a vectors txt file"""
    if not os.path.exists(txt_file):
        return 0
    
    with open(txt_file, 'r', encoding='utf-8') as f:
        return sum(1 for line in f if line.strip())


def check_job_completion(vectors_dir, job_numbers):
    """Check which jobs are actually complete
    
    Args:
        vectors_dir: Path to vectors directory
        job_numbers: List of job numbers to check
    
    Returns:
        tuple: (complete_jobs, incomplete_jobs) - lists of job numbers, or ([], job_numbers) if skipped
    """
    if not vectors_dir or not os.path.exists(vectors_dir):
        return [], list(job_numbers)
    
    if not HDDM_AVAILABLE:
        print("  Warning: hddm_s not available, skipping completion check")
        return [], list(job_numbers)
    
    # First pass: count HDDM files and estimate work
    hddm_files_to_check = []
    total_size_bytes = 0
    
    for job_num in job_numbers:
        hddm_pattern = os.path.join(vectors_dir, f"*_{job_num:04d}.hddm")
        hddm_files = glob.glob(hddm_pattern)
        if hddm_files:
            hddm_file = hddm_files[0]
            hddm_files_to_check.append((job_num, hddm_file))
            if os.path.exists(hddm_file):
                total_size_bytes += os.path.getsize(hddm_file)
    
    num_hddm_files = len(hddm_files_to_check)
    
    # Decide if we need to prompt
    # Threshold: >10 files OR total size > 500MB
    size_threshold_mb = 500
    count_threshold = 10
    total_size_mb = total_size_bytes / (1024 * 1024)
    
    needs_confirmation = (num_hddm_files > count_threshold or total_size_mb > size_threshold_mb)
    
    if needs_confirmation:
        print(f"\n  Completion check will read {num_hddm_files} HDDM file(s) (~{total_size_mb:.1f} MB)")
        print("  This may take a while depending on file sizes.")
        response = input("  Proceed with completion check? [Y/n]: ").strip().lower()
        
        if response and response not in ['y', 'yes']:
            print("  Skipping completion check. All jobs will be treated as incomplete.")
            return [], list(job_numbers)
    
    # Proceed with checking
    complete = []
    incomplete = []
    
    print(f"  Checking completion for {len(job_numbers)} job(s)...")
    
    for idx, job_num in enumerate(job_numbers, 1):
        txt_file = os.path.join(vectors_dir, f"vectors_{job_num:04d}.txt")
        hddm_pattern = os.path.join(vectors_dir, f"*_{job_num:04d}.hddm")
        hddm_files = glob.glob(hddm_pattern)
        
        # Show progress for every job when checking HDDM files
        if hddm_files or (idx % 5 == 0 or idx == len(job_numbers)):
            print(f"    [{idx}/{len(job_numbers)}] Checking job {job_num}...", end='\r')
        
        txt_events = 0
        hddm_events = 0
        
        if os.path.exists(txt_file):
            txt_events = count_txt_events(txt_file)
        
        if hddm_files:
            hddm_events = count_hddm_events(hddm_files[0], verbose=False)
        
        # Complete if both exist and match
        if txt_events > 0 and hddm_events > 0 and txt_events == hddm_events:
            complete.append(job_num)
        else:
            incomplete.append(job_num)
    
    # Clear progress line
    print(f"    Completed checking {len(job_numbers)} job(s).{' ' * 20}")
    
    return complete, incomplete


# Configuration for problem-specific resource adjustments
PROBLEM_RESOURCE_FIXES = {
    'SLURM_TIMEOUT': {'time': ('mult', 2)},
    'SLURM_OUT_OF_MEMORY': {'ram': ('mult', 2)},
    'SLURM_NODE_FAIL': {},  # Just retry, no modification needed
    'SLURM_DISK_QUOTA': {'disk': ('mult', 1.5)},
    # Add more problem types as needed
}


def modify_jobs_resources(workflow_name, attempt_ids, problem_type, dry_run=True):
    """Modify job resources based on problem type using swif2 modify-jobs
    
    This will resolve the problem and re-dispatch with adjusted resources.
    """
    if not attempt_ids:
        return False
    
    # Get resource modifications for this problem type
    modifications = PROBLEM_RESOURCE_FIXES.get(problem_type, {})
    
    if not modifications:
        # Unknown problem type or no modifications needed, use basic retry
        return None
    
    # Build modification arguments
    mod_args = []
    for resource, (operation, value) in modifications.items():
        # Format value as integer if it's a whole number, otherwise use the value
        if isinstance(value, float) and value.is_integer():
            formatted_value = str(int(value))
        else:
            formatted_value = str(value)
        mod_args.extend([f'-{resource}', operation, formatted_value])
    
    attempts_str = ' '.join(attempt_ids)
    cmd = f"swif2 modify-jobs {workflow_name} {' '.join(mod_args)} {attempts_str}"
    
    if dry_run:
        print(f"\nWould execute: {cmd}")
        return True
    else:
        result = run_command(cmd, capture=False)
        return result is not None


def bless_jobs(workflow_name, attempt_ids, dry_run=True):
    """Bless jobs using swif2 bless-jobs (mark as successful without re-running)"""
    if not attempt_ids:
        return False
    
    attempts_str = ' '.join(attempt_ids)
    cmd = f"swif2 bless-jobs {workflow_name} {attempts_str}"
    
    if dry_run:
        print(f"\nWould execute: {cmd}")
        return True
    else:
        result = run_command(cmd, capture=False)
        return result is not None


def retry_jobs(workflow_name, attempt_ids, dry_run=True):
    """Retry jobs using swif2 retry-jobs"""
    if not attempt_ids:
        return False
    
    attempts_str = ' '.join(attempt_ids)
    cmd = f"swif2 retry-jobs {workflow_name} {attempts_str}"
    
    if dry_run:
        print(f"\nWould execute: {cmd}")
        return False
    else:
        print(f"\nRetrying {len(attempt_ids)} job(s)...")
        result = run_command(cmd, capture=False)
        return result is not None


def main():
    parser = argparse.ArgumentParser(
        description='Retry problem jobs in SWIF2 workflow',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single workflow - dry run with completion checking (default):
  ./swif2_retry_problems.py gen_CPP_FFS1_BCFFN_2205_135DEG_IRADOFF_ee
  
  # Execute with smart completion checking and resource adjustments:
  ./swif2_retry_problems.py gen_CPP_FFS1_BCFFN_2205_135DEG_IRADOFF_ee --execute
  
  # All workflows starting with pattern:
  ./swif2_retry_problems.py gen_CPP --pattern --execute
  
  # Disable completion checking (retry all problem jobs):
  ./swif2_retry_problems.py gen_CPP_FFS1_BCFFN_2205_135DEG_IRADOFF_ee --no-check-complete --execute
  
  # With manual vectors path:
  ./swif2_retry_problems.py gen_CPP_FFS1_BCFFN_2205_135DEG_IRADOFF_ee \
    --vectors output/RBHG/CPP_FFS1/BC_FFN/2205_135DEG_FFN_IRADOFF/vectors --execute
        """
    )
    
    parser.add_argument('workflow', help='SWIF2 workflow name or pattern prefix')
    parser.add_argument('--pattern', action='store_true',
                       help='Treat workflow argument as a pattern prefix (match all workflows starting with it)')
    parser.add_argument('--execute', action='store_true',
                       help='Actually retry jobs (default: dry run)')
    parser.add_argument('--vectors', nargs='?', const=True, default=None,
                       help='Clean up partial output files. Optionally specify vectors directory path, '
                            'otherwise auto-detected from workflow name (single workflow only)')
    parser.add_argument('--no-cleanup', action='store_true',
                       help='Skip cleanup of partial output files')
    parser.add_argument('--no-check-complete', action='store_true',
                       help='Disable completion checking. By default, script checks if jobs are actually '
                            'complete (matching txt/hddm event counts) and blesses them instead of retrying. '
                            'Use this flag to skip checking and retry all problem jobs.')
    parser.add_argument('--no-smart-retry', action='store_true',
                       help='Disable smart resource modifications. All jobs will use basic retry '
                            'instead of adjusting time/RAM/disk based on problem type.')
    parser.add_argument('--time-mult', type=float, default=2.0,
                       help='Time multiplier for SLURM_TIMEOUT (default: 2.0)')
    parser.add_argument('--ram-mult', type=float, default=2.0,
                       help='RAM multiplier for SLURM_OUT_OF_MEMORY (default: 2.0)')
    parser.add_argument('--disk-mult', type=float, default=1.5,
                       help='Disk multiplier for SLURM_DISK_QUOTA (default: 1.5)')
    
    args = parser.parse_args()
    
    # Determine if we should check completion (default: True unless disabled)
    check_complete = not args.no_check_complete
    
    # Update resource fix configuration based on arguments
    if not args.no_smart_retry:
        PROBLEM_RESOURCE_FIXES['SLURM_TIMEOUT']['time'] = ('mult', args.time_mult)
        PROBLEM_RESOURCE_FIXES['SLURM_OUT_OF_MEMORY']['ram'] = ('mult', args.ram_mult)
        PROBLEM_RESOURCE_FIXES['SLURM_DISK_QUOTA']['disk'] = ('mult', args.disk_mult)
    else:
        # Disable all smart modifications
        for key in PROBLEM_RESOURCE_FIXES:
            PROBLEM_RESOURCE_FIXES[key] = {}
    
    args = parser.parse_args()
    
    # Check if completion checking is usable
    if check_complete and not HDDM_AVAILABLE:
        print("WARNING: Completion checking requires hddm_s module which is not available.")
        print("         All problem jobs will be retried without completion checking.")
        print("         Use --no-check-complete to suppress this warning.")
        print()
        check_complete = False
    
    dry_run = not args.execute
    
    if dry_run:
        print("="*70)
        print("DRY RUN - No changes will be made")
        print("Add --execute to actually retry jobs")
        print("="*70)
        print()
    
    # Determine workflows to process
    if args.pattern:
        workflows = get_matching_workflows(args.workflow)
        if not workflows:
            print(f"No workflows found matching pattern: {args.workflow}")
            return
        print(f"Found {len(workflows)} matching workflow(s):")
        for wf in workflows:
            print(f"  - {wf}")
        print()
    else:
        workflows = [args.workflow]
    
    # Process each workflow
    total_problems = 0
    total_retried = 0
    
    for workflow_name in workflows:
        if len(workflows) > 1:
            print("\n" + "="*70)
            print(f"Processing: {workflow_name}")
            print("="*70)
        
        print(f"Finding problem jobs in workflow: {workflow_name}")
        
        # Get problem job details (returns list of dicts with attempt_id, job_number, problem_type)
        problem_jobs = get_problem_attempts(workflow_name)
        
        if not problem_jobs:
            print("  No problem jobs found.")
            continue
        
        # Extract lists for easier handling
        attempt_ids = [job['attempt_id'] for job in problem_jobs]
        job_nums = [job['job_number'] for job in problem_jobs]
        
        total_problems += len(problem_jobs)
        
        # Determine vectors directory
        vectors_dir = None
        if args.vectors or check_complete:
            # Auto-detect vectors directory if --vectors used without path
            if args.vectors is True or check_complete:
                vectors_dir = derive_vectors_path(workflow_name)
                if vectors_dir:
                    print(f"  Auto-detected vectors directory: {vectors_dir}")
                else:
                    print("  Warning: Could not auto-detect vectors directory from workflow name")
            else:
                vectors_dir = args.vectors
        
        # Check completion status (enabled by default)
        complete_jobs = []
        incomplete_jobs = []
        
        if check_complete and vectors_dir:
            print("  Checking job completion status...")
            complete_job_nums, incomplete_job_nums = check_job_completion(vectors_dir, job_nums)
            
            # Separate complete from incomplete jobs
            complete_jobs = [job for job in problem_jobs if job['job_number'] in complete_job_nums]
            incomplete_jobs = [job for job in problem_jobs if job['job_number'] in incomplete_job_nums]
            
            if complete_jobs:
                complete_nums = ', '.join(str(j['job_number']) for j in complete_jobs)
                print(f"  ✓ {len(complete_jobs)} job(s) are actually COMPLETE (will bless): {complete_nums}")
            if incomplete_jobs:
                incomplete_nums = ', '.join(str(j['job_number']) for j in incomplete_jobs)
                print(f"  ✗ {len(incomplete_jobs)} job(s) are INCOMPLETE (will handle): {incomplete_nums}")
        else:
            # No completion check - treat all as incomplete
            incomplete_jobs = problem_jobs
        
        # Clean up vectors files for incomplete jobs only
        if vectors_dir and not args.no_cleanup and incomplete_jobs:
            incomplete_nums = [job['job_number'] for job in incomplete_jobs]
            print(f"  Checking for partial outputs in: {vectors_dir}")
            cleanup_vectors_files(vectors_dir, incomplete_nums, dry_run)
        elif args.vectors and not args.no_cleanup and not incomplete_jobs:
            print("  No incomplete jobs to clean up.")
        elif not args.vectors and not check_complete and not args.no_cleanup and len(workflows) == 1:
            print("  Note: No vectors directory specified. Skipping cleanup and completion checking.")
            print("  Use --vectors to enable auto-detection and completion verification.")
        
        # Bless complete jobs
        if complete_jobs:
            complete_attempt_ids = [job['attempt_id'] for job in complete_jobs]
            if bless_jobs(workflow_name, complete_attempt_ids, dry_run):
                print(f"  ✓ {len(complete_jobs)} complete job(s) blessed")
            elif not dry_run:
                print("  ✗ Failed to bless complete jobs")
        
        # Handle incomplete jobs - group by problem type and apply smart fixes
        if incomplete_jobs:
            # Group incomplete jobs by problem type
            incomplete_by_problem = {}
            for job in incomplete_jobs:
                ptype = job.get('problem_type', 'UNKNOWN')
                if ptype not in incomplete_by_problem:
                    incomplete_by_problem[ptype] = []
                incomplete_by_problem[ptype].append(job)
            
            print(f"\n  Processing {len(incomplete_jobs)} incomplete job(s):")
            
            jobs_modified = 0
            jobs_retried = 0
            
            for problem_type, pjobs in sorted(incomplete_by_problem.items()):
                attempt_ids = [job['attempt_id'] for job in pjobs]
                job_nums_display = ', '.join(str(j['job_number']) for j in pjobs[:5])
                if len(pjobs) > 5:
                    job_nums_display += f", ... ({len(pjobs)} total)"
                
                # Check if we have a smart resource fix for this problem type
                if problem_type in PROBLEM_RESOURCE_FIXES and PROBLEM_RESOURCE_FIXES[problem_type]:
                    # Use modify-jobs with resource adjustments
                    modifications = PROBLEM_RESOURCE_FIXES[problem_type]
                    mod_desc = ', '.join(f"{k} {op} {v}" for k, (op, v) in modifications.items())
                    print(f"    {problem_type}: Modifying {len(pjobs)} job(s) ({mod_desc})")
                    print(f"      Jobs: {job_nums_display}")
                    
                    result = modify_jobs_resources(workflow_name, attempt_ids, problem_type, dry_run)
                    if result:
                        jobs_modified += len(pjobs)
                        total_retried += len(pjobs)
                    elif result is False and not dry_run:
                        print("      \u2717 Failed to modify jobs")
                elif problem_type in PROBLEM_RESOURCE_FIXES:
                    # Problem type is known but no modifications needed, just retry
                    print(f"    {problem_type}: Retrying {len(pjobs)} job(s) (no resource changes needed)")
                    print(f"      Jobs: {job_nums_display}")
                    
                    if retry_jobs(workflow_name, attempt_ids, dry_run):
                        jobs_retried += len(pjobs)
                        total_retried += len(pjobs)
                    elif not dry_run:
                        print("      \u2717 Failed to retry jobs")
                else:
                    # Unknown problem type, use basic retry
                    print(f"    {problem_type}: Retrying {len(pjobs)} job(s) (unknown problem type)")
                    print(f"      Jobs: {job_nums_display}")
                    
                    if retry_jobs(workflow_name, attempt_ids, dry_run):
                        jobs_retried += len(pjobs)
                        total_retried += len(pjobs)
                    elif not dry_run:
                        print("      ✗ Failed to retry jobs")
            
            if jobs_modified > 0:
                print(f"\n  ✓ {jobs_modified} job(s) modified with resource adjustments")
            if jobs_retried > 0:
                print(f"  ✓ {jobs_retried} job(s) retried without modifications")
        else:
            print("  No incomplete jobs to retry.")
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Workflows processed: {len(workflows)}")
    print(f"Total problem attempts: {total_problems}")
    print(f"Total retried: {total_retried}")
    
    if dry_run:
        print("\nThis was a dry run. Use --execute to actually retry jobs.")
    elif total_retried > 0:
        print("\nJobs have been resubmitted. Check status with:")
        print("  swif2 status")
    print("="*70)


if __name__ == "__main__":
    main()
