#!/usr/bin/env python3
"""
Python-based HDDM file merger
Avoids the C++ hddm_merge_files bug that corrupts momentum_double data
"""

import sys
import os
import glob
import time
import argparse
import hddm_s

def format_time(seconds):
    """Format seconds into human-readable time"""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f}m"
    else:
        return f"{seconds/3600:.1f}h"

def merge_hddm_files(input_files, output_file, new_run_number=None, renumber_events=True):
    """Merge multiple HDDM files into a single output file"""
    
    print(f"Output file: {output_file}")
    print(f"Input files: {len(input_files)}")
    if new_run_number is not None:
        print(f"Changing run number to: {new_run_number}")
    if renumber_events:
        print(f"Renumbering events sequentially")
    
    # Open output stream
    fout = hddm_s.ostream(output_file)
    
    total_events = 0
    start_time = time.time()
    
    # Process each input file
    for i, input_file in enumerate(input_files):
        file_start = time.time()
        
        # Get just the filename for display
        filename = os.path.basename(input_file)
        
        print(f"\n[{i+1}/{len(input_files)}] Processing: {filename}")
        
        try:
            fin = hddm_s.istream(input_file)
        except Exception as e:
            print(f"  ERROR: Cannot open {input_file}: {e}")
            continue
        
        events_from_file = 0
        
        # Read and write each event
        for rec in fin:
            # Modify event properties if requested
            for phys_evt in rec.getPhysicsEvents():
                if new_run_number is not None:
                    phys_evt.runNo = new_run_number
                if renumber_events:
                    phys_evt.eventNo = total_events + 1
            
            fout.write(rec)
            events_from_file += 1
            total_events += 1
        
        file_time = time.time() - file_start
        print(f"  {events_from_file} events in {file_time:.1f}s ({events_from_file/file_time:.0f} events/s)")
        
        # Estimate remaining time
        elapsed = time.time() - start_time
        files_remaining = len(input_files) - (i + 1)
        if i > 0:  # Need at least one file to estimate
            avg_time_per_file = elapsed / (i + 1)
            est_remaining = avg_time_per_file * files_remaining
            print(f"  Progress: {i+1}/{len(input_files)} files, {total_events} total events")
            print(f"  Estimated time remaining: {format_time(est_remaining)}")
    
    total_time = time.time() - start_time
    
    print(f"\n{'='*70}")
    print(f"Merge complete!")
    print(f"  Total events: {total_events}")
    print(f"  Total time: {format_time(total_time)}")
    print(f"  Average rate: {total_events/total_time:.0f} events/s")
    print(f"  Output: {output_file}")
    print(f"{'='*70}")
    
    return total_events


def main():
    parser = argparse.ArgumentParser(
        description='Merge HDDM files (avoids C++ hddm_merge_files bug)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s merged.hddm file1.hddm file2.hddm file3.hddm
  %(prog)s merged.hddm /path/to/vectors/
  %(prog)s merged.hddm /path/to/vectors/ --run 40856
        """)
    
    parser.add_argument('output', help='Output HDDM file')
    parser.add_argument('inputs', nargs='+', help='Input HDDM file(s) or directory')
    parser.add_argument('--run', type=int, help='Change run number to this value')
    parser.add_argument('--no-renumber', action='store_true', 
                        help='Keep original event numbers (default: renumber sequentially)')
    
    args = parser.parse_args()
    
    output_file = args.output
    
    # Check if input is a directory
    if len(args.inputs) == 1 and os.path.isdir(args.inputs[0]):
        input_dir = args.inputs[0]
        print(f"Scanning directory: {input_dir}")
        # Find all .hddm files in the directory
        pattern = os.path.join(input_dir, "*.hddm")
        input_files = sorted(glob.glob(pattern))
        if not input_files:
            print(f"ERROR: No .hddm files found in {input_dir}")
            sys.exit(1)
        print(f"Found {len(input_files)} HDDM files")
    else:
        input_files = args.inputs
    
    print("="*70)
    print("Python HDDM Merger")
    print("="*70)
    
    merge_hddm_files(input_files, output_file, args.run, renumber_events=not args.no_renumber)


if __name__ == "__main__":
    main()
