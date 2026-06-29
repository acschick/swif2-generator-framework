#!/bin/bash

# Usage function
usage() {
    echo "Usage: $0 [OPTIONS] [directory]"
    echo ""
    echo "Merge ROOT files by run number using hadd"
    echo ""
    echo "Options:"
    echo "  -c, --clean         Delete original batch files after merging"
    echo "  --cleanup-only      Only clean up pre-merged files (don't merge)"
    echo "  -n, --dry-run       Show what would be merged without running hadd"
    echo "  -h, --help          Show this help message"
    echo ""
    echo "Arguments:"
    echo "  directory           Target directory (default: current directory)"
    exit 0
}

# Function to clean up pre-merged files
cleanup_premerged_files() {
    local dir="$1"
    shift  # Remove first argument, remaining are run numbers
    local runs_to_clean=("$@")
    
    echo "Cleaning up pre-merged files in $dir..."
    local count=0
    
    # If specific runs provided, only clean those
    if [ ${#runs_to_clean[@]} -gt 0 ]; then
        for run in "${runs_to_clean[@]}"; do
            for file in "$dir"/*_${run}_[0-9][0-9][0-9].root; do
                [ -e "$file" ] || continue
                echo "  Removing: $(basename "$file")"
                rm "$file"
                ((count++))
            done
        done
    else
        # Clean all batch files (only for --cleanup-only mode)
        for file in "$dir"/*_[0-9][0-9][0-9].root; do
            [ -e "$file" ] || continue
            echo "  Removing: $(basename "$file")"
            rm "$file"
            ((count++))
        done
    fi
    
    echo "Removed $count pre-merged files"
}

# Parse arguments
clean_after=false
cleanup_only=false
dry_run=false
target_dir="."

while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--clean)
            clean_after=true
            shift
            ;;
        --cleanup-only)
            cleanup_only=true
            shift
            ;;
        -n|--dry-run)
            dry_run=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        -*)
            echo "Error: Unknown option $1"
            usage
            ;;
        *)
            target_dir="$1"
            shift
            ;;
    esac
done

# Check if directory exists
if [ ! -d "$target_dir" ]; then
    echo "Error: Directory '$target_dir' does not exist"
    exit 1
fi

# Check if hadd command is available (unless in cleanup-only or dry-run mode)
if [ "$cleanup_only" = false ] && [ "$dry_run" = false ]; then
    if ! command -v hadd &> /dev/null; then
        echo "Error: hadd command not found. Please set up your ROOT environment."
        echo "  Example: source /path/to/root/bin/thisroot.sh"
        exit 1
    fi
fi

# If cleanup-only mode, just clean and exit
if [ "$cleanup_only" = true ]; then
    if [ "$dry_run" = true ]; then
        echo "DRY RUN: Would clean up pre-merged files in $target_dir"
        echo "Files that would be removed:"
        for file in "$target_dir"/*_[0-9][0-9][0-9].root; do
            [ -e "$file" ] || continue
            echo "  - $(basename "$file")"
        done
    else
        # Verify merged files exist before cleaning
        merged_count=0
        for file in "$target_dir"/*.root; do
            [ -e "$file" ] || continue
            basename_file=$(basename "$file")
            # Count files that don't have batch numbers (merged files)
            if ! echo "$basename_file" | grep -q '_[0-9]\{3\}\.root$'; then
                ((merged_count++))
            fi
        done
        
        if [ $merged_count -eq 0 ]; then
            echo "ERROR: No merged files found in $target_dir"
            echo "  Refusing to clean up - merge files first to be safe!"
            exit 1
        fi
        
        echo "Found $merged_count merged files, proceeding with cleanup..."
        cleanup_premerged_files "$target_dir"
    fi
    exit 0
fi

declare -A groups
declare -A successfully_merged_runs

for file in "$target_dir"/*.root; do
    # Skip if no .root files found
    [ -e "$file" ] || continue
    
    basename_file=$(basename "$file")
    
    # Extract run number by finding the digits before the final _XXX.root pattern
    # This works for any filename pattern: tree_epemmissprot__B2_041151_000.root or tree_thrown_041212_006.root
    run=$(echo "$basename_file" | sed -n 's/.*_\([0-9]\+\)_[0-9]\{3\}\.root$/\1/p')
    
    # Skip files that don't match the pattern (e.g., already merged files)
    [ -z "$run" ] && continue
    
    groups[$run]="${groups[$run]} $file"
done

if [ "$dry_run" = true ]; then
    echo "DRY RUN: Would merge the following groups in $target_dir:"
    echo "=========================================================="
fi

for run in "${!groups[@]}"; do
    # Get the first file in the group to construct output name
    first_file=$(echo ${groups[$run]} | awk '{print $1}')
    
    # Remove the batch number suffix (last _XXX before .root)
    basename_output=$(basename "$first_file" | sed 's/_[0-9]\{3\}\.root$/.root/')
    output_file="$target_dir/$basename_output"
    
    if [ "$dry_run" = true ]; then
        # Count input files
        input_count=$(echo ${groups[$run]} | wc -w)
        echo ""
        echo "Run: $run"
        echo "  Output: $basename_output"
        echo "  Input files: $input_count"
        echo "  Command: hadd \"$output_file\" ${groups[$run]}"
    else
        # Run hadd with output file followed by all input files
        if hadd "$output_file" ${groups[$run]}; then
            # Verify the output file exists and has non-zero size
            if [ -f "$output_file" ] && [ -s "$output_file" ]; then
                echo "Merged run $run into $output_file"
                successfully_merged_runs[$run]=1
            else
                echo "ERROR: Merge for run $run failed - output file is missing or empty"
            fi
        else
            echo "ERROR: hadd command failed for run $run"
        fi
    fi
done

if [ "$dry_run" = true ]; then
    echo ""
    echo "=========================================================="
    echo "Total groups to merge: ${#groups[@]}"
    exit 0
fi

# Print summary of merge operations
total_runs=${#groups[@]}
successful_runs=${#successfully_merged_runs[@]}
failed_runs=$((total_runs - successful_runs))

echo ""
echo "=========================================================="
echo "Merge Summary:"
echo "  Total runs: $total_runs"
echo "  Successful: $successful_runs"
echo "  Failed: $failed_runs"
echo "=========================================================="

# Clean up pre-merged files if requested
if [ "$clean_after" = true ]; then
    echo ""
    if [ "$dry_run" = true ]; then
        echo "DRY RUN: Would also clean up pre-merged files after merging"
    else
        # Only clean up batch files for successfully merged runs
        if [ ${#successfully_merged_runs[@]} -gt 0 ]; then
            cleanup_premerged_files "$target_dir" "${!successfully_merged_runs[@]}"
        else
            echo "WARNING: No runs were successfully merged, skipping cleanup"
        fi
    fi
fi

# Exit with error if any merges failed
if [ $failed_runs -gt 0 ]; then
    echo ""
    echo "WARNING: Some merges failed. Batch files for failed runs were NOT deleted."
    exit 1
fi
