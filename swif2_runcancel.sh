#!/bin/bash

# Parse arguments
workflow_prefix="$1"
action="$2"
shift 2  # Remove first two arguments, leaving any additional args

# Check if minimum arguments are provided
if [ -z "$workflow_prefix" ] || [ -z "$action" ]; then
    echo "Usage: $0 <workflow_prefix> <action> [additional_args...]"
    echo ""
    echo "Actions:"
    echo "  run              - Run all workflows matching prefix"
    echo "  cancel           - Cancel and delete all workflows matching prefix"
    echo "  status           - Show status of all workflows matching prefix"
    echo "  list             - List all workflows matching prefix"
    echo "  modify-jobs      - Modify jobs in all workflows matching prefix"
    echo "  retry-jobs       - Retry jobs in all workflows matching prefix"
    echo ""
    echo "For modify-jobs, pass additional swif2 modify-jobs arguments:"
    echo "  Example: $0 SIM_FFS1 modify-jobs -time add 8h -problems SLURM_TIMEOUT"
    echo "  Example: $0 SIM_FFS1 modify-jobs -ram add 2GB -jobs PROBLEM"
    echo "  Example: $0 SIM_FFS1 modify-jobs -disk set 15GB -jobs all"
    echo ""
    echo "For retry-jobs, pass additional swif2 retry-jobs arguments:"
    echo "  Example: $0 SIM_FFS1 retry-jobs -problems SLURM_TIMEOUT"
    echo "  Example: $0 SIM_FFS1 retry-jobs -problems SLURM_OUT_OF_MEMORY"
    echo "  Example: $0 SIM_FFS1 retry-jobs -jobs PROBLEM"
    exit 1
fi

# Validate the action argument
if [[ "$action" != "run" && "$action" != "cancel" && "$action" != "status" && "$action" != "list" && "$action" != "modify-jobs" && "$action" != "retry-jobs" ]]; then
    echo "Invalid action: $action"
    echo "Action must be 'run', 'cancel', 'status', 'list', 'modify-jobs', or 'retry-jobs'."
    exit 1
fi

# For modify-jobs or retry-jobs, ensure additional arguments are provided
if [ "$action" == "modify-jobs" ] && [ "$#" -eq 0 ]; then
    echo "Error: modify-jobs requires additional arguments"
    echo "Example: $0 SIM_FFS1 modify-jobs -time add 8h -problems SLURM_TIMEOUT"
    exit 1
fi

if [ "$action" == "retry-jobs" ] && [ "$#" -eq 0 ]; then
    echo "Error: retry-jobs requires additional arguments"
    echo "Example: $0 SIM_FFS1 retry-jobs -problems SLURM_TIMEOUT"
    exit 1
fi

# Perform the action based on the input
if [ "$action" == "run" ]; then
    swif2 list | awk -v prefix="$workflow_prefix" '$0 ~ "^"prefix {print $1}' | while read wf; do
        echo "Running workflow: $wf"
        swif2 run -workflow "$wf"
    done
elif [ "$action" == "cancel" ]; then
    swif2 list | awk -v prefix="$workflow_prefix" '$0 ~ "^"prefix {print $1}' | while read wf; do
        echo "Canceling workflow: $wf"
        swif2 cancel -delete -workflow "$wf"
    done
elif [ "$action" == "status" ]; then
    swif2 list | awk -v prefix="$workflow_prefix" '$0 ~ "^"prefix {print $1}' | while read wf; do
        echo "=== Status for workflow: $wf ==="
        swif2 status -workflow "$wf"
        echo ""
    done
elif [ "$action" == "list" ]; then
    echo "Workflows matching prefix '$workflow_prefix':"
    swif2 list | awk -v prefix="$workflow_prefix" '$0 ~ "^"prefix {print $1}';
elif [ "$action" == "modify-jobs" ]; then
    # Get list of matching workflows
    matching_workflows=$(swif2 list | awk -v prefix="$workflow_prefix" '$0 ~ "^"prefix {print $1}')
    
    if [ -z "$matching_workflows" ]; then
        echo "No workflows found matching prefix: $workflow_prefix"
        exit 1
    fi
    
    echo "Found workflows matching prefix '$workflow_prefix':"
    echo "$matching_workflows"
    echo ""
    echo "Will run: swif2 modify-jobs -workflow <workflow> $@"
    echo ""
    read -p "Continue? (y/n) " -n 1 -r
    echo ""
    
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 0
    fi
    
    echo "$matching_workflows" | while read wf; do
        echo "Modifying jobs in workflow: $wf"
        swif2 modify-jobs -workflow "$wf" "$@"
        if [ $? -eq 0 ]; then
            echo "  ✓ Successfully modified $wf"
        else
            echo "  ✗ Failed to modify $wf"
        fi
        echo ""
    done
    echo "Done modifying jobs for all workflows."
elif [ "$action" == "retry-jobs" ]; then
    # Get list of matching workflows
    matching_workflows=$(swif2 list | awk -v prefix="$workflow_prefix" '$0 ~ "^"prefix {print $1}')
    
    if [ -z "$matching_workflows" ]; then
        echo "No workflows found matching prefix: $workflow_prefix"
        exit 1
    fi
    
    echo "Found workflows matching prefix '$workflow_prefix':"
    echo "$matching_workflows"
    echo ""
    echo "Will run: swif2 retry-jobs -workflow <workflow> $@"
    echo ""
    read -p "Continue? (y/n) " -n 1 -r
    echo ""
    
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 0
    fi
    
    echo "$matching_workflows" | while read wf; do
        echo "Retrying jobs in workflow: $wf"
        swif2 retry-jobs -workflow "$wf" "$@"
        if [ $? -eq 0 ]; then
            echo "  ✓ Successfully retried jobs in $wf"
        else
            echo "  ✗ Failed to retry jobs in $wf"
        fi
        echo ""
    done
    echo "Done retrying jobs for all workflows."
fi
