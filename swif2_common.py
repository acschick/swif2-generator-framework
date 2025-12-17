#!/usr/bin/env python3

"""
Reusable SWIF2 workflow management module
Extracted from swif2_RBHG.py for use across the analysis pipeline
"""

import subprocess
import os
from datetime import datetime

def create_swif2_workflow(workflow_name):
    """
    Create a SWIF2 workflow
    
    Args:
        workflow_name (str): Name of the workflow to create
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        print(f"Creating SWIF2 workflow: {workflow_name}")
        subprocess.run(
            ["swif2", "create", "-workflow", workflow_name], 
            capture_output=True, text=True, check=True
        )
        print(f"Workflow created successfully: {workflow_name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error creating workflow {workflow_name}: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        return False

def add_swif2_job(workflow_name, job_name, command, **job_options):
    """
    Add a job to a SWIF2 workflow
    
    Args:
        workflow_name (str): Name of the workflow
        job_name (str): Name of the job
        command (str): Shell command to execute
        **job_options: Additional job options (cores, memory, disk, etc.)
        
    Returns:
        bool: True if successful, False otherwise
    """
    # Default job options
    defaults = {
        'cores': 1,
        'memory': '4GB',
        'disk': '10GB',
        'time': '2h'
    }
    defaults.update(job_options)
    
    # Build the SWIF2 command
    swif_cmd = [
        "swif2", "add-job", 
        "-workflow", workflow_name,
        "-name", job_name,
        "-cores", str(defaults['cores']),
        "-ram", defaults['memory'],
        "-disk", defaults['disk'],
        "-time", defaults['time']
    ]
    
    # Add optional parameters
    if 'stdout' in defaults:
        swif_cmd.extend(["-stdout", defaults['stdout']])
    if 'stderr' in defaults:
        swif_cmd.extend(["-stderr", defaults['stderr']])
    if 'input' in defaults:
        for input_file in defaults['input']:
            swif_cmd.extend(["-input", input_file])
    if 'output' in defaults:
        for output_file in defaults['output']:
            swif_cmd.extend(["-output", output_file])
    
    # Add the command
    swif_cmd.append(command)
    
    try:
        print(f"Adding job: {job_name} to workflow: {workflow_name}")
        subprocess.run(swif_cmd, capture_output=True, text=True, check=True)
        print(f"Job added successfully: {job_name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error adding job {job_name}: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        return False

def submit_workflow(workflow_name):
    """
    Submit a SWIF2 workflow for execution
    
    Args:
        workflow_name (str): Name of the workflow to submit
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        print(f"Submitting workflow: {workflow_name}")
        subprocess.run(
            ["swif2", "run", "-workflow", workflow_name], 
            capture_output=True, text=True, check=True
        )
        print(f"Workflow submitted successfully: {workflow_name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error submitting workflow {workflow_name}: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        return False

def get_workflow_status(workflow_name):
    """
    Get the status of a SWIF2 workflow
    
    Args:
        workflow_name (str): Name of the workflow
        
    Returns:
        dict: Workflow status information
    """
    try:
        result = subprocess.run(
            ["swif2", "status", "-workflow", workflow_name], 
            capture_output=True, text=True, check=True
        )
        
        # Parse the status output
        status_info = {
            'workflow_name': workflow_name,
            'status': 'unknown',
            'jobs_total': 0,
            'jobs_pending': 0,
            'jobs_running': 0,
            'jobs_successful': 0,
            'jobs_failed': 0,
        }
        
        # Simple parsing - could be enhanced based on actual SWIF2 output format
        lines = result.stdout.split('\n')
        for line in lines:
            if 'STATUS' in line.upper():
                parts = line.split()
                if len(parts) >= 2:
                    status_info['status'] = parts[1].lower()
            elif 'TOTAL' in line.upper():
                # Extract job counts if available
                pass
        
        return status_info
        
    except subprocess.CalledProcessError as e:
        print(f"Error getting status for workflow {workflow_name}: {e}")
        return {'workflow_name': workflow_name, 'status': 'error'}

def cancel_workflow(workflow_name):
    """
    Cancel a SWIF2 workflow
    
    Args:
        workflow_name (str): Name of the workflow to cancel
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        print(f"Cancelling workflow: {workflow_name}")
        subprocess.run(
            ["swif2", "cancel", "-workflow", workflow_name], 
            capture_output=True, text=True, check=True
        )
        print(f"Workflow cancelled successfully: {workflow_name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error cancelling workflow {workflow_name}: {e}")
        return False

def generate_workflow_name(step, study_name, dataset_name="", timestamp=True):
    """
    Generate a standardized workflow name
    
    Args:
        step (str): Pipeline step (gen, prep, sim, dsel, tmva, hist)
        study_name (str): Study name
        dataset_name (str): Dataset name (optional)
        timestamp (bool): Include timestamp for uniqueness
        
    Returns:
        str: Generated workflow name
    """
    parts = [step, study_name]
    if dataset_name:
        parts.append(dataset_name.replace('/', '_'))
    
    if timestamp:
        # Add timestamp for uniqueness
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        parts.append(ts)
    
    workflow_name = "_".join(parts)
    # Ensure workflow name is valid (no special characters)
    workflow_name = workflow_name.replace(' ', '_').replace('-', '_')
    
    return workflow_name

def create_hddm_merge_job(workflow_name, input_dir, output_file, job_name=None):
    """
    Create a SWIF2 job to merge HDDM files using hddm_merge_files
    
    Args:
        workflow_name (str): SWIF2 workflow name
        input_dir (str): Directory containing vectors_*.hddm files
        output_file (str): Output merged HDDM file path
        job_name (str): Optional job name
        
    Returns:
        bool: True if job created successfully
    """
    if job_name is None:
        job_name = f"merge_hddm_{os.path.basename(input_dir)}"
    
    # Build the merge command
    # First, we need to set up environment and find all HDDM files
    command = f"""#!/bin/bash
cd {input_dir}
source /group/halld/builds/Linux_CentOS7-x86_64-gcc4.8.5/setup.csh || source /group/halld/setup_jlab.csh
gxenv

# Find all vectors_*.hddm files
hddm_files=$(ls vectors_*.hddm 2>/dev/null)
if [ -z "$hddm_files" ]; then
    echo "ERROR: No vectors_*.hddm files found in {input_dir}"
    exit 1
fi

echo "Found HDDM files: $hddm_files"

# Count events before merging
total_events_before=0
for file in $hddm_files; do
    events=$(python {os.path.dirname(__file__)}/../template_MCWrapper/count_hddm_events.py $file | grep "Total events:" | cut -d: -f2 | tr -d ' ')
    total_events_before=$((total_events_before + events))
    echo "File $file: $events events"
done
echo "Total events before merge: $total_events_before"

# Merge the files
echo "Merging HDDM files into {output_file}"
hddm_merge_files -o {output_file} $hddm_files

# Verify the merge
if [ -f "{output_file}" ]; then
    echo "Merge completed. Verifying event count..."
    events_after=$(python {os.path.dirname(__file__)}/../template_MCWrapper/count_hddm_events.py {output_file} | grep "Total events:" | cut -d: -f2 | tr -d ' ')
    echo "Total events after merge: $events_after"
    
    if [ "$total_events_before" -eq "$events_after" ]; then
        echo "SUCCESS: Event count verified ($events_after events)"
        echo "$events_after" > {output_file}.event_count
    else
        echo "WARNING: Event count mismatch! Before: $total_events_before, After: $events_after"
        exit 1
    fi
else
    echo "ERROR: Merge failed - output file not created"
    exit 1
fi
"""
    
    # Job options
    job_options = {
        'cores': 1,
        'memory': '8GB',
        'disk': '20GB',
        'time': '1h',
        'stdout': f"{input_dir}/{job_name}.out",
        'stderr': f"{input_dir}/{job_name}.err",
        'output': [output_file, f"{output_file}.event_count"]
    }
    
    return add_swif2_job(workflow_name, job_name, command, **job_options)

# Workflow management utilities
class SwifWorkflowManager:
    """
    High-level workflow management class
    """
    
    def __init__(self, base_name, user="acschick"):
        self.base_name = base_name
        self.user = user
        self.workflows = {}
    
    def create_step_workflow(self, step, study_name, dataset_name=""):
        """Create a workflow for a specific pipeline step"""
        workflow_name = generate_workflow_name(step, study_name, dataset_name)
        
        if create_swif2_workflow(workflow_name):
            self.workflows[step] = workflow_name
            return workflow_name
        return None
    
    def submit_all(self):
        """Submit all created workflows"""
        for step, workflow_name in self.workflows.items():
            print(f"Submitting {step} workflow: {workflow_name}")
            submit_workflow(workflow_name)
    
    def get_all_status(self):
        """Get status of all workflows"""
        status = {}
        for step, workflow_name in self.workflows.items():
            status[step] = get_workflow_status(workflow_name)
        return status

if __name__ == "__main__":
    # Test the module
    print("Testing SWIF2 common module")
    
    # Test workflow name generation
    test_name = generate_workflow_name("prep", "TestStudy", "qDATAq_FFN")
    print(f"Generated workflow name: {test_name}")
    
    print("SWIF2 common module loaded successfully")