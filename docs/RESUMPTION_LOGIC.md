# Resumption Logic for Generator Framework

## Overview

The framework now supports **automatic resumption** when re-running the same configuration. This prevents duplicate work and allows you to safely restart after crashes or interruptions.

## How It Works

### Configuration Option

```
SKIP_EXISTING_OUTPUTS        True    # Skip dataset/polarization combos that already have output
```

- **`True` (default)**: Checks for existing output before processing each dataset/polarization combination
- **`False`**: Always regenerate everything (original behavior)

### What Gets Checked

For each dataset/polarization combination, the framework checks:

1. **Directory existence**: Does the output directory exist?
   - Example: `output/RBHG/CPP_FFS1/BC_FFN/2205_135DEG_FFN_DBLRAD/`

2. **Config file presence**: Does `rbhg_config.json` exist in that directory?
   - This file is created after successful setup and compilation

3. **Job submission status**: Does the config show `jobs_submitted: true`?
   - Added to config after ALL jobs successfully submitted to SWIF2
   - Prevents skipping datasets where setup completed but submission failed

### Behavior

When `SKIP_EXISTING_OUTPUTS = True`:

####  **Will Skip**
- Any dataset/polarization with existing output directory AND config file
- Prints clear message showing what was skipped and why
- Example output:
  ```
  ⏭  SKIPPING: BC_FFN / 2205_135DEG
     Reason: Found existing output: .../2205_135DEG_FFN_DBLRAD
     Config: .../2205_135DEG_FFN_DBLRAD/rbhg_config.json
     (Set SKIP_EXISTING_OUTPUTS = False to regenerate)
  ```

####  **Will Process**
- Dataset/polarizations with no output directory
- Directories that exist but have no config file (incomplete setup)
- New configurations not previously run

## Workflow Conflict Detection

**NEW**: The framework now checks if a SWIF2 workflow already exists before creating a new one.

### When Conflicts Occur

A workflow conflict happens when:
- You re-run the same configuration
- Previous run was interrupted but workflows still exist in SWIF2
- Jobs were already submitted and are running/completed

### Interactive Resolution

When a conflict is detected, you'll be prompted with 4 options:

```
WARNING: Workflow 'gen_CPP_FFS1_BCFFN_2205_135DEG_ee' already exists in SWIF2

Options:
  1) SKIP this dataset/polarization (recommended if jobs are running)
  2) Create NEW workflow with different name (appends _1, _2, etc.)
  3) CANCEL existing workflow and create fresh one (CAUTION: kills running jobs)
  4) Abort entire script
```

**Option 1: SKIP** (Recommended if jobs running)
- Leaves existing workflow untouched
- Skips to next dataset/polarization
- Safe for active jobs

**Option 2: NEW workflow** (Safe, creates parallel submission)
- Appends `_1`, `_2`, etc. to workflow name
- Both workflows can coexist
- Useful if you want to regenerate with different settings

**Option 3: CANCEL existing** (CAUTION)
- Cancels existing workflow (kills any running jobs!)
- Creates fresh workflow with original name
- Use only if previous run failed badly

**Option 4: ABORT**
- Exits script immediately
- Lets you manually inspect/cleanup

### Job Submission Tracking

After ALL jobs are successfully submitted, the config file is updated:

```json
{
  "submission_status": {
    "jobs_submitted": true,
    "workflow_name": "gen_CPP_FFS1_BCFFN_2205_135DEG_ee",
    "total_jobs": 8000,
    "submission_time": "2026-02-16T14:23:45.123456"
  },
  ...
}
```

This marker ensures:
-  Skip truly complete submissions (config exists + jobs_submitted = true)
-  Regenerate incomplete setups (config exists but jobs_submitted = false/missing)
-  Regenerate if submission failed midway

## Use Cases

### Scenario 1: Crash During FFS Study

You're running Form Factor Study (FFS) mode with 4 datasets × 3 polarizations = 12 combinations:
-  qDATAq_FFN completed (3 polarizations, all jobs submitted)
-  BC_FFN completed (3 polarizations, all jobs submitted) 
-  BC_FF1 crash after 1 polarization (setup complete, NO jobs submitted)
-  SIM_FF1 not started

**Solution**: Re-run with `SKIP_EXISTING_OUTPUTS = True`
- Skips 6 completed combinations (have jobs_submitted=true)
- **Regenerates BC_FF1/135DEG** (config exists but jobs_submitted=false)
- Completes 2 remaining BC_FF1 polarizations  
- Generates all 3 SIM_FF1 polarizations

### Scenario 2: Server Issues

SWIF2 server returns errors partway through submission:
- Some jobs submitted successfully
- Script crashes on server error
- Workflows may still exist in SWIF2

**Solution**: 
1. Fix applied to handle server errors gracefully
2. Re-run with `SKIP_EXISTING_OUTPUTS = True`
3. **Workflow check will detect existing workflows and prompt you**
   - Choose option 1 (SKIP) if jobs are running
   - Choose option 2 (NEW) to create parallel submission
4. Only processes remaining work

### Scenario 3: Accidental Re-run

You forgot you already submitted jobs and run the script again:

**What happens:**
- Skip check sees existing outputs with jobs_submitted=true → skips them ✓
- **Workflow check catches any conflicts**:
  ```
  WARNING: Workflow 'gen_CPP_FFS1_BCFFN_2205_135DEG_ee' already exists
  Options: [1-4]
  ```
- Choose option 1 to safely skip and keep existing jobs running

### Scenario 4: Intentional Regeneration

You modified settings and want to regenerate specific datasets:

**Option A**: Delete specific output directories
- Remove directories you want to regenerate
- Keep directories you want to preserve
- Re-run with `SKIP_EXISTING_OUTPUTS = True`

**Option B**: Regenerate everything
- Set `SKIP_EXISTING_OUTPUTS = False`
- Re-run (will regenerate all combinations)

## Technical Details

### Detection Logic

The `check_output_exists()` function:

```python
def check_output_exists(base_output_dir, study_name, nametag, form_factor, 
                       run_period, polarization):
    """
    Check if output already exists for a dataset/polarization combination.
    
    Returns:
        (exists: bool, config_path: str, reason: str)
    """
    # 1. Look for dataset directory: {study_name}/{nametag}_{form_factor}/
    # 2. Search for subdirectories matching: {run_period}_{polarization}_*
    # 3. Check for rbhg_config.json in matching directories
    # 4. Return existence status with explanation
```

### Integration Points

The check happens in the main loop:

```python
for dataset in ffs_datasets:
    for run_period, polarizations in actual_run_configs.items():
        for polarization, events in polarizations.items():
            # ← Resumption check happens here
            if SKIP_EXISTING_OUTPUTS:
                exists, config_path, reason = check_output_exists(...)
                if exists:
                    print(f"⏭️  SKIPPING: {dataset}/{run_period}_{polarization}")
                    continue  # Skip to next iteration
            
            # Proceed with directory creation, compilation, submission...
```

### Portability to Other Generators

This logic is designed to be easily adapted for SPIZG and future generators:

1. Add `SKIP_EXISTING_OUTPUTS` config option
2. Copy the `check_output_exists()` helper function
3. Modify to check for generator-specific markers (e.g., `spizg_config.json`)
4. Add check in main processing loop

## Limitations & Considerations

###  What's NOT Checked

- **Job completion status**: Skips based on setup completion, not job success
- **Workflow status**: Doesn't check if SWIF2 jobs finished successfully
- **Output file integrity**: Doesn't validate that generated files are correct

### 🔍 When to Manually Verify

After resuming, you may want to check:
- SWIF2 workflow status: `swif2 status <workflow_name>`
- Output file counts: Verify expected number of HDDM files
- Job logs: Check for runtime errors

###  Forcing Regeneration

To regenerate specific combinations:

```bash
# Option 1: Remove entire dataset
rm -rf output/RBHG/CPP_FFS1/BC_FF1/

# Option 2: Remove specific polarization
rm -rf output/RBHG/CPP_FFS1/BC_FF1/2205_135DEG_FF1_DBLRAD/

# Option 3: Disable skipping (regenerates everything)
# In config file: SKIP_EXISTING_OUTPUTS = False
```

## Examples

### Example 1: Fresh Run
```bash
# First run - nothing exists yet
python swif2_gen.py RBHG_CPP_FFS1.config

# Output: Processes all 12 combinations
# Result: All directories created, jobs submitted
```

### Example 2: Resuming After Crash
```bash
# Second run - some work already done
python swif2_gen.py RBHG_CPP_FFS1.config

# Output:
#   ⏭  SKIPPING: qDATAq_FFN / 2205_135DEG (7 more skips...)
#   Processing: BC_FF1 / 2205_45DEG
#   Processing: BC_FF1 / 2205_AMO
#   Processing: SIM_FF1 / 2205_135DEG (3 more new...)
```

### Example 3: Selective Regeneration
```bash
# Remove one dataset to regenerate
rm -rf output/RBHG/CPP_FFS1/BC_FFN/

# Re-run
python swif2_gen.py RBHG_CPP_FFS1.config

# Output:
#   ⏭  SKIPPING: qDATAq_FFN / ... (skips completed)
#   Processing: BC_FFN / ... (regenerates removed dataset)
#   ⏭  SKIPPING: BC_FF1 / ... (skips completed)
#   ⏭  SKIPPING: SIM_FF1 / ... (skips completed)
```

## Best Practices

1. **Keep `SKIP_EXISTING_OUTPUTS = True` by default** for safety
2. **Check terminal output** for skip messages to verify expected behavior
3. **Verify SWIF2 workflows** after resuming to ensure jobs are running
4. **Document intentional regenerations** (e.g., why you deleted certain outputs)
5. **Back up critical outputs** before regenerating with `SKIP_EXISTING_OUTPUTS = False`

## Troubleshooting

### "Not skipping despite existing output"

**Check:**
- Does `rbhg_config.json` exist in the output directory?
- Was setup interrupted before config file was written?

**Solution:** Remove incomplete directory and re-run

### "Skipping when I want to regenerate"

**Options:**
1. Set `SKIP_EXISTING_OUTPUTS = False` in config
2. Delete specific output directories manually
3. Move existing outputs to backup location

### "How do I check what will be skipped?"

**Preview mode** (currently not implemented, but could add):
- Could add `--dry-run` flag to show what would be processed/skipped
- Currently: Check output directories manually before running

## Future Enhancements

Potential improvements:
- `--dry-run` mode to preview what will be done
- `--force-regenerate` flag for specific datasets
- Check SWIF2 workflow completion status
- Validate output file integrity
- Resume from master config files (for multi-dataset coordination)
