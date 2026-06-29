# DSelector Templates and Configuration

This directory contains DSelector files and configuration templates for the swif2-generator-framework DSelector analysis stage.

## Contents

- **DSelector_2eMissingProton_Systematics.C/h** - Default DSelector for e+e- missing proton analysis
- **jobs_root_analysis.config** - Template configuration file for DSelector farm jobs  
- **Run_Selector.C** - ROOT macro for running DSelector analysis

## Usage

### Single Configuration Processing

Process a single individual configuration:

```bash
python3 swif2_RBHG_DSelector.py \
    /path/to/output/FFS1/BC_FFN/1801_0DEG_FFN_DBLRAD/rbhg_config.json \
    --dry-run
```

### Hierarchical Batch Processing

Process all configurations in a run period master config (like `prepareSimulation.py`):

```bash
python3 swif2_RBHG_DSelector.py \
    /path/to/output/FFS1/BC_FFN/run_period_master_config.json \
    --dry-run
```

Process all datasets in an FFS master config:

```bash
python3 swif2_RBHG_DSelector.py \
    /path/to/output/FFS1/ffs_master_config.json \
    --dry-run
```

### Using Custom DSelector

Specify a different DSelector file:

```bash
python3 swif2_RBHG_DSelector.py config.json \
    --dselector /path/to/custom/DSelector.C
```

## Features

### Automatic Configuration

The framework automatically:

1. **Detects tree names** - Opens a sample ROOT file to determine the tree name
2. **Calculates time limits** - Estimates job time based on average file size:
   - < 10 MB: 30 minutes
   - 10-30 MB: 1 hour
   - 30-60 MB: 2 hours
   - > 60 MB: 4 hours
3. **Resolves paths** - Finds actual volatile directory paths even if JSON configs have incorrect paths
4. **Generates workflow names** - Creates unique workflow names from study parameters

### Path Resolution

The script handles the discrepancy between JSON config paths and actual volatile paths:

- **JSON may say**: `/volatile/.../FFS1/BC_FFN/0DEG_DBLRAD/simulation/trees`  
- **Actual path is**: `/volatile/.../FFS1/BCFFN/1801_0DEG_FFN_DBLRAD/Simulation/root/trees`

The framework automatically constructs and finds the correct paths.

## Output

DSelector output is placed in:
```
/volatile/halld/home/acschick/RBHG/{STUDY}/{NAMETAG+FF}/{CONFIG}/DSelector/
```

For example:
```
/volatile/halld/home/acschick/RBHG/FFS1/BCFFN/1801_0DEG_FFN_DBLRAD/DSelector/
```

## DSelector Filename Handling

### Automatic (Default)
By default, the framework modifies DSelector output filenames to include unique identifiers:
- `Hist_{study}_{nametag}_{ff}_{period}_{pol}_{lepton}.root`
- `Flat_{study}_{nametag}_{ff}_{period}_{pol}_{lepton}.root`

### Tokens Mode (DSelector_2eMissingProton_Systematics)
The included DSelector has a "tokens" mode that automatically generates output names from input filenames. Set `dOutputFileName = "tokens"` in the DSelector to enable this.

## Configuration Template

See `jobs_root_analysis.config` for the template that gets auto-generated with appropriate values for:
- Job resources (cores, RAM, disk, time)
- Environment (2026launch infrastructure)
- Workflow name
- Input/output paths
- Tree name
- DSelector path

## Infrastructure

Uses the 2026launch batch submission infrastructure:
- Launch script: `/work/halld/home/acschick/channels/batch_submission/2026launch/launch_scripts/launch/launch.py`
- Environment: `version_7.5.0.xml`
- Script: `script.sh`
- Run_Selector: `Run_Selector.C`
