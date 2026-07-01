# DSelector Workflow Tools

This repository has two DSelector submission wrappers:

- `swif2_Data_DSelector.py`: runs a DSelector over real data ROOT trees.
- `swif2_RBHG_DSelector.py`: runs a DSelector over RBHG simulation trees produced by this framework.

They both submit SWIF2 workflows, but they take different kinds of configuration files.

Run the framework setup once after cloning:

```bash
./setup_framework.py
```

Setup detects `$USER`, writes `framework_settings.json`, creates safe user-owned base directories when possible, and creates `workflows_root_analysis.config` if it does not already exist. The scripts also have `$USER`-aware defaults, so they no longer depend on hardcoded personal paths.

## Real Data: `swif2_Data_DSelector.py`

Use this when you already have real data analysis trees, usually in `/mss/halld/...` or another shared tree directory, and want to submit one DSelector workflow per run-period/polarization/target combination.

This tool does accept `--config`:

```bash
python3 swif2_Data_DSelector.py --config my_dselector_data.config --dry-run
```

The config file is a simple `KEY=VALUE` file. It is not JSON.

### Minimal Config Example

```ini
RUN_PERIODS=2205
POLARIZATIONS=45DEG,135DEG,AMO
TARGETS=FULL,EMPTY

TREE_TYPE=gpb208_epemmisspb208__B2
TREE_NAME=gpb208_epemmisspb208__B2_Tree
ANALYSIS_VERSION_2205=ver01

DSELECTOR_PATH=/Users/andrew/Documents/swif2-generator-framework/DSelector_Templates/2eMissingProton_Systematics/DSelector_2eMissingProton_Systematics.C

OUTPUT_BASE_DIR=/volatile/halld/home/{USERNAME}/DSelectorData/{RunPeriod}
LOG_BASE_DIR=/farm_out/{USERNAME}/DSelectorData/{RunPeriod}
WORKFLOW_PREFIX=DSELECTOR_DATA
EXISTING_OUTPUT_MODE=fail

TIMELIMIT=8hrs
CORES=1
RAM=3GB
DISK=100MB
SWIF2_ACCOUNT=halld
SWIF2_PARTITION=production
```

`RUN_PERIODS`, `POLARIZATIONS`, and `TARGETS` may also be overridden from the command line:

```bash
python3 swif2_Data_DSelector.py \
  --config my_dselector_data.config \
  --run-period 2205 \
  --polarizations 45DEG,135DEG \
  --targets FULL,EMPTY \
  --dry-run
```

Remove `--dry-run` to submit.

### Output Directories And Existing Output

The data tool builds output directories from `OUTPUT_BASE_DIR`, then appends the polarization and, when active, the target type.

The default base paths are:

```ini
OUTPUT_BASE_DIR=/volatile/halld/home/{USERNAME}/RealDataAnalysis/{RunPeriod}
LOG_BASE_DIR=/farm_out/{USERNAME}/DSelector_logs/RealData/{RunPeriod}
```

With those defaults, a GlueX run over `1801` and `0DEG` writes under:

```text
/volatile/halld/home/$USER/RealDataAnalysis/1801/0DEG/
```

For CPP with target splitting, `2205`, `45DEG`, and `FULL` writes under:

```text
/volatile/halld/home/$USER/RealDataAnalysis/2205/45DEG/FULL/
```

You can place the placeholders directly in the config if you want a different layout:

```ini
OUTPUT_BASE_DIR=/volatile/halld/home/{USERNAME}/RealDataAnalysis/{RunPeriod}/{Polarization}/{TargetType}
LOG_BASE_DIR=/farm_out/{USERNAME}/DSelector_logs/RealData/{RunPeriod}/{Polarization}/{TargetType}
```

If `{Polarization}` is not present in `OUTPUT_BASE_DIR`, the tool appends the polarization directory automatically. If `{TargetType}` is not present and target splitting is active, it appends the target directory automatically.

The tool now checks for a non-empty output directory before writing the generated `jobs_root_analysis.config` or submitting jobs. Configure this with:

```ini
EXISTING_OUTPUT_MODE=fail
```

Valid modes are:

| Mode | Behavior |
| --- | --- |
| `fail` | Default. Do not submit that combination if the output directory already has contents. |
| `skip` | Skip that combination and continue with the rest. |
| `allow` | Reuse the existing directory. The generated config may be overwritten, and matching DSelector output filenames may be overwritten by jobs. |

You can override the mode on the command line:

```bash
python3 swif2_Data_DSelector.py \
  --config my_dselector_data.config \
  --existing-output skip \
  --dry-run
```

The tool does not delete old output for you. If you want a clean rerun after changing a DSelector, the safest pattern is to change `OUTPUT_BASE_DIR` or `WORKFLOW_PREFIX`, for example:

```ini
OUTPUT_BASE_DIR=/volatile/halld/home/{USERNAME}/RealDataAnalysis_dselectorV2/{RunPeriod}
WORKFLOW_PREFIX=DSELECTOR_DATA_V2
```

### `RunPeriods.json` Requirements

The data tool reads run metadata from `RunPeriods.json`. Each run period needs:

- `run_range`: numeric run range scanned by `launch.py`.
- `RCDB_query`: query template used to select runs.
- `Polarizations`: allowed polarization labels.
- `tape_analysisTrees_paths.DATA`: input tree directories.
- `TargetTypes`: optional; use this for CPP FULL/EMPTY splitting.

Example CPP-style production layout:

```json
"2205": {
  "description": "2022 CPP data",
  "experiment": "CPP/NPP",
  "run_range": "100500-101622",
  "RCDB_query": "@is_cpp_production and (status==1 or status==6) and (polarization_angle=={polDeg}) and target_type==\"{targetType}\"",
  "TargetTypes": ["FULL", "EMPTY"],
  "tape_analysisTrees_paths": {
    "DATA": {
      "ver01": {
        "tree_gpb208_epemmisspb208__B2": "/mss/halld/RunPeriod-2022-05/analysis/ver01/tree_gpb208_epemmisspb208__B2/merged"
      }
    }
  },
  "Polarizations": {
    "135DEG": {},
    "45DEG": {},
    "AMO": {}
  }
}
```

For GlueX running, omit `TargetTypes` and use a GlueX RCDB query, for example:

```json
"RCDB_query": "@is_2018production and @status_approved and (polarization_angle=={polDeg})"
```

### How Selection Works

`POLARIZATIONS=ALL` means the polarizations listed under that run period in `RunPeriods.json`. For `2205`, that currently means `135DEG`, `45DEG`, and `AMO`.

AMO is translated to:

```text
polarization_angle==-1.0
```

CPP target splitting uses `{targetType}` in the RCDB query. If `TargetTypes` is `["FULL", "EMPTY"]`, the tool submits separate workflows for each target type:

```text
DSELECTOR_DATA_2205_45DEG_FULL
DSELECTOR_DATA_2205_45DEG_EMPTY
```

For special skim directories that are separated by polarization, the tool can infer a matching analysis path key such as `epemMLskim_45DEG`. For normal future production where all polarizations live under one `ver01` tree directory, use `ANALYSIS_VERSION_2205=ver01` or let it fall back to `ver01`.

### Fresh ROOT Trees

If you have freshly made ROOT trees, point `RunPeriods.json` at their directory under `tape_analysisTrees_paths.DATA`. The filenames must include six-digit run numbers because `DSelector_Templates/launch.py` finds files by run number, for example:

```text
tree_gpb208_epemmisspb208__B2_101582_000.root
tree_gpb208_epemmisspb208__B2_101582.root
```

Then dry-run first:

```bash
python3 swif2_Data_DSelector.py --config my_dselector_data.config --dry-run
```

## Simulation: `swif2_RBHG_DSelector.py`

Use this for RBHG simulation output generated by the framework. This tool does not use `--config`; it takes a JSON config path as a required positional argument:

```bash
python3 swif2_RBHG_DSelector.py path/to/config.json --dry-run
```

It understands three generated config levels:

- Individual config: one generated dataset directory, usually `rbhg_config.json`.
- Run-period master config: a group of individual configs.
- FFS master config: a full form-factor study with several datasets.

Examples:

```bash
# One simulated dataset
python3 swif2_RBHG_DSelector.py \
  output/FFS1/qDATAq_FFN/1801_0DEG_FFN_DBLRAD/rbhg_config.json \
  --dry-run

# One run-period master
python3 swif2_RBHG_DSelector.py \
  output/FFS1/qDATAq_FFN/run_period_master_config.json \
  --dry-run

# Entire FFS study
python3 swif2_RBHG_DSelector.py \
  output/FFS1/ffs_master_config.json \
  --dry-run
```

Remove `--dry-run` to submit.

### Simulation Options

Use a different DSelector:

```bash
python3 swif2_RBHG_DSelector.py output/FFS1/ffs_master_config.json \
  --dselector /path/to/MyDSelector.C \
  --dry-run
```

Copy the DSelector into each output directory and rewrite output filenames:

```bash
python3 swif2_RBHG_DSelector.py output/FFS1/ffs_master_config.json \
  --modify-output-names
```

By default the tool uses the central DSelector at:

```text
DSelector_Templates/2eMissingProton_Systematics/DSelector_2eMissingProton_Systematics.C
```

It auto-detects the TTree name from the first `tree_*.root` file when ROOT Python is available. If detection fails, it falls back to:

```text
epemmissprot__B2_Tree
```

### Important Difference From Data

`swif2_RBHG_DSelector.py` discovers input trees from the generated simulation JSON and expected output directory structure. It does not use `RunPeriods.json` for RCDB filtering, and it does not split FULL/EMPTY targets. That split belongs to the real-data tool.

The simulation wrapper writes an `RBHG_dselector.config` and submits it through the root-analysis launch script. It is meant to be run after simulation output trees already exist.

## Quick Decision Table

| Situation | Tool | Config Type |
| --- | --- | --- |
| Real GlueX/CPP data trees | `swif2_Data_DSelector.py` | `KEY=VALUE` config plus `RunPeriods.json` |
| Fresh real data ROOT trees in a custom directory | `swif2_Data_DSelector.py` | Add/update `RunPeriods.json`, then use `--config` |
| RBHG simulation output from this framework | `swif2_RBHG_DSelector.py` | Generated JSON config |
| Entire form-factor simulation study | `swif2_RBHG_DSelector.py` | `ffs_master_config.json` |
