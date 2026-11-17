# SWIF2 Generator Framework

Automated workflow system for exotic event generators in the GlueX experiment at Jefferson Lab.

## Overview

This framework provides a unified swif2-based automation system for event generators that don't fit cleanly into the standard MCWrapper pipeline. Features include:

- **Automated job submission** via SWIF2
- **JSON-based chronicling** of all workflow parameters
- **Event validation** and quality checks  
- **Integrated MC simulation** chain

## Current Generators

### RBHG (Rory's Bethe-Heitler Generator)
Electron-positron pair production via Bethe-Heitler process for exotic meson studies.

### SPIZG (Sergey's π⁰ Generator)
Single π⁰ photoproduction generator for background studies and cross-section measurements.

**Future generators welcome!** This framework is designed to be extensible.

## One-time framework setup

After cloning, run this from the framework root on your JLab Linux host:

```bash
python3 setup_framework.py
```

Setup detects your username, framework location, and shell (prompting if needed).
It creates `output/RBHG` and `output/SPIZG`, attempts to create the configured
`/volatile` output and real-data DSelector log directories, updates the example
configs with your username, normalizes cobrems paths in `RunPeriods.json`, creates
`data_dselector.config` if absent, and saves `framework_settings.json` for
downstream scripts. No separate `mkdir` is needed after successful setup.
Check any directory-creation warnings: setup continues if a directory cannot be
created. Generator-specific study and log directories are created when jobs run.

Setup configures the framework; it does not install dependencies or load the
GlueX environment. Load that environment in each new shell as shown below.

## Quick Start: SWIF2 submission

```bash
# 1. Load environment (after one-time setup)
source /group/halld/Software/gluex_env.sh
gxenv

# 2. Copy and customize configuration
cp GeneratorConfigExamples/RBHG.config my_study.config
# Edit my_study.config for your physics parameters

# 3. Run generator phase
python swif2_gen.py my_study.config

# 4. Prepare MC simulation
python prepareSimulation.py <study_directory>

# 5. Run DSelector analysis
python swif2_RBHG_DSelector.py <study_directory>
```

## Interactive/local run: 10,000 events on one core

After the one-time setup, from the framework root on a JLab Linux host:

```bash
source /group/halld/Software/gluex_env.sh
gxenv
python swif2_gen.py GeneratorConfigExamples/RBHG_10k_local.config --local
```

`--interactive` is an equivalent alias; use either command:

```bash
python swif2_gen.py GeneratorConfigExamples/RBHG_10k_local.config --interactive
```

The [10k example](GeneratorConfigExamples/RBHG_10k_local.config) sets
`FIXED_EVENT_COUNT True`, `NEVENTS_TOTAL 10000`, and `NEVENTS_PERFILE 10000`.
It leaves `RUN_PERIOD` blank to select only `POL_DEG 90`, with `FLUX_STUDY`
and `TENX_FLAG` disabled. Physics settings otherwise match `RBHG.config`:
GlueX electron-positron production with internal radiation enabled.

This compiles RBHG and runs one serial generator process, followed by ASCII-to-HDDM
conversion. No core-count option is needed. Counts are configured in the file;
there is no `--nevents` option. Omitting `--local`/`--interactive` submits to SWIF2,
even when the config filename contains `local`. MC simulation and DSelector
analysis are separate steps and are not started by this command.

Results are written below `output/RBHG/interactive_10k/`, including `vectors/`,
`hists/`, and `FortranFiles/` within the generated study structure. For a fresh
rerun, copy the config and change `STUDY_NAME`; existing outputs may be skipped
or existing directories rejected.

The current local runner still creates directories under
`/farm_out/<USERNAME>/RBHG/`, and requires GNU `sed`, `gfortran`, and a `python`
environment with NumPy and `hddm_s`. It uses the current shell environment;
`ENV_FILE` does not initialize local execution. `USERNAME {USER}` uses your login
name. RBHG currently fixes generation output and log locations in code;
`GENERATOR_OUTPUT_DIR_BASE` and `LOG_DIR_BASE` do not override them.

## Pipeline Phases

1. **Generator** - Event production (`swif2_gen.py`)
2. **MCWrapper** - GEANT4 simulation (`prepareSimulation.py`)
3. **DSelector** - Physics analysis (`swif2_RBHG_DSelector.py`)
4. **TMVA** - ML selection (future)
5. **Histograms** - Final results (future)

## Directory Structure

```
swif2-generator-framework/
├── setup_framework.py            # One-time directories and user configuration
├── swif2_gen.py                  # Main dispatcher
├── swif2_common.py               # Shared utilities
├── prepareSimulation.py          # MC setup
├── RunPeriods.json               # GlueX run definitions
├── GeneratorConfigExamples/      # Example configs with docs
├── hddm_scripts/                 # Event file utilities
├── generators/                   # Generator-specific code
│   ├── RBHG/
│   ├── SPIZG/
│   └── CobremsDistributions/
└── template_MCWrapper/           # MC simulation templates
    └── jana_configs/
```

## Requirements

- Python 3.7+ (`python` must select Python 3)
- NumPy
- HDDM Python bindings
- Fortran compiler (gfortran)
- JLab Linux environment with GNU `sed` and writable framework and `/farm_out` directories
- SWIF2 for batch submission (not needed for `--local`)
- GlueX software stack

## Documentation

See `docs/` for detailed guides and `GeneratorConfigExamples/README.md` for config help.

## Contributing

To add a new generator:
1. Create `generators/YOURGEN/` directory
2. Implement `swif2_YOURGEN.py` (follow RBHG pattern)
3. Add templates in `generators/YOURGEN/template_YOURGEN/`
4. Update `swif2_gen.py` dispatcher

## Contact

A. Schick (acschick@jlab.org)
