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

## Quick Start

```bash
# 1. Setup environment
source /group/halld/Software/gluex_env.sh
gxenv

# 2. Create configuration
cp GeneratorConfigExamples/RBHG.config my_study.config
# Edit my_study.config for your physics parameters

# 3. Run generator phase
python swif2_gen.py my_study.config

# 4. Prepare MC simulation
python prepareSimulation.py <study_directory>

# 5. Run DSelector analysis
python swif2_RBHG_DSelector.py <study_directory>
```

## Pipeline Phases

1. **Generator** - Event production (`swif2_gen.py`)
2. **MCWrapper** - GEANT4 simulation (`prepareSimulation.py`)
3. **DSelector** - Physics analysis (`swif2_RBHG_DSelector.py`)
4. **TMVA** - ML selection (future)
5. **Histograms** - Final results (future)

## Directory Structure

```
swif2-generator-framework/
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
    └── JANAConfigs/
```

## Requirements

- Python 3.6+
- NumPy
- HDDM Python bindings
- Fortran compiler (gfortran)
- JLab computing farm access
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
