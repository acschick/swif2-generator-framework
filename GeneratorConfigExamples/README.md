# Generator Configuration Examples

Example configuration files for event generators with detailed documentation.

## Quick Start

Run `python3 setup_framework.py` once from the framework root after cloning.
It creates the output directories and configures the examples for your username.
See [one-time setup](../README.md#one-time-framework-setup) for details.

1. **Copy a config file to your working directory:**
   ```bash
   cd /path/to/swif2-generator-framework
   cp GeneratorConfigExamples/RBHG.config my_study.config
   ```

2. **Edit the config:**
   - Set `USERNAME` to your JLab username
   - Leave `FRAMEWORK_HOME` blank (auto-detects framework location)
   - Or set it explicitly: `FRAMEWORK_HOME /path/to/swif2-generator-framework`

3. **Run the generator:**
   ```bash
   python swif2_gen.py my_study.config
   ```

## Path Configuration

### FRAMEWORK_HOME

This should point to the **framework root directory** (containing `generators/`, `hddm_scripts/`, `RunPeriods.json`):

```
swif2-generator-framework/         ← FRAMEWORK_HOME should point here
├── generators/
│   ├── RBHG/
│   │   ├── swif2_RBHG.py
│   │   └── template_RBHG/
│   ├── SPIZG/
│   │   ├── swif2_SPIZG.py
│   │   └── template_SPIZG/
│   └── CobremsDistributions/
├── hddm_scripts/
├── RunPeriods.json
└── ...
```

**Recommended**: Leave `FRAMEWORK_HOME` blank in your config. The script will auto-detect the framework location.

**If you set it explicitly**, use the full path to the framework root:
```
FRAMEWORK_HOME    /work/halld/home/YOUR_USERNAME/swif2-generator-framework
```

**Note**: Both RBHG and SPIZG use the same `FRAMEWORK_HOME` variable, since they share the same framework infrastructure (RunPeriods.json, CobremsDistributions, hddm_scripts, etc.).

### Output Directories

RBHG currently writes generation results below `<FRAMEWORK_HOME>/output/RBHG/`
and creates log directories below `/farm_out/<USERNAME>/RBHG/`, including in
local mode. `GENERATOR_OUTPUT_DIR_BASE` and `LOG_DIR_BASE` are not currently
honored by the RBHG directory creation code. `setup_framework.py` creates
`output/RBHG` during one-time setup.

`MCWRAPPER_OUTPUT_DIR_BASE` configures the later simulation stage, typically
`/volatile/halld/home/{USERNAME}/RBHG/`. Downstream settings can use `NESTED`
(subdirectory of the preceding stage) or `SIBLING` (alongside it) where supported.
See the annotated configs for stage-specific options.

## Interactive/local 10,000-event example

After one-time setup, from the framework root with the GlueX environment loaded
on a JLab Linux host:

```bash
python swif2_gen.py GeneratorConfigExamples/RBHG_10k_local.config --local
# Equivalent alternative:
python swif2_gen.py GeneratorConfigExamples/RBHG_10k_local.config --interactive
```

Run either command, once. Both flags select the same local execution mode.
Omitting the flag submits to SWIF2, regardless of the config filename.
The example runs one serial 10,000-event generator job, then converts the vectors
to HDDM. It retains the physics settings of `RBHG.config` and sets:

```text
STUDY_NAME          interactive_10k
FIXED_EVENT_COUNT   True
NEVENTS_TOTAL       10000
NEVENTS_PERFILE     10000
RUN_PERIOD
TENX_FLAG           False
FLUX_STUDY         False
```

A blank `RUN_PERIOD` selects only the configured `POL_DEG` (90 in this example).
There is no event-count or core-count command-line option. Counts go in the config.
Change `STUDY_NAME` in a copy for a fresh rerun. See the
[main README](../README.md#interactivelocal-run-10000-events-on-one-core) for
output locations and local environment requirements, including `/farm_out`
write access, GNU `sed`, `gfortran`, NumPy, and `hddm_s`.

## RBHG Configurations

- **RBHG.config** - Fully annotated with all parameters explained
- **RBHG_FluxStudy.config** - CPP 2022 counter-by-counter flux-study example
- **RBHG_10k_local.config** - One 10,000-event job; run with `--local` or `--interactive`
- **RBHG_CPP_FFS1.config** - CPP electron-pair form-factor study
- **RBHG_CPP_FFS1_muons.config** - CPP muon-pair form-factor study

## SPIZG Configurations

- **SPIZG.config** - Single π⁰ photoproduction configuration

## Configuration Format

Use uppercase keys and whitespace-separated values, with no equals sign or quotes.
Only the first two whitespace-separated fields are read; `#` starts a comment.
A key on its own sets a blank value.

```text
# Comments start with #
NEVENTS_TOTAL       10000
GLUEX_PHOTON_ENERGY_MIN 7.0  # GeV
RUN_PERIOD                  # blank = one configured polarization
```

## RBHG Run Controls

| Parameter | Description | Example |
|-----------|-------------|---------|
| `RUN_PERIOD` | Run-period preset; blank selects one polarization | `1808` |
| `FIXED_EVENT_COUNT` | Use explicit counts instead of preset statistics | `True` |
| `NEVENTS_TOTAL` | Events per dataset/polarization in fixed-count mode | `10000` |
| `NEVENTS_PERFILE` | Events per generator job/file | `10000` |
| `POL_DEG` | Polarization when `RUN_PERIOD` is blank | `90` |
| `STUDY_NAME` | Study output directory name | `interactive_10k` |

## Generator-Specific Parameters

See individual config files for physics parameters unique to each generator.

## Usage

```bash
# 1. Copy example
cp GeneratorConfigExamples/RBHG.config my_study.config

# 2. Edit parameters
vim my_study.config

# 3. Run
python swif2_gen.py my_study.config
```

See main README.md for full workflow.
