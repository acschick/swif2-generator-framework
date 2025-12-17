# Generator Configuration Examples

Example configuration files for event generators with detailed documentation.

## Quick Start

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

Edit these to match your JLab paths. You can use placeholders that will be automatically replaced:
- `{GENERATOR_TYPE}` - Replaced with RBHG, SPIZG, etc.
- `{USERNAME}` - Replaced with your username
- `{FRAMEWORK_HOME}` - Replaced with framework root path

**Recommended paths** (automatically organized by generator):
```
GENERATOR_OUTPUT_DIR_BASE    {FRAMEWORK_HOME}/output/{GENERATOR_TYPE}/
MCWRAPPER_OUTPUT_DIR_BASE    /volatile/halld/home/{USERNAME}/{GENERATOR_TYPE}/
LOG_DIR_BASE                 /farm_out/{USERNAME}/{GENERATOR_TYPE}/
```

This creates separate directories for RBHG and SPIZG outputs:
- `/volatile/halld/home/user/RBHG/` for RBHG simulations
- `/volatile/halld/home/user/SPIZG/` for SPIZG simulations

You can also use explicit paths without placeholders:
```
GENERATOR_OUTPUT_DIR_BASE    /work/halld/home/YOUR_USERNAME/generators/output/
MCWRAPPER_OUTPUT_DIR_BASE    /volatile/halld/home/YOUR_USERNAME/custom_dir/
LOG_DIR_BASE                 /farm_out/YOUR_USERNAME/my_logs/
```

**Directory organization modes**:
- **EXPLICIT**: Uses base directory + study structure  
- **NESTED**: Each step creates subdirectory in previous step's output
- **SIBLING**: Steps create directories alongside each other at same level

## RBHG Configurations

- **RBHG.config** - Fully annotated with all parameters explained
- **RBHG_debug.config** - Small test run for debugging
- **RBHG_Full2018amo.config** - Full 2018 Spring production

## SPIZG Configurations

- **SPIZG.config** - Single π⁰ photoproduction configuration

## Configuration Format

Simple key-value format:

```bash
# Comments start with #
parameter_name = value
beam_energy = 8.4  # GeV
```

## Common Parameters

Standard across all generators:

| Parameter | Description | Example |
|-----------|-------------|---------|
| `run_period` | GlueX run period | `"2018-01"` |
| `num_events` | Events per job | `50000` |
| `num_jobs` | Parallel jobs | `100` |
| `beam_energy` | Beam energy (GeV) | `8.4` |
| `output_dir` | Output directory | `"MyStudy"` |

## Generator-Specific Parameters

See individual config files for physics parameters unique to each generator.

## Usage

```bash
# 1. Copy example
cp RBHG.config my_study.config

# 2. Edit parameters
vim my_study.config

# 3. Run
python swif2_gen.py my_study.config
```

See main README.md for full workflow.
