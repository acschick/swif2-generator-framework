# JANA Configuration Presets

This directory contains preset JANA configuration files for common MCWrapper analysis scenarios.

## Available Configurations

### `jana_CPPlaunch1_epem.config`
- **Use for**: CPP Launch 1 (2022-05, run period 2205) e+e- analysis
- **Plugins**: ReactionFilter, monitoring_hists
- **Features**: Primary e+e- reaction plus additional EM reactions (pi+pi-, pi0pi+pi-, Dalitz) for ML training
- **Target**: Pb208
- **Beam Bunches**: ±2

### `jana_CPPlaunch1_mupmum.config`
- **Use for**: CPP Launch 1 (2022-05, run period 2205) mu+mu- analysis
- **Plugins**: ReactionFilter, monitoring_hists
- **Features**: Primary mu+mu- reaction plus additional reactions (pi+pi-, pi0pi+pi-) for ML training
- **Target**: Pb208
- **Beam Bunches**: ±2

### `jana_2205_mupmum_ReaFil.config`
- **Use for**: 2022-05 mu+mu- basic analysis
- **Plugins**: ReactionFilter, mcthrown_tree, danarest, monitoring_hists
- **Features**: Simple mu+mu- ReactionFilter
- **Target**: Pb208

### `2018_epem_standard.config`
- **Use for**: 2018 Spring (1801) and Fall (1808) e+e- analysis
- **Plugins**: ReactionFilter, danarest, monitoring_hists
- **Features**: Standard E/p cuts, BCAL/FCAL verbose output
- **Background**: Matches to recon-2018_01-ver02 or recon-2018_08-ver02

### `2017_epem_recon_only.config`
- **Use for**: 2017 data, reconstruction only (no ReactionFilter)
- **Plugins**: danarest, monitoring_hists
- **Features**: FCAL bad block removal, basic reconstruction
- **Note**: Does not run ReactionFilter due to 2017 trigger limitations

### `2205_cpp_epem.config`
- **Use for**: CPP (2022-05, run period 2205) e+e- analysis
- **Plugins**: ReactionFilter, danarest, CPP trigger emulation
- **Features**: CPP trigger, Kalman vertex at 1cm, pi/mu neural net support
- **Target**: Pb208

## Usage

In `RBHG.config`, set:
```
MCWRAPPER_JANA_CONFIG        2018_epem_standard.config
```

Or use an absolute path to a custom configuration:
```
MCWRAPPER_JANA_CONFIG        /path/to/my/custom_jana.config
```

## Creating Custom Configs

Copy one of the existing configs and modify as needed. Key sections:

1. **PLUGINS** - Which plugins to load
2. **Plugin-specific options** - Parameters for each plugin
3. **E/p cuts** - If using ReactionFilter
4. **Trigger settings** - Bypass or emulation
5. **Verbose output** - BCAL/FCAL extra data

See individual config files for detailed examples.
