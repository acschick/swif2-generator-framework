# Environment Configuration Presets

This directory contains preset environment setup scripts for different run periods.

## Available Configurations

### 2018 Spring (1801)
- `2018_01.csh` / `2018_01.sh`
- Reconstruction: recon-2018_01-ver02_29.xml
- Analysis: analysis-2018_01-ver03.xml

### 2018 Fall (1808)
- `2018_08.csh` / `2018_08.sh`
- Reconstruction: recon-2018_08-ver02_28.xml
- Analysis: version_4.20.0.xml

### 2017 Data
- `2017.csh` / `2017.sh`
- Reconstruction: TBD based on newest 2017 reconstruction

### CPP (2022-05, run period 2205)
- `2205_cpp.csh` / `2205_cpp.sh`
- CPP-specific geometry and calibration context

## Usage

The correct environment file is automatically selected based on:
1. Run period from your RBHG generation
2. User's shell (detected from $SHELL environment variable)

You can also manually specify in `RBHG.config`:
```
MCWRAPPER_RECON_ENV          /path/to/specific/version.xml
```

## Notes

- `.csh` files are for tcsh/csh users
- `.sh` files are for bash/sh users
- Username (HOME directory) is automatically configured
