# RBHG Bremsstrahlung File Naming Convention

This document describes the naming convention and organization for bremsstrahlung 
distribution files used by the RBHG event generator.

## Current File Structure (template_RBHG/)

### Legacy Fixed-Width Files
- `CobremsDistribution.dat` - GlueX legacy (for backward compatibility)  
- `CPP_Cobrems.dat` - CPP legacy

### New Variable-Width Files
- `CobremsDistribution_varbin.dat` - GlueX variable-width default
- `CPP_Cobrems_varbin.dat` - CPP variable-width default

## Future Polarization-Specific Files

### Planned GlueX Files
```
CobremsDistribution_varbin_0DEG.dat    # 0° polarization
CobremsDistribution_varbin_45DEG.dat   # 45° polarization  
CobremsDistribution_varbin_90DEG.dat   # 90° polarization
CobremsDistribution_varbin_135DEG.dat  # 135° polarization
CobremsDistribution_varbin_AMO.dat     # Amorphous (unpolarized)
```

### Planned CPP Files
```
CPP_Cobrems_varbin_45DEG.dat    # 45° polarization
CPP_Cobrems_varbin_135DEG.dat   # 135° polarization  
CPP_Cobrems_varbin_AMO.dat      # Amorphous (unpolarized)
```

## Configuration Options

### Current Usage
```
# RBHG.config
COBREMS_VARBIN = False    # Uses legacy fixed-width files
COBREMS_VARBIN = True     # Uses new variable-width files
```

### Future Usage (Planned)
```
# RBHG.config - Manual file override
COBREMS_FILE_GLUEX = CobremsDistribution_run42559_varbin.dat  # Custom file
COBREMS_FILE_CPP = CPP_Cobrems_run101582_varbin.dat          # Custom file

# Automatic polarization selection (future)
COBREMS_POLARIZATION_SPECIFIC = True   # Use pol-specific files when available
```

## File Format

### Legacy Format (2 columns)
```
# Energy(GeV)  Rate
7.00  1234.5
7.04  1456.7
...
```

### Variable-Width Format (3 columns + metadata)
```
# RBHG Variable-Width Bremsstrahlung Distribution  
# SOURCE: BGRate_42559.root : dRtdkH1
# ENERGY_LIMITS (GeV):
# E_MIN_AVAILABLE 2.907218
# E_MAX_AVAILABLE 11.406374
# Format: lowEdge(GeV)  highEdge(GeV)  counts
2.907218  2.907314  1234.56
2.907314  2.908102  2345.67
...
```

## Implementation Status

### ✅ Completed
- Variable-width bin support with binary search
- Automatic energy range validation  
- Backward compatibility with fixed-width files
- Smart default file selection framework

### 🚧 Future Features  
- Polarization-specific file auto-selection
- Batch extraction script for multiple polarizations
- Configuration file override options
- Run-specific file naming conventions

## Usage Examples

### Extract Single Distribution
```bash
python3 extract_varbin_cobrems.py BGRate_42559.root dRtdkH1 CobremsDistribution_varbin.dat
```

### Extract Multiple Polarizations (Future)
```bash
python3 extract_varbin_cobrems.py --batch BGRate_all_pol.root dRtdkH1_
# Creates: CobremsDistribution_varbin_0DEG.dat, etc.
```

### Use Custom File  
```bash
# RBHG.config
COBREMS_FILE_GLUEX = CobremsDistribution_run42559_0DEG_varbin.dat
```