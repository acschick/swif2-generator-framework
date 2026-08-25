# RBHG Energy-Binning Schemes

This directory stores reusable photon-energy schemes for RBHG flux studies and
future counter-by-counter generation tests.

Each scheme is a whitespace-delimited text file. Blank lines are ignored,
comments begin with `#`, and data rows have exactly three fields:

```text
DETECTOR  COUNTER  E_CENTER_GEV
TAGH      1        11.338816
TAGM      1         5.911535
```

Flux schemes must list TAGH counters 1-274 and TAGM counters 1-102 exactly
once. No counters are excluded by the loader. Add provenance as `# key: value`
metadata so it is copied into each generated `rbhg_config.json`.

The default CPP 2022 scheme uses counter centers calculated at the mean raw
endpoint energy of 11.549508 GeV. Flux jobs generate uniformly within +/-0.05
MeV of each center and evaluate the integrated cross section at the center.

