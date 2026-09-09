# RBHG utility scripts

## Collecting flux-study cross sections

`collect_flux_cross_sections.py` recursively scans SWIF2 stdout files and
writes one CSV row per flux-counter job. It records the last integrated
`total cross section nb` estimate. `CROSS_MAX` is retained only as a diagnostic
sampling weight and is not used as the reported physical cross section.

Run it on the common parent of the TAGH and TAGM log directories. The scan is
recursive, so this one command reads both polarizations and both detectors:

```bash
python3 generators/RBHG/scripts/collect_flux_cross_sections.py \
  /farm_out/$USER/RBHG/CPP_FluxStudy3/generation \
  -o CPP_FluxStudy3_cross_sections.csv
```

You can also pass several directories or individual stdout files. Rows with a
complete four-stage integration end with `status=cross_section_complete`.
Missing or unfinished calculations remain in the CSV with a diagnostic status
so failed and timed-out counters are visible. The workflow, run period, and
polarization columns keep identically numbered counters from different
workflows distinct.
