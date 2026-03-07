#!/usr/bin/env python3
"""
Update all rbhg_config.json files in CPP_FFS1 with a 10x reduction in events.
Run this BEFORE running the regeneration bash script.
"""

import json
import os
import sys

BASE = "/work/halld/home/acschick/channels/swif2-generator-framework/output/RBHG/CPP_FFS1"

# All 12 simulation directories with their dataset type
DATASETS = {
    "BC_FF1":    ["2205_135DEG_FF1_IRADOFF", "2205_45DEG_FF1_IRADOFF", "2205_AMO_FF1_IRADOFF"],
    "BC_FFN":    ["2205_135DEG_FFN_IRADOFF", "2205_45DEG_FFN_IRADOFF", "2205_AMO_FFN_IRADOFF"],
    "SIM_FF1":   ["2205_135DEG_FF1_IRADOFF", "2205_45DEG_FF1_IRADOFF", "2205_AMO_FF1_IRADOFF"],
    "qDATAq_FFN":["2205_135DEG_FFN_IRADOFF", "2205_45DEG_FFN_IRADOFF", "2205_AMO_FFN_IRADOFF"],
}

dry_run = "--dry-run" in sys.argv

if dry_run:
    print("DRY RUN MODE - no changes will be written")
print()

total_updated = 0
errors = []

for dataset, dirs in DATASETS.items():
    for sim_dir in dirs:
        config_path = os.path.join(BASE, dataset, sim_dir, "rbhg_config.json")
        
        if not os.path.exists(config_path):
            print(f"  WARNING: Not found: {config_path}")
            errors.append(config_path)
            continue
        
        with open(config_path, "r") as f:
            config = json.load(f)
        
        ec = config["rbhg_config"]["event_counts"]
        old_total  = ec["total_events"]
        old_jobs   = ec["total_jobs"]
        old_epf    = ec["events_per_file"]
        
        new_total  = old_total  // 10
        new_jobs   = old_jobs   // 10
        
        print(f"{dataset}/{sim_dir}")
        print(f"  total_events: {old_total:>12,} → {new_total:>11,}")
        print(f"  total_jobs:   {old_jobs:>12,} → {new_jobs:>11,}")
        print(f"  events_per_file: {old_epf} (unchanged)")
        
        if not dry_run:
            ec["total_events"] = new_total
            ec["total_jobs"]   = new_jobs
            with open(config_path, "w") as f:
                json.dump(config, f, indent=2)
            print(f"  ✓ Written")
        else:
            print(f"  (dry run - not written)")
        print()
        total_updated += 1

print(f"{'='*60}")
if dry_run:
    print(f"DRY RUN: Would update {total_updated} config files")
else:
    print(f"Updated {total_updated} config files")

if errors:
    print(f"\nWARNING: {len(errors)} files not found:")
    for e in errors:
        print(f"  {e}")
