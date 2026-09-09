#!/usr/bin/env python3
"""Collect integrated RBHG flux-study cross sections from SWIF2 stdout files."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Match, Optional, Pattern, Set, Tuple


FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"

COMMAND_RE = re.compile(
    rf"Command-line args:\s*seed=\s*(?P<seed>\d+)\s+"
    rf"job=\s*(?P<job>\d+)\s+nevents=\s*(?P<events>\d+)"
)
ENERGY_RE = re.compile(
    rf"Runtime energies:\s*E_lo=\s*(?P<low>{FLOAT})\s+"
    rf"E_hi=\s*(?P<high>{FLOAT})\s+E_coherent=\s*(?P<center>{FLOAT})"
)
X_RANGE_RE = re.compile(
    rf"Configured x range:\s*(?P<minimum>{FLOAT})\s+to\s+(?P<maximum>{FLOAT})"
)
INTEGRAL_RE = re.compile(
    rf"test events\s+(?P<tests>\d+)\s+Egamma\s+(?P<energy>{FLOAT})\s+"
    rf"total cross section nb\s+(?P<cross_section>{FLOAT})\s+"
    rf"percent change\s+(?P<change>{FLOAT})"
)
CROSS_MAX_RE = re.compile(rf"CROSS_MAX FOUND:\s*(?P<cross_max>{FLOAT})")
COUNTER_RE = re.compile(r"(?P<detector>TAG[HM])_(?P<counter>\d{1,3})(?:\D|$)", re.I)
EXIT_RE = re.compile(r"Fortran executable completed with exit code:\s*(?P<code>-?\d+)")

FIELDNAMES = [
    "study",
    "workflow",
    "run_period",
    "polarization",
    "detector",
    "counter",
    "job_index",
    "seed",
    "requested_events",
    "E_lo_GeV",
    "E_hi_GeV",
    "E_coherent_GeV",
    "x_min",
    "x_max",
    "cross_section_nb",
    "integration_test_events",
    "final_percent_change",
    "integration_estimate_count",
    "cross_max_sampling_weight",
    "event_generation_started",
    "fortran_exit_code",
    "status",
    "stdout_file",
]

POLARIZATION_RE = re.compile(r"(?:^|_)(?P<polarization>(?:0|45|90|135)DEG|AMO)(?:_|$)", re.I)
RUN_PERIOD_RE = re.compile(r"(?:^|_)(?P<run_period>\d{4})(?:_|$)")


def as_float(value: str) -> float:
    """Parse either a C/Python E exponent or a Fortran D exponent."""
    return float(value.replace("D", "E").replace("d", "e"))


def last_match(pattern: Pattern[str], text: str) -> Optional[Match[str]]:
    matches = list(pattern.finditer(text))
    return matches[-1] if matches else None


def counter_from_path(path: Path) -> Tuple[str, Optional[int]]:
    match = last_match(COUNTER_RE, path.name)
    if match is None:
        match = last_match(COUNTER_RE, str(path.parent))
    if match is None:
        return "", None
    return match.group("detector").upper(), int(match.group("counter"))


def workflow_metadata(path: Path) -> Tuple[str, str, str, str]:
    """Extract study, workflow, run period, and polarization from farm paths."""
    workflow_path = next(
        (parent for parent in path.parents if parent.name.upper().startswith("FLUX_TAG")),
        None,
    )
    workflow = workflow_path.name if workflow_path else path.parent.name
    polarization_match = POLARIZATION_RE.search(workflow)
    run_period_match = RUN_PERIOD_RE.search(workflow)
    study = ""
    if workflow_path is not None:
        generation_path = workflow_path.parent
        if generation_path.name == "generation":
            study = generation_path.parent.name
    return (
        study,
        workflow,
        run_period_match.group("run_period") if run_period_match else "",
        polarization_match.group("polarization").upper() if polarization_match else "",
    )


def parse_stdout(path: Path) -> Dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="replace")
    command = last_match(COMMAND_RE, text)
    energy = last_match(ENERGY_RE, text)
    x_range = last_match(X_RANGE_RE, text)
    integrals = list(INTEGRAL_RE.finditer(text))
    final_integral = integrals[-1] if integrals else None
    cross_max = last_match(CROSS_MAX_RE, text)
    exit_match = last_match(EXIT_RE, text)
    detector, counter = counter_from_path(path)
    study, workflow, run_period, polarization = workflow_metadata(path)

    # Flux job indices are zero-based and counters are one-based. This fallback
    # supports logs that were renamed and no longer contain TAGH_###/TAGM_###.
    if counter is None and command is not None:
        counter = int(command.group("job")) + 1

    row: Dict[str, object] = {field: "" for field in FIELDNAMES}
    row.update(
        {
            "study": study,
            "workflow": workflow,
            "run_period": run_period,
            "polarization": polarization,
            "detector": detector,
            "counter": counter if counter is not None else "",
            "job_index": int(command.group("job")) if command else "",
            "seed": int(command.group("seed")) if command else "",
            "requested_events": int(command.group("events")) if command else "",
            "E_lo_GeV": as_float(energy.group("low")) if energy else "",
            "E_hi_GeV": as_float(energy.group("high")) if energy else "",
            "E_coherent_GeV": as_float(energy.group("center")) if energy else "",
            "x_min": as_float(x_range.group("minimum")) if x_range else "",
            "x_max": as_float(x_range.group("maximum")) if x_range else "",
            "cross_section_nb": (
                as_float(final_integral.group("cross_section")) if final_integral else ""
            ),
            "integration_test_events": (
                int(final_integral.group("tests")) if final_integral else ""
            ),
            "final_percent_change": (
                as_float(final_integral.group("change")) if final_integral else ""
            ),
            "integration_estimate_count": len(integrals),
            "cross_max_sampling_weight": (
                as_float(cross_max.group("cross_max")) if cross_max else ""
            ),
            "event_generation_started": "Starting event generation" in text,
            "fortran_exit_code": int(exit_match.group("code")) if exit_match else "",
            "stdout_file": str(path.resolve()),
        }
    )

    if final_integral is None:
        row["status"] = "missing_cross_section"
    elif int(final_integral.group("tests")) < 100_000_000:
        row["status"] = "partial_integration"
    elif len(integrals) < 4:
        row["status"] = "incomplete_iterations"
    else:
        row["status"] = "cross_section_complete"
    return row


def find_stdout_files(inputs: List[Path]) -> List[Path]:
    files: Set[Path] = set()
    for item in inputs:
        if item.is_file():
            files.add(item)
        elif item.is_dir():
            files.update(path for path in item.rglob("stdout_*.out") if path.is_file())
        else:
            print(f"warning: input does not exist: {item}", file=sys.stderr)
    return sorted(files)


def sort_key(row: Dict[str, object]) -> Tuple[object, ...]:
    detector_order = {"TAGH": 0, "TAGM": 1}
    return (
        str(row["run_period"]),
        str(row["polarization"]),
        detector_order.get(str(row["detector"]), 2),
        int(row["counter"]) if row["counter"] != "" else 10_000,
        str(row["stdout_file"]),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Recursively parse RBHG stdout files and tabulate the final "
            "integrated cross section from each flux-counter job."
        )
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="stdout file(s), workflow log directory, or parent study directory",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("rbhg_flux_cross_sections.csv"),
        help="output CSV path (default: ./rbhg_flux_cross_sections.csv)",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    stdout_files = find_stdout_files(args.inputs)
    if not stdout_files:
        print("error: no stdout_*.out files found", file=sys.stderr)
        return 1

    rows = sorted((parse_stdout(path) for path in stdout_files), key=sort_key)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    complete = sum(row["status"] == "cross_section_complete" for row in rows)
    incomplete = len(rows) - complete
    print(f"Wrote {len(rows)} rows to {args.output}")
    print(f"Complete cross sections: {complete}; incomplete/missing: {incomplete}")
    if incomplete:
        print("Review rows whose status is not 'cross_section_complete'.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
