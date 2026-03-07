#!/bin/bash
# =============================================================================
# hddm_batch_worker.sh
#
# Worker script for swif2 jobs that batch-convert txt → HDDM files using
# ascii2hddm.py.  This script is generated/called by submit_hddm_regen_swif2.py.
#
# Arguments:
#   $1  = ENVFILE   - path to GlueX version XML (e.g. version.xml)
#   $2  = ASCII2HDDM_PATH - absolute path to ascii2hddm.py
#   $3  = PARTICLE  - final state tag, e.g. "epem"
#   $4  = TARGET    - target tag, e.g. "pb208"
#   $5  = RUN_NUMBER
#   $6+ = alternating pairs: txt_file hddm_file txt_file hddm_file ...
#
# Usage (called by swif2):
#   bash hddm_batch_worker.sh <envfile> <ascii2hddm.py> epem pb208 101582 \
#       vectors_0.txt vectors_0.hddm vectors_1.txt vectors_1.hddm ...
# =============================================================================

set -u

# --- Environment setup ---
source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
ENVFILE="$1"
gxenv "$ENVFILE"

# --- Parse fixed arguments ---
ASCII2HDDM="$2"
PARTICLE="$3"
TARGET="$4"
RUN_NUMBER="$5"
shift 5   # remaining args: optional --vertex, then alternating txt/hddm pairs

# --- Optional --vertex flag ---
VERTEX_ARGS=()
if [[ "${1:-}" == "--vertex" ]]; then
    shift   # consume --vertex
    # consume exactly 4 values: vx vy zmin zmax
    VERTEX_ARGS=("--vertex" "$1" "$2" "$3" "$4")
    shift 4
fi

TOTAL=0
FAILED=0
CONVERTED=0

echo "======================================================="
echo "  hddm_batch_worker.sh"
echo "  ENVFILE:    $ENVFILE"
echo "  ascii2hddm: $ASCII2HDDM"
echo "  particle:   $PARTICLE  target: $TARGET  run: $RUN_NUMBER"
if [[ ${#VERTEX_ARGS[@]} -gt 0 ]]; then
echo "  vertex:     ${VERTEX_ARGS[*]}"
fi
echo "======================================================="

# --- Convert each pair ---
while [ "$#" -ge 2 ]; do
    TXT_FILE="$1"
    HDDM_FILE="$2"
    shift 2
    TOTAL=$((TOTAL + 1))

    if [ ! -f "$TXT_FILE" ]; then
        echo "  WARNING: input file not found, skipping: $TXT_FILE"
        FAILED=$((FAILED + 1))
        continue
    fi

    echo "  [$TOTAL] Converting: $(basename "$TXT_FILE") → $(basename "$HDDM_FILE")"
    python3 "$ASCII2HDDM" "$PARTICLE" "$TARGET" "$TXT_FILE" "$HDDM_FILE" --run "$RUN_NUMBER" "${VERTEX_ARGS[@]}"
    EXIT_CODE=$?

    if [ $EXIT_CODE -ne 0 ]; then
        echo "  ERROR: ascii2hddm.py failed (exit $EXIT_CODE) for $TXT_FILE"
        FAILED=$((FAILED + 1))
    else
        CONVERTED=$((CONVERTED + 1))
    fi
done

echo "========================================================"
echo "  Done. Converted: $CONVERTED / $TOTAL  (failed: $FAILED)"
echo "========================================================"

if [ "$FAILED" -gt 0 ]; then
    exit 1
fi
exit 0
