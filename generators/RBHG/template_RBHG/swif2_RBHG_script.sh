#!/bin/bash                                                                                                                             
source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh

# Arguments are split by swif2, so we receive them separately:
# $1 = ENVFILE
# $2 = executable path
# $3 = seed
# $4 = job number
# $5 = nevents
# Flux jobs add: $6=--energy-override, $7=E_lo, $8=E_hi,
# $9=E_coherent, $10=vectorspath, and $11+=HDDM command.
# Normal jobs retain $6=vectorspath and $7+=HDDM command.

ENVFILE="$1"
FORTRAN_EXE="$2"
SEED="$3"
JOBNUM="$4"
NEVENTS="$5"
ENERGY_ARGS=()
if [ "$6" = "--energy-override" ]; then
    ENERGY_ARGS=("$7" "$8" "$9")
    VECTORSPATH="${10}"
    shift 10
else
    VECTORSPATH="$6"
    shift 6
fi
HDDM_CMD=("$@")

# Set up environment
gxenv "$ENVFILE"

# Run Fortran executable with command-line arguments
echo "About to run Fortran: $FORTRAN_EXE $SEED $JOBNUM $NEVENTS ${ENERGY_ARGS[*]}"
"$FORTRAN_EXE" "$SEED" "$JOBNUM" "$NEVENTS" "${ENERGY_ARGS[@]}"
exitcode=$?
echo "Fortran executable completed with exit code: $exitcode"

# Check if Fortran succeeded
if [ $exitcode -ne 0 ]; then
    echo "ERROR: Fortran executable failed with exit code $exitcode"
    exit $exitcode
fi

# Change to vectors directory for HDDM conversion
echo "Changing to vectors directory: $VECTORSPATH"
cd "$VECTORSPATH" || exit 1

# Run HDDM conversion
echo "Running HDDM conversion: ${HDDM_CMD[@]}"
"${HDDM_CMD[@]}"
