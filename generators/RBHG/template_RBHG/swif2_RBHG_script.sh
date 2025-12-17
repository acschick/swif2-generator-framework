#!/bin/bash                                                                                                                             
source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh

# Arguments are split by swif2, so we receive them separately:
# $1 = ENVFILE
# $2 = executable path
# $3 = seed
# $4 = job number
# $5 = nevents
# $6 = vectorspath
# $7+ = python hddm conversion command (multiple args)

ENVFILE="$1"
FORTRAN_EXE="$2"
SEED="$3"
JOBNUM="$4"
NEVENTS="$5"
VECTORSPATH="$6"
# Everything from $7 onward is the HDDM conversion command
shift 6
HDDM_CMD=("$@")

# Set up environment
gxenv "$ENVFILE"

# Run Fortran executable with command-line arguments
echo "About to run Fortran: $FORTRAN_EXE $SEED $JOBNUM $NEVENTS"
"$FORTRAN_EXE" "$SEED" "$JOBNUM" "$NEVENTS"
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

