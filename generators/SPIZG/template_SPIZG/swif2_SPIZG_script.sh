#!/bin/bash
source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh

# Arguments are split by swif2, so we receive them separately:
# $1 = ENVFILE
# $2 = executable path (C++ compiled binary)
# $3 = job number
# $4 = nevents
# $5 = vectorspath
# $6+ = python hddm conversion command (multiple args)

ENVFILE="$1"
CPP_EXE="$2"
JOBNUM="$3"
NEVENTS="$4"
VECTORSPATH="$5"
# Everything from $6 onward is the HDDM conversion command
shift 5
HDDM_CMD=("$@")

# Set up environment
gxenv "$ENVFILE"

# Run C++ executable with command-line arguments: jobnum nevents
echo "About to run C++ generator: $CPP_EXE $JOBNUM $NEVENTS"
"$CPP_EXE" "$JOBNUM" "$NEVENTS"
exitcode=$?
echo "C++ executable completed with exit code: $exitcode"

# Check if C++ executable succeeded
if [ $exitcode -ne 0 ]; then
    echo "ERROR: C++ executable failed with exit code $exitcode"
    exit $exitcode
fi

# Change to vectors directory for HDDM conversion
echo "Changing to vectors directory: $VECTORSPATH"
cd "$VECTORSPATH" || exit 1

# Run HDDM conversion
echo "Running HDDM conversion: ${HDDM_CMD[@]}"
"${HDDM_CMD[@]}"
exitcode=$?

if [ $exitcode -ne 0 ]; then
    echo "ERROR: HDDM conversion failed with exit code $exitcode"
    exit $exitcode
fi

echo "SPIZG job completed successfully"
exit 0
