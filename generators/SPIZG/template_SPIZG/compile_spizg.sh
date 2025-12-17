#!/bin/tcsh
# Compilation script for SPIZG generator (C++ with ROOT and GSL)
# Usage: compile_spizg.sh <source_file.cpp> <output_executable>

# Check arguments
if ( $#argv != 2 ) then
    echo "Usage: compile_spizg.sh <source_file.cpp> <output_executable>"
    exit 1
endif

set SOURCE_FILE = $1
set OUTPUT_EXEC = $2

# Check if source file exists
if ( ! -f $SOURCE_FILE ) then
    echo "Error: Source file $SOURCE_FILE not found"
    exit 1
endif

# ROOT configuration
# root-config should be available if ROOT environment is loaded
set ROOT_CFLAGS = `root-config --cflags`
set ROOT_LIBS = `root-config --libs`

# GSL configuration  
# Assuming GSL is installed in standard location or via module
set GSL_CFLAGS = "-I/usr/include/gsl"
set GSL_LIBS = "-lgsl -lgslcblas -lm"

# Compile command
echo "Compiling $SOURCE_FILE -> $OUTPUT_EXEC"
echo "ROOT flags: $ROOT_CFLAGS"
echo "ROOT libs: $ROOT_LIBS"
echo "GSL flags: $GSL_CFLAGS"
echo "GSL libs: $GSL_LIBS"

g++ -o $OUTPUT_EXEC $SOURCE_FILE \
    $ROOT_CFLAGS $ROOT_LIBS \
    $GSL_CFLAGS $GSL_LIBS \
    -O2 -std=c++17

if ( $status == 0 ) then
    echo "Compilation successful: $OUTPUT_EXEC"
    exit 0
else
    echo "Compilation failed"
    exit 1
endif
