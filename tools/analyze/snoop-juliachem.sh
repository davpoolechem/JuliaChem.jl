#/bin/bash

#perform snoopcompile work
echo "----------------------------------------"
echo "    PERFORMING SNOOPCOMPILE ANALYSIS    "
echo "----------------------------------------"
if [ -d src/core/snoop ]; then
   rm -rf snoop
   rm snoop.csv
fi
JULIA_NUM_THREADS=1 julia tools/snoop.jl example_scripts/minimal-rhf.jl example_inputs/S22-01.json
