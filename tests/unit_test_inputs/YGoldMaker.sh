#!/bin/bash
base_name=$(basename $4)

$1 -np $2 $3 $4 --supress_beg_end_timelog > tests/unit_test_inputs/GoldenOutputs/${base_name}.gold

exit $?