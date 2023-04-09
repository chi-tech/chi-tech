#!/bin/bash
base_name=$(basename $4)

$1 -np $2 $3 $4 --supress_beg_end_timelog --suppress_color > ZUTOut.data

printf "diff command\n"
diff ZUTOut.data tests/unit_test_inputs/GoldenOutputs/${base_name}.gold 1>&2

exit $?