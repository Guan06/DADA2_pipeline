#!/bin/bash

set -e

scripts_dir=`dirname $0`

$scripts_dir/step1_demultiplex.sh
$scripts_dir/step2_dada2.sh
$scripts_dir/step3_get_ASV_table.sh
