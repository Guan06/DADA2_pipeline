#!/bin/bash

# scripts for demultiplexing raw sequencing data for 16S
# exits whenever a function returns 1
set -e

# load functions
scripts_dir=`dirname $0`
source $scripts_dir/activate.sh
source $scripts_dir/config.sh

log() {
    echo $(date -u)": "$1 >> $logfile
}

mkdir -p $working_dir

working_dir_s1=$working_dir/01.split_fq
working_dir=$working_dir/02.dada2

output=$working_dir/"output.txt"
logfile=$working_dir/"log.txt"

mkdir -p $working_dir
#rm -f -r $working_dir/*

for l in $l_list_roche
do
    # initialize lib. results directory
    rm -f -r $working_dir/"$l"
    mkdir $working_dir/"$l"

    $scripts_dir/dada2_bacteria_454.R $working_dir_s1/$l/ $working_dir/$l/

done
