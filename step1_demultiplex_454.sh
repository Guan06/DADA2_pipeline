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

working_dir=$working_dir/01.split_fq

output=$working_dir/"output.txt"
logfile=$working_dir/"log.txt"

mkdir -p $working_dir
rm -f -r $working_dir/*

# demultiplexing reads1 and reads2 files
for l in $l_list_roche
do
    # initialize lib. results directory
    log "["$l"] initializing the working_dir for roche step 1..."
    rm -f -r $working_dir/"$l"
    mkdir $working_dir/"$l"

       log "["$l"] demultiplexing Bacteria reads..."
    bc_len=`less $data_dir/$l\_mapping.txt |tail -n1 |
            awk '{print $2}' |wc -c`
    let "bc_len=$bc_len - 1"

    split_libraries.py -f $data_dir/"$l".fasta \
        -q $data_dir/"$l".qual \
        -m $data_dir/"$l"_mapping.txt \
        -b $bc_len \
        -s 0 \
        -d \
        -o $working_dir/$l \
        &>> $output

    convert_fastaqual_fastq.py -f $working_dir/$l/seqs.fna \
        -q $working_dir/$l/seqs_filtered.qual \
        -o $working_dir/$l \
        &>> $output

    split_sequence_file_on_sample_ids.py --file_type fastq \
        -i $working_dir/$l/seqs.fastq \
        -o $working_dir/$l/out \
        &>> $output

    log "["$l"] step 0 finished."
done
