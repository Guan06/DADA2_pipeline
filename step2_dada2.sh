#!/bin/bash

# scripts for demultiplexing raw sequencing data for 16S
# exits whenever a function returns 1
set -e

# load functions
scripts_dir=`dirname $0`
source $scripts_dir/activate
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

for l in $l_list_miseq
do
    # initialize lib. results directory
    rm -f -r $working_dir/"$l"
    mkdir $working_dir/"$l"

    if ($bacteria) then
        mkdir -p $working_dir/$l/Bacteria
        rm -f -r $working_dir_s1/$l/Bac_forward/out/filtN
        rm -f -r $working_dir_s1/$l/Bac_forward/out/cutadapt
        rm -f -r $working_dir_s1/$l/Bac_reverse/out/filtN
        rm -f -r $working_dir_s1/$l/Bac_reverse/out/cutadapt
        $scripts_dir/dada2_bacteria.R $working_dir_s1/$l/ $working_dir/$l/
    fi

    if ($fungi)
    then
        mkdir -p $working_dir/$l/Fungi
        rm -f -r $working_dir_s1/$l/Fun_forward/out/filtN
        rm -f -r $working_dir_s1/$l/Fun_forward/out/cutadapt
        rm -f -r $working_dir_s1/$l/Fun_reverse/out/filtN
        rm -f -r $working_dir_s1/$l/Fun_reverse/out/cutadapt
        $scripts_dir/dada2_fungi.R $working_dir_s1/$l/ $working_dir/$l/ \
            $fungi_FWD $fungi_REV
    fi

    if ($oomycetes)
    then
        mkdir -p $working_dir/$l/Oomycetes
        rm -f -r $working_dir_s1/$l/Oom_forward/out/filtN
        rm -f -r $working_dir_s1/$l/Oom_forward/out/cutadapt
        rm -f -r $working_dir_s1/$l/Oom_reverse/out/filtN
        rm -f -r $working_dir_s1/$l/Oom_reverse/out/cutadapt
        $scripts_dir/dada2_oomycetes.R $working_dir_s1/$l/ $working_dir/$l/ \
            $oomycetes_FWD $oomycetes_REV
    fi

done
