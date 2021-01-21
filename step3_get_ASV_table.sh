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
working_dir_s2=$working_dir/02.dada2
working_dir=$working_dir/03.ASV_tab

output=$working_dir/"output.txt"
logfile=$working_dir/"log.txt"

mkdir -p $working_dir
#rm -f -r $working_dir/*

if ($bacteria)
then
    mkdir -p $working_dir/Bacteria
    $scripts_dir/get_ASV_tab_bacteria.R \
        $working_dir_s2/ $working_dir/Bacteria/
    `less $working_dir/Bacteria/ASV_map.txt |
        sed 's/ASV/>ASV/; s/\t/\n/' > $working_dir/Bacteria/ASV.fasta`
fi

if ($fungi)
then
    mkdir -p $working_dir/Fungi
    $scripts_dir/get_ASV_tab_fungi.R \
        $working_dir_s2/ $working_dir/Fungi/
    `less $working_dir/Fungi/ASV_map.txt |
        sed 's/ASV/>ASV/; s/\t/\n/' > $working_dir/Fungi/ASV.fasta`
fi

if ($oomycetes)
then
    mkdir -p $working_dir/Oomycetes
    $scripts_dir/get_ASV_tab_oomycetes.R \
        $working_dir_s2/ $working_dir/Oomycetes/
    `less $working_dir/Oomycetes/ASV_map.txt |
    sed 's/ASV/>ASV/; s/\t/\n/' >$working_dir/Oomycetes/ASV.fasta`

    java -Xmx1g -jar \
        /biodata/dep_psl/grp_psl/thiergart/RDPTools/classifier.jar \
        -c 0.5 -t \
       $scripts_dir/tax_oomycetes/rRNAClassifier.properties \
        -o  $working_dir/Oomycetes/ASV_taxonomy.txt \
        $working_dir/Oomycetes/ASV.fasta
fi
