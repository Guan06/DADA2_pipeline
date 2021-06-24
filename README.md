Pipeline for data processing for amplicon sequencing data, suitable for
bacterial, fungal and oomycetal community profiling.
Usearch is needed for preprocessing the data, specifically, demultiplexing.  

To run the workflow:

1) Prepare your input files and put it in the data_dir (later defined in config.sh)

Name sequence fastq files as:  

    prefix_forward_reads.fastq.gz  
    prefix_reverse_reads.fastq.gz  
    prefix_barcodes.fastq.gz  

Prepare mapping file and name them as:

    prefix_mapping.txt (for Bacteria)
    prefix_mapping_ITSf.txt (for Fungi if applicable)
    prefix_mapping_ITSo.txt (for Oomycetes if applicable)

2) Validate the prepared mapping file

    ./activate  
    validate_mapping_file.py -m $data_dir/prefix_mapping.txt -o ./

3) Edit config file:

    Define the working_dir and data_dir in config file (./config.sh)  
    Define the profiled kingdom and amplification primer set if applicable(for fungi and oomycetes)  


4) Run scripts step by step or all together (take MiSeq data as example here)
    ./step1_demultiplex.sh

    Checking the output in the folder (take Bacteria as example here):

    less $working_dir/01.split_fq/$l_list_miseq/Bac_forward/split_library_log.txt

    In this file, the number of reads were shown and for those who has less than
    10 (or 20 to be more strict) samples, remove the demultiplexed file, where
    you could find in subfolder under the same folder(e.g. sample OD1):

    rm ./out/OD1.fastq                  ## forward read file of this sample
    rm ../Bac_reverse/out/OD1.fastq     ## reverse read file of the same sample

    Then run the later steps as follow:

    ./step2_dada2.sh
    ./step3_get_ASV_table.sh

    or

    ./all_in_one.sh
