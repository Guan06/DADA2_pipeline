Pipeline for data processing for amplicon sequencing data, suitable for
bacterial, fungal and oomycetal community profiling.
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

4) Run scripts step by step
