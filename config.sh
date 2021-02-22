#!/bin/bash

# config file for amplicon data analysis
# guan@mpipz.mpg.de, modified from garridoo@mpipz.mpg.de

## set the directory for the results to be saved
working_dir="/netscratch/somewhere"

## set the directory where the sequencing data is stored
data_dir="/biodata/somewhere"

## define the profiled kingdom
bacteria=false

fungi=true

## for ITS2
fungi_FWD="GTGARTCATCGAATCTTTG"
fungi_REV="TCCTCCGCTTATTGATATGC"

## for ITS1
#fungi_FWD="CTTGGTCATTTAGAGGAAGTAA"
#fungi_REV="GCTGCGTTCTTCATCGATGC"

oomycetes=false

## for ITS1o
oomycetes_FWD="CGGAAGGATCATTACCAC"
oomycetes_REV="AGCCTAGACATCCACTGCTG"

## set the list of sequencing data
## should be the same as "prefix" in README.md
l_list_miseq="prefix"
