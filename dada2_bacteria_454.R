#!//biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(dada2)
library(ggplot2)
library(ShortRead)
library(Biostrings)

args = commandArgs(trailingOnly = TRUE)
path <- args[1]
output <- args[2]

########################################### trim Roche 454 data
## change primer sequence for different amplicon region
FWD <- "AGAGTTTGATCMTGGCTCAG"
REV <- "GWATTACCGCGGCKGCTG"

allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    # Convert to Biostrings object
    dna <- DNAString(primer)

    orients <- c(Forward = dna, Complement = complement(dna),
                 Reverse = reverse(dna), RevComp = reverseComplement(dna))

    # Convert back to character vector
    return(sapply(orients, toString))
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


path1 <- paste0(path, "/out")
output_bac <- output

fnFs <- sort(list.files(path1, pattern=".fastq", full.names = TRUE))

# count the number of times the primers appear
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

rbind(FWD.Reads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]),
      REV.Reads = sapply(REV.orients, primerHits, fn = fnFs[[1]]))

# change this path to your local path where cutadapt installed
cutadapt <- "/opt/share/software/bin/cutadapt"
path.cut <- file.path(path1, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)

fnFs.cut <- file.path(path.cut, basename(fnFs))

FWD.RC <- dada2:::rc(FWD)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
flag <- paste("-g", REV, "-a", FWD.RC)

# Run Cutadapt
for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(flag, " -o", fnFs.cut[i], fnFs[i]))
}

rbind(FWD.READs = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.READs = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]))

sample.names <- sapply(strsplit(basename(fnFs), ".fastq"), `[`, 1)

# filter the reads
p1 <- plotQualityProfile(fnFs[1:2])
p2 <- plotQualityProfile(fnFs.cut[1:2])
ggsave(paste0(output_bac, "quality_before_cut.pdf"), p1)
ggsave(paste0(output_bac, "quality_after_cut.pdf"), p2)

filtFs <- file.path(paste0(path1, "_filtered"),
                    paste0(sample.names, "_filtF.fastq.gz"))
names(filtFs) <- sample.names

saveRDS(filtFs, paste0(output_bac, "filtFs.rds"))

out <- filterAndTrim(fnFs.cut, filtFs, maxLen = 540, minLen = 200,
                     maxN=0, maxEE = 2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=40)

write.table(out, paste0(output_bac, "filter.txt"), quote = F, sep = "\t")
p3 <- plotQualityProfile(filtFs[1:2])
ggsave(paste0(output_bac, "quality_after_filter.pdf"), p3)

# error learning
filtFs <- filtFs[file.exists(filtFs)]
errF <- learnErrors(filtFs, multithread=40, MAX_CONSIST = 20)

saveRDS(errF, paste0(output_bac, "errF.rds"))

p1 <- plotErrors(errF, nominalQ = TRUE)
ggsave(paste0(output_bac, "errorF.pdf"), p1)

filtFs <- derepFastq(filtFs, n = 1e+06)

dadaFs <- dada(filtFs, err = errF, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32,
               multithread = 40, pool = TRUE)

saveRDS(dadaFs, paste0(output_bac, "dadaFs.rds"))
