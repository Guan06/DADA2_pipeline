#!//biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(dada2)
library(ggplot2)
library(ShortRead)
library(Biostrings)

args = commandArgs(trailingOnly = TRUE)
path <- args[1]
output <- args[2]

########################################### trim Roche 454 data
REV1 <- "ACGTCATCCCCACCTTCC"
REV2 <- "ACGTCATCCCCACCTTCT"
REV3 <- "GTCATCCCCACCTTCC"

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
REV1.orients <- allOrients(REV1)
REV2.orients <- allOrients(REV2)
REV3.orients <- allOrients(REV3)

################################ for Bacteria
path1 <- paste0(path, "/out")
output_bac <- output

fnFs <- sort(list.files(path1, pattern=".fastq", full.names = TRUE))

# count the number of times the primers appear
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

rbind(REV1.Reads = sapply(REV1.orients, primerHits, fn = fnFs[[1]]),
      REV2.Reads = sapply(REV2.orients, primerHits, fn = fnFs[[1]]),
      REV3.Reads = sapply(REV3.orients, primerHits, fn = fnFs[[1]]))

cutadapt <- "/opt/share/software/bin/cutadapt"
path.cut <- file.path(path1, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)

fnFs.cut <- file.path(path.cut, basename(fnFs))

REV1.RC <- dada2:::rc(REV1)
flag1 <- paste("-a", REV1.RC)
#flag2 <- paste("-a", REV2.RC)
#flag3 <- paste("-a", REV3.RC)

# Run Cutadapt
for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(flag1, "-o", fnFs.cut[i], fnFs[i]))
    #system2(cutadapt, args = c(flag2, "-o", fnFs.cut[i], fnFs[i]))
    #system2(cutadapt, args = c(flag3, "-o", fnFs.cut[i], fnFs[i]))
}
rbind(REV1.Reads = sapply(REV1.orients, primerHits, fn = fnFs.cut[[1]]),
      REV2.Reads = sapply(REV2.orients, primerHits, fn = fnFs.cut[[1]]),
      REV3.Reads = sapply(REV3.orients, primerHits, fn = fnFs.cut[[1]]))

sample.names <- sapply(strsplit(basename(fnFs), ".fastq"), `[`, 1)

# filter the reads
filtFs <- file.path(paste0(path1, "_filtered"),
                    paste0(sample.names, "_filtF.fastq.gz"))
names(filtFs) <- sample.names

saveRDS(filtFs, paste0(output_bac, "filtFs.rds"))

out <- filterAndTrim(fnFs.cut, filtFs, maxLen = 440,
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=40)

write.table(out, paste0(output_bac, "filter.txt"), quote = F, sep = "\t")

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
