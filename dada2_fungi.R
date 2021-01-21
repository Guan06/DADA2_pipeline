#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(dada2)
library(ggplot2)
library(ShortRead)
library(Biostrings)

args = commandArgs(trailingOnly = TRUE)
path <- args[1]
output <- args[2]
FWD <- args[3]
REV <- args[4]

########################################### for Fungi
path1 <- paste0(path, "/Fun_forward/out")
path2 <- paste0(path, "/Fun_reverse/out")
output_fun <- paste0(output, "Fungi/")

fnFs <- sort(list.files(path1, pattern=".fastq", full.names = TRUE))
fnRs <- sort(list.files(path2, pattern=".fastq", full.names = TRUE))

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

# pre-filter the reads with Ns
fnFs.filtN <- file.path(path1, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path2, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, compress = FALSE,
              multithread = 40)

# count the number of times the primers appear
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

cutadapt <- "/opt/share/software/bin/cutadapt"
#system2(cutadapt, args = "--version")

path1.cut <- file.path(path1, "cutadapt")
if(!dir.exists(path1.cut)) dir.create(path1.cut)
path2.cut <- file.path(path2, "cutadapt")
if(!dir.exists(path2.cut)) dir.create(path2.cut)

fnFs.cut <- file.path(path1.cut, basename(fnFs))
fnRs.cut <- file.path(path2.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

# Run Cutadapt
for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags,
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i],
                               fnFs.filtN[i], fnRs.filtN[i]))
}

# count the presence of primers after cutadapt
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

cutFs <- sort(list.files(path1.cut, pattern = ".fastq", full.names = TRUE))
cutRs <- sort(list.files(path2.cut, pattern = ".fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq")[[1]][1]

sample.names <- unname(sapply(fnFs, get.sample.name))

# filter the reads
filtFs <- file.path(paste0(path1, "_filtered"),
                    paste0(sample.names, "_filtF.fastq.gz"))
filtRs <- file.path(paste0(path2, "_filtered"),
                    paste0(sample.names, "_filtR.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                     maxN = 0, maxEE = c(2, 2),
                     multithread = 40)

saveRDS(filtFs, paste0(output_fun, "filtFs.rds"))
saveRDS(filtRs, paste0(output_fun, "filtRs.rds"))

write.table(out, paste0(output_fun, "filter.txt"), quote = F, sep = "\t")

# error learning
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]

errF <- learnErrors(filtFs, multithread = 40, MAX_CONSIST = 20)
errR <- learnErrors(filtRs, multithread = 40, MAX_CONSIST = 20)

saveRDS(errF, paste0(output_fun, "errF.rds"))
saveRDS(errR, paste0(output_fun, "errR.rds"))

p1 <- plotErrors(errF, nominalQ = TRUE)
p2 <- plotErrors(errR, nominalQ = TRUE)
ggsave(paste0(output_fun, "errorF.pdf"), p1)
ggsave(paste0(output_fun, "errorR.pdf"), p2)

filtFs <- derepFastq(filtFs, n = 1e+06)
filtRs <- derepFastq(filtRs, n = 1e+06)

dadaFs <- dada(filtFs, err = errF, multithread = 40, pool = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = 40, pool = TRUE)

saveRDS(dadaFs, paste0(output_fun, "dadaFs.rds"))
saveRDS(dadaRs, paste0(output_fun, "dadaRs.rds"))

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
saveRDS(mergers, paste0(output_fun, "mergers.rds"))

seqtab <- makeSequenceTable(mergers, orderBy = "abundance")
saveRDS(seqtab, paste0(output_fun, "seqtab.rds"))
