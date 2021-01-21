#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(dada2)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)
path <- args[1]
output <- args[2]

################################ for Bacteria
path1 <- paste0(path, "/Bac_forward/out")
path2 <- paste0(path, "/Bac_reverse/out")
output_bac <- paste0(output, "Bacteria/")

fnFs <- sort(list.files(path1, pattern=".fastq", full.names = TRUE))
fnRs <- sort(list.files(path2, pattern=".fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), ".fastq"), `[`, 1)

# filter the reads
filtFs <- file.path(paste0(path1, "_filtered"),
                    paste0(sample.names, "_filtF.fastq.gz"))
filtRs <- file.path(paste0(path2, "_filtered"),
                    paste0(sample.names, "_filtR.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

saveRDS(filtFs, paste0(output_bac, "filtFs.rds"))
saveRDS(filtRs, paste0(output_bac, "filtRs.rds"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,240),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=40)
write.table(out, paste0(output_bac, "filter.txt"), quote = F, sep = "\t")

# error learning
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]
errF <- learnErrors(filtFs, multithread=40, MAX_CONSIST = 20)
errR <- learnErrors(filtRs, multithread=40, MAX_CONSIST = 20)

saveRDS(errF, paste0(output_bac, "errF.rds"))
saveRDS(errR, paste0(output_bac, "errR.rds"))

p1 <- plotErrors(errF, nominalQ = TRUE)
p2 <- plotErrors(errR, nominalQ = TRUE)
ggsave(paste0(output_bac, "errorF.pdf"), p1)
ggsave(paste0(output_bac, "errorR.pdf"), p2)

filtFs <- derepFastq(filtFs, n = 1e+06)
filtRs <- derepFastq(filtRs, n = 1e+06)

dadaFs <- dada(filtFs, err = errF, multithread = 40, pool = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = 40, pool = TRUE)

saveRDS(dadaFs, paste0(output_bac, "dadaFs.rds"))
saveRDS(dadaRs, paste0(output_bac, "dadaRs.rds"))

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
saveRDS(mergers, paste0(output_bac, "mergers.rds"))

seqtab <- makeSequenceTable(mergers, orderBy = "abundance")
saveRDS(seqtab, paste0(output_bac, "seqtab.rds"))
