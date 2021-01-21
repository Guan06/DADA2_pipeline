#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(dada2)

args = commandArgs(trailingOnly = TRUE)
path <- args[1]
output <- args[2]

path <- list.dirs(path, full.names = T, recursive = F)
path <- paste0(path, "/Fungi/")

seqtabs <- list.files(path, pattern = "seqtab.rds", full.names = TRUE)

if (length(seqtabs) < 2) {
    asv <- readRDS(seqtabs)
} else {
    asv <- mergeSequenceTables(tables = seqtabs, repeats = "sum",
                               orderBy = "abundance")
}

## remove chimeras
asv.nochime <- removeBimeraDenovo(asv, method = "consensus",
                                  multithread = 40, verbose = TRUE)
dim(asv.nochime)
sum(asv.nochime) / sum(asv)

asv <- asv.nochime
nc <- ncol(asv)
map <- data.frame(ASV_ID = paste0("ASV_", seq(1 : nc)),
                  Sequence = colnames(asv))

## assign taxonomy
unite.ref <- "./Databases/sh_general_release_dynamic_04.02.2020.fasta"
taxa <- assignTaxonomy(asv, unite.ref, multithread = 40, tryRC = TRUE)
rownames(taxa) <- map$ASV_ID[match(rownames(taxa), map$Sequence)]

write.table(taxa, paste0(output, "ASV_taxonomy_unite.txt"),
            quote = F, sep = "\t")

## re-format the ASV table and write it to file
colnames(asv) <- map$ASV_ID[match(colnames(asv), map$Sequence)]
asv <- t(asv)
stat <- data.frame(Sample_ID = colnames(asv),
                   Present = colSums(asv > 0),
                   Abundance = colSums(asv))

write.table(asv, paste0(output, "ASV_table.txt"),
            quote = F, sep = "\t")
write.table(stat, paste0(output, "ASV_lib_stat.txt"),
            quote = F, sep = "\t", row.names = F)
write.table(map, paste0(output, "ASV_map.txt"),
            quote = F, sep = "\t", row.names = F, col.names = F)
