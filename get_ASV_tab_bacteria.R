#!/biodata/dep_psl/grp_rgo/guan/bin/Rscript

library(dada2)

args = commandArgs(trailingOnly = TRUE)
path <- args[1]
output <- args[2]

path <- list.dirs(path, full.names = T, recursive = F)
path <- paste0(path, "/Bacteria/")
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
taxtrain1 <- "./Databases/silva_nr_v132_train_set.fa.gz"
taxa1 <- assignTaxonomy(asv, taxtrain1, multithread = 40)
rownames(taxa1) <- map$ASV_ID[match(rownames(taxa1), map$Sequence)]

taxtrain2 <- "./Databases/rdp_train_set_16.fa.gz"
taxa2 <- assignTaxonomy(asv, taxtrain2, multithread = 40)
rownames(taxa2) <- map$ASV_ID[match(rownames(taxa2), map$Sequence)]

taxtrain3 <- "./Databases/gg_13_8_train_set_97.fa.gz"
taxa3 <- assignTaxonomy(asv, taxtrain3, multithread = 40)
rownames(taxa3) <- map$ASV_ID[match(rownames(taxa3), map$Sequence)]

write.table(taxa1, paste0(output, "/ASV_taxonomy_silva.txt"),
            quote = F, sep = "\t")
write.table(taxa2, paste0(output, "/ASV_taxonomy_rdp.txt"),
            quote = F, sep = "\t")
write.table(taxa3, paste0(output, "/ASV_taxonomy_gg.txt"),
            quote = F, sep = "\t")

## re-format the ASV table and write it to file
colnames(asv) <- map$ASV_ID[match(colnames(asv), map$Sequence)]
asv <- t(asv)
stat <- data.frame(Sample_ID = colnames(asv),
                   Present = colSums(asv > 0),
                   Abundance = colSums(asv))

write.table(asv, paste0(output, "/ASV_table.txt"),
            quote = F, sep = "\t")
write.table(stat, paste0(output, "/ASV_lib_stat.txt"),
            quote = F, sep = "\t", row.names = F)
write.table(map, paste0(output, "/ASV_map.txt"),
            quote = F, sep = "\t", row.names = F, col.names = F)
