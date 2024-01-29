#by Viktor Mamontov, 2022
### good phyloseq script
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library("tibble") 
library(RColorBrewer)
library(dplyr); packageVersion("dplyr")
library(tidyr)

# database
all_OTU_phyl = read.csv2('work/results_21_10_22/OTU_new/all_OTUs_phylogeny.tsv',
                        header=TRUE, check.names = F, sep='\t', row.names = 1)
### taxonomy table
tax_mat = read.csv2('work/results_21_10_22/OTU_new/all_phylogeny.tsv',
                    header=TRUE, check.names = T, sep='\t')
tax_mat[1] = paste("OTU_", seq(nrow(tax_mat)), sep = '')

tax_mat$OTU_seq = rownames(all_OTU_phyl)

rownames(tax_mat) = tax_mat$X
tax_mat$X = NULL

tax_mat$OTUs = tax_mat$X


add_table = cbind(data.frame(tax_mat$OTUs), data.frame(tax_mat$OTU_seq))

colnames(add_table) = c("OTUs", "OTU_seq")


OTU_control = read.csv2('work/results_21_10_22/tables/OTU_phyl_control.txt',
                        header=TRUE, check.names = F, sep='\t')




df3 = merge(OTU_control, add_table, by = "OTUs")

df3$OTU_num <- as.numeric(gsub("OTU_", "", df3$OTUs))
df_sorted <- df3[order(df3$OTU_num), ]
df_sorted$OTU_num <- NULL

write.table(df_sorted, 'work/results_21_10_22/results/OTU_control_seq.tsv', 
            sep = "\t", row.names = FALSE, col.names = TRUE)


