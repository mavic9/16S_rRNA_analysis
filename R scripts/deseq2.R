#by Viktor Mamontov, 2022
library("DESeq2")

otu_df = read.csv2("work/results_21_10_22/OTU_grouped_by_control_all_vs_all/OTUs_decontam_either_0.5_grouped_by_control_all_vs_all.tsv",
                   header=TRUE, check.names = F, sep='\t')

otu_mat_filtered <- as.matrix(otu_df)
OTU_f <- otu_table(otu_mat_filtered, taxa_are_rows = TRUE)

current_samples <- soil_s

carbom <- phyloseq(OTU_f, TAX, current_samples)
carbom
sample_variables(carbom)
sample_sums(carbom)
carbom

newOrder <- c("Stool", "Microbiome", "Fecal", "B&T", 
              "PowerSoil", "PureLink", "Monarch", "Soil")

sample_data(carbom)$Kits_type <- factor(sample_data(carbom)$Kits_type, levels = newOrder)
head(sample_data(carbom)$Kits_type, 24)


ds = phyloseq_to_deseq2(carbom, ~ Kits_type)
ds = DESeq(ds)

alpha = 0.01
res = results(ds, contrast = c("Kits_type", "Microbiome", "B&T"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig

res_sig = cbind(as(res_sig, "data.frame"), as(tax_table(carbom)[rownames(res_sig), ], "matrix"))
ggplot(res_sig, aes(x=Order, y=log2FoldChange, color=Phylum)) +
  geom_jitter(size=3, width = 0.2) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))





