#by Viktor Mamontov, 2022
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")


current_samples = water_all
carbom <- phyloseq(OTU, TAX, current_samples)
carbom


sample_data(carbom)$is.neg <- sample_data(carbom)$Sample_type == "control";
sample_data(carbom)

?isContaminant
sample_variables(carbom)

contamdf.either <- isContaminant(carbom, method="either", neg="is.neg",
                                 conc="DNA.concentration..ng.mkl",
                                 threshold=0.5, normalize = T)

table(contamdf.either$contaminant)
head(which(contamdf.either$contaminant))


contamdf.prev <- isContaminant(carbom, method="prevalence", neg="is.neg", 
                               threshold=0.5, normalize = T)
table(contamdf.prev$contaminant)

head(which(contamdf.prev$contaminant))


# Make phyloseq object of presence-absence in negative controls and true samples
carbom.pa <- transform_sample_counts(carbom, function(abund) 1*(abund>0))
sample_sums(carbom.pa)
carbom.pa.neg <- prune_samples(sample_data(carbom.pa)$Sample_type == "control", carbom.pa)
carbom.pa.pos <- prune_samples(sample_data(carbom.pa)$Sample_type == "sample", carbom.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(carbom.pos=taxa_sums(carbom.pa.pos), carbom.neg=taxa_sums(carbom.pa.neg),
                    contaminant=contamdf.either$contaminant)
ggplot(data=df.pa, aes(x=carbom.neg, y=carbom.pos, color=contaminant)) + geom_point(size=3, alpha=0.4) +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


?prune_taxa

carbom.pos <- prune_samples(sample_data(carbom)$Sample_type == "sample", carbom)
carbom.noncontam <- prune_taxa(!contamdf.prev$contaminant, carbom.pos)
carbom.noncontam

noncontam_otu = data.frame(otu_table(carbom.noncontam), check.names = F)

### Transform data to relative
percentFraction = function(x) (x/sum(x)*100)
carbom.noncontam = transform_sample_counts(carbom.noncontam, percentFraction)
sample_sums(carbom.noncontam)

carbom_abund <- filter_taxa(carbom.noncontam, function(x) sum(x > 2) > 0, TRUE)
sample_sums(carbom_abund)
carbom_abund <- transform_sample_counts(carbom_abund, percentFraction)
sample_sums(carbom_abund)
carbom_abund

### plot bars
plot_bar(carbom_abund, fill = "Order") + 
  geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack") + labs(x='Samles', y = 'Abundance') +
  theme(legend.title = element_text(size=16), legend.text = element_text(size=12),
        axis.text = element_text(size = 12), axis.title = element_text(size = 24),
        axis.text.x = element_text(angle = 65, hjust = 1), axis.text.y = element_text(size=12))


?isContaminant

carbom <- carbom.noncontam
carbom

### only contamination OTUs

carbom.pos <- prune_samples(sample_data(carbom)$Sample_type == "sample", carbom)
sample_sums(carbom.pos)
otu_sums = as.vector(sample_sums(carbom.pos))
carbom.contam <- prune_taxa(contamdf.either$contaminant, carbom.pos)
carbom.contam

otu_contam = data.frame(otu_table(carbom.contam), check.names = F)
new_otu = otu_contam

for (i in seq(1, length(otu_sums))){
  for (j in seq(1, length(otu_contam[,1]))){
    new_otu[j, i] = otu_contam[j, i]/as.numeric(otu_sums[i]) * 100
  }
  print(i)
}

write.table(new_otu, "/home/viktor/Desktop/results_21_10_22/tables_decontam/otu_water_contam_relative.tsv", 
            sep = "\t", row.names = TRUE, col.names = TRUE)

otu_mat = new_otu

current_samples = water_all

# create phyloseq object
otu_mat <- as.matrix(otu_mat)

OTU_new = otu_table(otu_mat, taxa_are_rows = TRUE)
carbom <- phyloseq(OTU_new, TAX, current_samples)
carbom

carbom.contam <- carbom
sample_sums(carbom.contam)

carbom_abund <- filter_taxa(carbom.contam, function(x) sum(x > 0.5) > 0, TRUE)
sample_sums(carbom_abund)
carbom_abund

### plot bars
plot_bar(carbom_abund, fill = "Order") + 
  geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack") + labs(x='Samles', y = 'Abundance') +
  theme(legend.title = element_text(size=16), legend.text = element_text(size=12),
        axis.text = element_text(size = 12), axis.title = element_text(size = 24),
        axis.text.x = element_text(angle = 65, hjust = 1), axis.text.y = element_text(size=12))




### plot bars
plot_bar(carbom.contam, fill = "Order") + 
  geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack") + labs(x='Samles', y = 'Abundance') +
  theme(legend.title = element_text(size=16), legend.text = element_text(size=12),
        axis.text = element_text(size = 12), axis.title = element_text(size = 24),
        axis.text.x = element_text(angle = 65, hjust = 1), axis.text.y = element_text(size=12))



### Transform data to relative
percentFraction = function(x) (x/sum(x)*100)
carbom.contam = transform_sample_counts(carbom.contam, percentFraction)
sample_sums(carbom.contam)

carbom_abund <- filter_taxa(carbom.contam, function(x) sum(x > 0.1) > 0, TRUE)
sample_sums(carbom_abund)
carbom_abund <- transform_sample_counts(carbom_abund, percentFraction)
sample_sums(carbom_abund)
carbom_abund

### plot bars
plot_bar(carbom_abund, fill = "Order") + 
  geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack") + labs(x='Samles', y = 'Abundance') +
  theme(legend.title = element_text(size=16), legend.text = element_text(size=12),
        axis.text = element_text(size = 12), axis.title = element_text(size = 24),
        axis.text.x = element_text(angle = 65, hjust = 1), axis.text.y = element_text(size=12))




