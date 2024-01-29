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
samples_all = read.csv2('work/results_21_10_22/Metadata_DA.tsv',
                        header=TRUE, check.names = F, sep='\t', row.names = 1)

# otu tables
otu_mat = read.csv2('work/results_21_10_22/OTU_new/all_OTU_frequency.tsv',
                    header=TRUE, check.names = F, sep='\t')


#rename otu_table
sample_id = colnames(otu_mat)
id_names = c()
for (i in sample_id){
  id_names = c(id_names, substr(i, 2, 6))
}
colnames(otu_mat) = id_names

### taxonomy table
tax_mat = read.csv2('work/results_21_10_22/OTU_new/all_phylogeny.tsv',
                    header=TRUE, check.names = F, sep='\t')
tax_mat[1] = paste("OTU_", seq(nrow(tax_mat)), sep = '')

# filter metadata with IDs 
sample_id <- colnames(otu_mat)
samples_df = data.frame()
for (i in sample_id){
  print(i)
  row = samples_all[samples_all$Sample_ID == i, ]
  samples_df = rbind(samples_df, row)
}


OTU_name = tax_mat[1]
colnames(OTU_name) = 'OTU'


rownames(tax_mat) = OTU_name$OTU
tax_mat[1] = NULL
rownames(otu_mat) = OTU_name$OTU
otu_mat[1] = NULL


# rename rows in metatable as Sample.common.name
# rownames(samples_df) = samples_df$Sample.common.name
# samples_df$Sample_ID = NULL

# create phyloseq object
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)



####
carbom <- phyloseq(OTU, TAX, samples)
carbom
sample_sums(carbom)



### subset for soil
soil_s = subset(samples, (samples$Expedition.station.number == "76" & samples$Expedition.type == "Marine"))
soil_s$Source.sample.ID = as.character(soil_s$Source.sample.ID)

### soil and animal control
soil_control = subset(samples, (samples$Expedition.station.number == "76" & samples$Sample_type == "control"))

soil_all = subset(samples, samples$Expedition.station.number == "76" &
                    samples$Sample_type != "")

animal_control = soil_control


### subset for water
water_s = subset(samples, (samples$Expedition.station.number == "Technopark" & samples$Sample_type == "sample"))
water_s$Source.sample.ID = as.character(water_s$Source.sample.ID)

### water control
water_control = subset(samples, (samples$Expedition.station.number == "Technopark" & samples$Sample_type == "control"))
water_all = subset(samples, samples$Expedition.station.number == "Technopark" &
                     samples$Sample_type != "")

### all_control
all_samples = subset(samples, (samples$Sample_type == "sample"))
all_control = subset(samples, (samples$Sample_type == "control"))
all_data = subset(samples, (samples$Sample_type != ""))


### subset for animals
animal_s = subset(samples, (samples$Expedition.station.number == "M2" & samples$Sample_type == "sample"))
animal_s$Source.sample.ID = as.character(animal_s$Source.sample.ID)
animal_all = subset(samples, (samples$Expedition.station.number == "M2" &
                                samples$Sample_type != "") |
                      (samples$Expedition.station.number == "76" & 
                         samples$Sample_type == "control"))


current_samples <- all_control

carbom <- phyloseq(OTU, TAX, current_samples)
carbom
sample_variables(carbom)
sample_sums(carbom)
carbom


### Transform data to relative
percentFraction = function(x) (x/sum(x)*100)
carbom = transform_sample_counts(carbom, percentFraction)
sample_sums(carbom)

carbom_abund <- filter_taxa(carbom, function(x) sum(x > 1) > 0, TRUE)
sample_sums(carbom_abund)
carbom_abund <- transform_sample_counts(carbom_abund, percentFraction)
sample_sums(carbom_abund)
carbom_abund

?plot_bar

tax_table_df = as.data.frame(tax_table(carbom_abund))


orders_01 <- tax_table_df[["Order"]]

orders_02 <- tax_table_df[["Order"]]

full_orders <- c(orders_01, orders_02)
unique_orders_2 <- unique(full_orders)
# Print the list of unique orders

print(unique_orders_2)
print(length(unique_orders_2))

custom_colors_64 <- c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#808080", "#911eb4", "#46f0f0", "#f032e6",
                      "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9A6324", "#fffac8", "#800000", "#aaffc3",
                      "#808000", "#ffd8b1", "#000075", "#f58231", "#4B0082", "#7F00FF", "#ADFF2F", "#FFD700",
                      "#FFA07A", "#20B2AA", "#9370DB", "#00FA9A", "#7B68EE", "#D8BFD8", "#48D1CC", "#87CEFA",
                      "#FF4500", "#2E8B57", "#8A2BE2", "#5F9EA0", "#DC143C", "#D2691E", "#FFDEAD", "#98FB98",
                      "#00FFFF", "#FFC0CB", "#FFA500", "#C71585", "#32CD32", "#FF6347", "#4169E1", "#8B008B",
                      "#556B2F", "#00CED1", "#6A5ACD", "#BA55D3", "#CD5C5C", "#B0E0E6", "#F08080", "#4682B4",
                      "#7FFF00", "#40E0D0", "#EE82EE", "#FF69B4", "#D2B48C", "#66CDAA", "#FA8072", "#00BFFF")




orders_1 <- tax_table_df[["Order"]]

orders_2 <- tax_table_df[["Order"]]

full_orders <- c(orders_1, orders_2)
unique_orders <- unique(full_orders)
# Print the list of unique orders

print(unique_orders)
print(length(unique_orders))

order_colors <- setNames(custom_colors_64, unique_orders_2)
print(order_colors)




custom_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#7f7f7f",
                   "#8c564b", "#e377c2", "#9467bd", "#bcbd22", "#17becf",
                   "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
                   "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
                   "#c5e8a7", "#e8abab", "#e8c5da", "#fac3ac")

custom_panel_names <- c("76" = "Soil",
                        "M2" = "Gut",
                        "Technopark" = "Water")





# Custom labelling function
custom_labeller <- function(variable, value) {
  return(custom_panel_names[value])
}

### plot bars
p = plot_bar(carbom_abund, fill = "Order") + 
  geom_bar(aes(color=Family, fill=Order), stat="identity", position="stack", color=NA) + labs(x='', y = '') + 
  theme_classic() + 
  theme(legend.title = element_text(size=14), legend.text = element_text(size=12),
        axis.text = element_text(size = 11), axis.title = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1), axis.text.y = element_text(size=12)) + 
  facet_wrap(~ Expedition.station.number, scales="free_x", nrow = 3, labeller = custom_labeller) + 
  scale_fill_manual(values = order_colors)
p + guides(fill = guide_legend(ncol = 2))

### plot bars
plot_bar(carbom_abund, fill = "Family") + 
  geom_bar(aes(color=family, fill=Family), stat="identity", position="stack") + labs(x='', y = '') +
  theme(legend.title = element_text(size=16), legend.text = element_text(size=12),
        axis.text = element_text(size = 12), axis.title = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1), axis.text.y = element_text(size=12))

### plot bars
plot_bar(carbom_abund, fill = "Genus") + 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + labs(x='', y = '') +
  theme(legend.title = element_text(size=16), legend.text = element_text(size=12),
        axis.text = element_text(size = 12), axis.title = element_blank(),
        axis.text.x = element_text(angle = 65, hjust = 1), axis.text.y = element_text(size=12))


### 

otu_df = read.csv2("work/results_21_10_22/OTU_grouped_by_control_all_vs_all/OTUs_decontam_either_0.5_grouped_by_control_all_vs_all.tsv",
                   header=TRUE, check.names = F, sep='\t')

otu_mat_filtered <- as.matrix(otu_df)
OTU_f <- otu_table(otu_mat_filtered, taxa_are_rows = TRUE)

current_samples <- animal_s

carbom <- phyloseq(OTU_f, TAX, current_samples)
carbom
sample_variables(carbom)
sample_sums(carbom)
carbom


####

#### Ordination plot

### Ordination: do multivariate analysis based on Bray-Curtis distance and NMDS ordination.
### Normalize number of reads in each sample using median sequencing depth
### Normalize number of reads in each sample using median sequencing depth
current_samples <- animal_s

carbom <- phyloseq(OTU, TAX, current_samples)
carbom
sample_variables(carbom)
sample_sums(carbom)
carbom

distance_matrix <- phyloseq::distance(carbom, method = "bray")
metadata_table <- data.frame(sample_data(carbom))
permanova_result <- vegan::adonis2(distance_matrix ~ Kits_type, 
                                  data = metadata_table, permutations = 1000)
print(permanova_result)

permanova_result$`Pr(>F)`

library(reshape2)
long_format <- melt(permanova_result, id.vars = "Pr(>F)", 
                    variable.name = "Category", value.name = "PValue")






metadf <- data.frame(sample_data(carbom))

unifrac.dist <- UniFrac(carbom, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

permanova <- adonis(unifrac.dist ~ scientific_name, data = metadf)

permanova



current_samples <- water_s

carbom <- phyloseq(OTU_f, TAX, current_samples)
carbom
sample_variables(carbom)
sample_sums(carbom)
carbom


carbom = transform_sample_counts(carbom, percentFraction)
sample_sums(carbom)

carbom_abund <- filter_taxa(carbom, function(x) sum(x > 0.1) > 0, TRUE)
sample_sums(carbom_abund)
carbom_abund <- transform_sample_counts(carbom_abund, percentFraction)
sample_sums(carbom_abund)
carbom_abund

set.seed(1)
carbom.ord <- ordinate(carbom_abund, "NMDS", "bray", k=3)

sample_variables(carbom_abund)

ord_p <- plot_ordination(carbom_abund, carbom.ord, axes=c(1, 2), type="Sample.common.name", 
                color="Kits_type") + 
  geom_point(size=3) + scale_shape_manual(values=c(15, 16, 17, 18, 3, 4, 8, 9, 10)) + theme_light() + 
  theme(legend.title = element_text(size=12), legend.text = element_text(size=12),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 16))


library(factoextra)
library(ggforce)

ord_p + ggforce::geom_mark_ellipse(aes(fill = Kits_type,
                                       color = Kits_type)) +
  scale_x_continuous(limits = c(-1.2, 0.65)) + scale_y_continuous(limits = c(-0.72, 0.72))
  
  
  
plot_ordination(carbom_abund, carbom.ord, axes=c(3, 2), type="Sample.common.name", 
                color="Kits_type") + 
  geom_point(size=3) + scale_shape_manual(values=c(15, 16, 17, 18, 3, 4, 8, 9, 10)) + theme_light() + 
  theme(legend.title = element_text(size=12), legend.text = element_text(size=12),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 16))


plot_ordination(carbom_abund, carbom.ord, axes=c(3, 2), type="Sample.common.name", 
                color="Kits_type", shape = "Sample.type") + 
  geom_point(size=3) + scale_shape_manual(values=c(15, 16, 17, 18, 3, 4, 8, 9, 10)) + theme_light() + 
  theme(legend.title = element_text(size=12), legend.text = element_text(size=12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 16))



plot_ordination(carbom_abund, carbom.ord, axes=c(3, 2), type="Sample.common.name", 
                color="Метод.выделения", shape="Sample.type") + 
  geom_point(size=3) + scale_shape_manual(values=c(16, 17, 15, 18)) + theme_light() + 
  theme(legend.title = element_text(size=18), legend.text = element_text(size=16),
        axis.text.x = element_text(size = 12), 
        axis.title = element_text(size = 16))

set.seed(1)
carbom.ord <- ordinate(carbom_abund, "NMDS", "bray", weighted=T, k=3)
ordu = plot_ordination(carbom_abund, carbom.ord, axes=c(1, 2), type="Sample.common.name", 
                color="Kits_type") + 
  geom_point(size=3) + scale_shape_manual(values=c(15, 16, 17, 18, 3, 4, 8, 9, 10)) + theme_light() + 
  theme(legend.title = element_text(size=12), legend.text = element_text(size=12),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 16))

ordu + stat_ellipse()



legend(-5.5, 2.5, legend=c("2w","8w","1yr"), col=c("green","red","blue"), pch=20)

color_list <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
kit_list = unique(as.vector(soil_s$Kits_type))

plot(carbom.ord, type="n", main="Bray-Curtis",
     xlim = c(-0.5, 0.5), ylim = c(-0.7, 0.7))
legend(x = 1, y = 0.8, legend = kit_list,
       col = color_list,
       pch = 1,                # Plot symbols
       lty = 1,                # Line types
       cex = 1)
#Add an ellipse for 2w
for (i in seq(1:length(kit_list))){
  ordiellipse(carbom.ord, groups=soil_s$Kits_type, 
              display="sites", kind="se", conf=0.99, 
              label=FALSE, col=color_list[i], draw="polygon", 
              alpha=200, show.groups = kit_list[i], border=FALSE)
}

?legend()




library("DESeq2")
packageVersion("DESeq2")



log_otu_matrix <- log2(otu_mat_filtered + 1)


OTU_log <- otu_table(log_otu_matrix, taxa_are_rows = TRUE)

current_samples <- animal_s
carbom <- phyloseq(OTU_log, TAX, current_samples)
carbom
sample_variables(carbom)
sample_sums(carbom)
carbom


carbom = transform_sample_counts(carbom, percentFraction)
sample_sums(carbom)

carbom_abund <- filter_taxa(carbom, function(x) sum(x > 0.1) > 0, TRUE)
sample_sums(carbom_abund)
carbom_abund <- transform_sample_counts(carbom_abund, percentFraction)
sample_sums(carbom_abund)
carbom_abund


kit_dep = phyloseq_to_deseq2(carbom, ~ Kits_type)
kit_dep = DESeq(kit_dep, test="Wald", fitType="parametric")

res = results(kit_dep, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = as.data.frame(sigtab)

sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(carbom)[rownames(sigtab), ], "matrix"))
head(sigtab)





library("vegan")

metadf <- data.frame(sample_data(carbom_abund))
unifrac.dist <- UniFrac(carbom_abund, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)



### contamination
otu_df = as.data.frame(otu_table(carbom))
contam_relative = contam_otu

for (i in seq(1:length(contam_otu))){
  total = sum(otu_df[i])
  contam_relative[i] = round(contam_otu[i]/total, 4) * 100
  print(i)
}


otu_mat_contam <- as.matrix(contam_relative)
OTU_c <- otu_table(otu_mat_contam, taxa_are_rows = TRUE)

current_samples <- animal_s

carbom <- phyloseq(OTU_c, TAX, current_samples)
carbom
sample_variables(carbom)
sample_sums(carbom)
carbom

### Filter data

carbom_abund <- filter_taxa(carbom, function(x) sum(x > 0.1) > 0, TRUE)
sample_sums(carbom_abund)
sample_sums(carbom_abund)
carbom_abund


### plot bars
plot_bar(carbom_abund, fill = "Genus") + 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + labs(x='', y = '') + 
  scale_y_continuous(limits = c(0, 100)) + 
  theme(legend.title = element_text(size=16), legend.text = element_text(size=12),
        axis.text = element_text(size = 12), axis.title = element_blank(), panel.background = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 65, hjust = 1), axis.text.y = element_text(size=12))

##########
#################
##########################
#######################################

conamin_value = sample_sums(carbom)

contam_mean = rep(0, each = 24)
for (j in rep(0:23)){
  i = j * 3
  start = i+1
  print(start)
  end = i+3
  contam_mean[j+1] = mean(conamin_value[start:end])
}

contam_stdev = rep(0, each = 24)
for (j in rep(0:23)){
  i = j * 3
  start = i+1
  print(start)
  end = i+3
  contam_stdev[j+1] = sd(conamin_value[start:end])
}

df_contam_dev <-  data.frame(Contamin_mean = contam_mean, Contamin_sd = contam_stdev,
                             Type_of_kit = rep(c("Stool", "Microbiome", "Fecal", "B&T", "PwrSoil", 
                                             "PureLink", "Monarch", "Soil"), 3),
                    Sample.type = rep(c('Soil', 'Water', 'Organism'), each=8))

p = ggplot(df_contam_dev, aes(x = Type_of_kit, y = Contamin_mean, fill = Sample.type)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) +
  geom_errorbar(aes(ymin = Contamin_mean - Contamin_sd, ymax = Contamin_mean + Contamin_sd), 
                width = 0.2, position = position_dodge(0.6)) +
  theme(legend.title = element_text(size=16), legend.text = element_text(size=12),
        axis.text = element_text(size = 12), axis.title = element_blank(), panel.background = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 65, hjust = 1)) +
  labs(x = "Kits", y = "Mean Value") +
  theme(legend.position = "right")

p = ggplot(df_contam_dev, aes(x = Type_of_kit, y = Contamin_mean, fill = Sample.type)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) +
  geom_errorbar(aes(ymin = Contamin_mean - Contamin_sd, ymax = Contamin_mean + Contamin_sd), 
                width = 0.2, position = position_dodge(0.6)) +
  theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size=12), legend.title = element_text(size=16),
        axis.text = element_text(size = 12), axis.title = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"))

newOrder = c("Stool", "Microbiome", "Fecal", "B&T", "PwrSoil", 
             "PureLink", "Monarch", "Soil")
p$data$Type_of_kit = as.character(p$data$Type_of_kit)
p$data$Type_of_kit <- factor(p$data$Type_of_kit, levels=newOrder)
p + facet_wrap( ~ Sample.type) + theme(legend.title = element_text(size=16), legend.text = element_text(size=12),
                                       axis.text = element_text(size = 12), axis.title = element_blank(), panel.background = element_blank(),
                                       axis.ticks = element_line(color = "black"),
                                       axis.line = element_line(color = "black"),
                                       axis.text.x = element_text(angle = 65, hjust = 1))






