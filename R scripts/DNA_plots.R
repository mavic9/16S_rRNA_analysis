#by Viktor Mamontov, 2022
library("phyloseq")
library("ggplot2")
library("vegan"); packageVersion("vegan")
library("ape")
library('ggpubr')
library(scales)
library(ggradar)
library(ggbeeswarm)

my_df <- read.csv2('../Desktop/Work/Metagenomes/TMT438_all.txt', header=T,
                   sep='\t', check.names = F)
my_df$Values <- as.numeric(my_df$Values)

df1 <- read.csv2('../Desktop/Work/Metagenomes/TMT438_stats.tsv', header=T,
                 sep='\t', check.names = F)

df1$Value <- as.numeric(df1$Value)
df1$Statistics <- as.factor(df1$Statistics)
type(df1$Statistics)


df <- read.csv2('work/results_21_10_22/tables/Pure_DNA.tsv', header=T,
                sep = '\t', check.names = F)
color_list <- c("PowerSoil" = "#ff5a5f", "PowerFecal" = "#ffb400", "Stool" = "#007a87",  
                "B&T" = "#8ce071", "LSBio" = "#7b0051", "PureLink" = "#00d1c1", 
                "Microbiome" = "#ffaa91", "Monarch" = "#b4a76c")

df$`Ct(18S)-Ct(16S)` <-  as.numeric(as.vector(df$`Ct(18S)-Ct(16S)`))

df$"DNA, ng" = as.numeric(as.vector(df$`DNA, ng`))
df$`Pure DNA, %` = as.numeric(as.vector(df$`Pure DNA, %`))

df$"DNA, ng (log scale)" = as.numeric(as.vector(df$"DNA, ng (log scale)"))
df$DIN = as.numeric(as.vector(df$DIN))
df$`260/280` = as.numeric(as.vector(df$`260/280`))
df$`260/230` = as.numeric(as.vector(df$`260/230`))
df$Type <- factor(df$Type, levels=c("Soil", "Water", "Guts"))

### facet_wrap
p <- ggplot(df, aes(x = Type, y = `Pure DNA, %`, color=Type)) + 
  scale_color_manual(values = c("#6394eb", "#ee684f","#cc940a", "#688c21")) + 
  geom_point(size=6, alpha=1, show.legend = F) + facet_wrap(~ Kit, ncol = 2) + 
  ylim(0, 150) + theme_classic(base_size = 13) + 
  theme(axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, 
                                   vjust = 0.5, colour = 'black'),
        strip.background = element_blank(), 
        strip.text = element_text(size = 12, face = 'bold'))

newOrder = c("AMPure", "Evrogen", "Thermo", "NEB", "PowerSoil")
p$data$Type <- factor(p$data$Type, levels=newOrder)
p + xlab("") + ylab("Estimated DNA retention, %")

?geom_jitter()


### enrichment_1
df$Type <- factor(df$Type, levels=c("Input", "NEB Enrichment"))
df$Kit <- factor(df$Kit, levels=c("Stool", "PowerFecal", "PowerSoil"))
### facet_wrap
p <- ggplot(df, aes(x = Kit, y = `Ct(18S)-Ct(16S)`, color=Kit)) + 
  scale_shape_manual(values=c(18, 16, 17)) + 
  scale_color_manual(values = c("#f8766d", "#7cae00","#00bfc4")) +
  geom_point(size=6, alpha=1, show.legend = FALSE) + facet_wrap(~ Type, ncol = 3) + 
  theme_classic(base_size = 13) + theme(axis.text.y = element_text(colour = 'black'),
                                        axis.text.x = element_text(angle = 90, hjust = 1, 
                                                                   vjust = 0.5, colour = 'black'),
                                        strip.background = element_blank(),
                                        strip.text = element_text(size = 12, face = 'bold'))

p + xlab("") + ylab("Ct(18S)-Ct(16S)")

#enrichment_2
df$`DNA, ng` <-  as.numeric(as.vector(df$`DNA, ng`))
df$`Estimated DNA retention, %` <- as.numeric(as.vector(df$`Estimated DNA retention, %`))

df$Type <- factor(df$Type, levels=c("before NEB Enrichment", "after NEB Enrichment"))
df$Kit <- factor(df$Kit, levels=c("Stool", "PowerFecal", "PowerSoil"))
### facet_wrap
p <- ggplot(df, aes(x = Kit, y = `DNA, ng`, color=Kit)) + 
  scale_shape_manual(values=c(18, 16, 17)) + 
  scale_color_manual(values = c("#f8766d", "#7cae00","#00bfc4")) +
  geom_point(size=6, alpha=1) + facet_wrap(~ Type, ncol = 3) + 
  theme_classic(base_size = 13) + theme(axis.text.y = element_text(colour = 'black'),
                                        axis.text.x = element_text(angle = 90, hjust = 1, 
                                                                   vjust = 0.5, colour = 'black'),
                                        strip.background = element_blank(),
                                        strip.text = element_text(size = 12, face = 'bold'))

p + xlab("") + ylab("DNA, ng")

#enrichment_3
df <- read.csv2('work/results_21_10_22/results/enrichment_3.txt', header=T,
                sep = '\t', check.names = F)
df$`Estimated DNA retention, %` <- as.numeric(as.vector(df$`Estimated DNA retention, %`))

df$Type <- factor(df$Type, levels=c("before NEB Enrichment", "after NEB Enrichment"))
df$Kit <- factor(df$Kit, levels=c("Stool", "PowerFecal", "PowerSoil"))


p <- ggplot(df, aes(x = Kit, y = `Estimated DNA retention, %`, color=Kit)) + 
  scale_shape_manual(values=c(18, 16, 17)) + 
  scale_color_manual(values = c("#f8766d", "#7cae00","#00bfc4")) +
  geom_point(size=6, alpha=1) + facet_wrap(~ Type, ncol = 3) + 
  theme_classic(base_size = 13) + theme(axis.text.y = element_text(colour = 'black'),
                                        axis.text.x = element_text(angle = 90, hjust = 1, 
                                                                   vjust = 0.5, colour = 'black'),
                                        strip.background = element_blank(),
                                        strip.text = element_text(size = 12, face = 'bold'))

p + xlab("") + ylab("Estimated DNA retention, %")

### common plot
p <- ggplot(df, aes(x = Kit, y = `18S/16S`, color=Kit, shape=Type, 
                    group=interaction(Kit, Type))) + 
  scale_shape_manual(values=c(18, 16, 17)) +
  geom_point(size=3, position = position_dodge(width=0.35), alpha=0.8) +
  theme_classic(base_size = 13) + theme(axis.text.y = element_text(colour = 'black'),
                                        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1,
                                                                  colour = 'black'),
                                        strip.background = element_blank(),
                                        strip.text = element_text(size = 12, face = 'bold'))
p$data$Kit = as.character(p$data$Kit)
p$data$Kit <- factor(p$data$Kit, levels=newOrder)
p

### common plot
p <- ggplot(df, aes(x = Kit, y = `18S/16S`, color=Type, 
                    group=interaction(Kit, Type))) + 
  scale_shape_manual(values=c(18, 16, 17)) + 
  scale_color_manual(values = c("#f8766d", "#00a9ff","#00be67")) +
  geom_point(size=3, position = position_dodge(width=0.35), alpha=0.8) +
  theme_classic(base_size = 13) + theme(axis.text.y = element_text(colour = 'black'),
                                        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1,
                                                                   colour = 'black'),
                                        strip.background = element_blank(),
                                        strip.text = element_text(size = 12, face = 'bold'))

p$data$Kit = as.character(p$data$Kit)
p$data$Kit <- factor(p$data$Kit, levels=newOrder)
p

### TMT438
p <- ggplot(my_df, aes(x = Samples, y = Values, fill = Stats)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("#f8766d", "#00a9ff")) +
  theme_minimal() +
  labs(x = "Samples", y = "Values", fill = "Statistics")
p + xlab("Samples") + ylab("Number of Reads")


p <- ggplot(df1, aes(x = Sample, y = Value, fill = Statistics)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Statistics, scales = "free_y") +
  theme_minimal() +
  labs(x = "Samples", y = "Values", fill = "Statistics")
p$data$Statistics = factor(p$data$Statistics, 
                              levels = c('Total Length, bp', 'Length >1 kb',
                                         'Total Contigs', 'Contigs >1 kb', 'Contigs >5 kb',
                                         'N50', 'L50', 'Largest contig, bp',
                                         'gDNA, ng', 'Raw Reads/gDNA', 'Trimmed Reads/gDNA',
                                         'Passed Reads, %'))
p$data$Sample = factor(p$data$Sample, levels=c('TMT438', 'PTU308'))
p + theme(strip.text = element_text(size = 10, margin = margin(t = 10, b = 10), face = 'bold'))



