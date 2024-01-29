#by Viktor Mamontov, 2022
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")

carbom
head(sample_data(carbom))

df <- as.data.frame(sample_data(carbom))
df$LibrarySize <- sample_sums(carbom)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_type)) + geom_point()


#Identify Contaminants - Frequency
?isContaminant
contamdf.freq <- isContaminant(carbom, method="frequency", conc=sample_sums(carbom))
head(contamdf.freq)
table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant))

plot_frequency(carbom, taxa_names(carbom)[c(2,5,127)], conc=sample_sums(carbom)) + 
  xlab("DNA Concentration")

set.seed(1)
plot_frequency(carbom, taxa_names(ps)[sample(which(contamdf.freq$contaminant),3)], conc="quant_reading") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

carbom.noncontam <- prune_taxa(!contamdf.freq$contaminant, carbom)
carbom.noncontam
carbom <- carbom.noncontam
carbom
