#by Viktor Mamontov, 2022
library("phyloseq")
library("ggplot2")
library("vegan"); packageVersion("vegan")
library("ape")
library('ggpubr')
library(scales)
library(ggradar)

library(dada2); packageVersion("dada2")


carbom
tab <- as.matrix(as.data.frame(otu_table(carbom)))

rarecurve(t(tab), step=500, cex=0.5, col = 'steelblue4', label = FALSE)


p = ggplot(mydf, aes(x = Contamination, y = Shennon, color=Type_of_kit, shape = Type_of_kit, alpha=Type_of_kit)) + 
  geom_point(size =6) + 
  scale_fill_manual("Type_of_kit", labels = c("Stool", "Microbiome", "Fecal", "B&T", "PwrSoil", 
                                              "PureLink", "Monarch", "Soil"),
                    values = c(0.3, 0.3, 1, 1, 1, 0.25, 0.25, 0.25)) +
  scale_alpha_manual(name = "Type_of_kit", values = c(0.3, 0.3, 1, 1, 1, 0.25, 0.25, 0.25))
newOrder = c("Stool", "Microbiome", "Fecal", "B&T", "PwrSoil", 
             "PureLink", "Monarch", "Soil")
p$data$Type_of_kit = as.character(p$data$Type_of_kit)
p$data$Type_of_kit <- factor(p$data$Type_of_kit, levels=newOrder)
p + xlab("Contamination, %") + ylab("Shannon") +  
  scale_shape_manual(values=c(15, 16, 17, 18, 4, 10, 8, 9, 13)) + theme_classic(base_size = 16)





sample_table = t(otu_table(carbom))
head(sample_table)
?specnumber
S <- specnumber(sample_table) # observed number of species
(raremax <- min(rowSums(sample_table)))
?rarefy
Srare <- rarefy(sample_table, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", col = "turquoise4")
abline(0, 0.5, col = 'red')

abline(0, 1, col = 'orange')
rarecurve(sample_table, step = 20, sample = raremax, col = "blue", cex = 0.6, label = TRUE)

abline(v = 17000, col = 'red')
rarecurve(sample_table, step = 20, sample = raremax, col = "blue", cex = 0.6, 
          xlim = c(0,21000), ylim = c(0,150))
plot(S, Srare, xlab = "Observed No. of Species", 
     ylab = "Rarefied No. of Species", col = "turquoise4", ylim = c(0,15), xlim = c(0,30))
text(S, Srare, labels = row.names(sample_table),
     cex = 0.8, pos = 4, col = "lightblue4")
?text



??ggrare

ggrare(carbom, step = 50, label = 'Sample', color = 'Sample')

samples_all <- subset(samples, samples$Sample.type != "control") 

current_samples <- samples_all

carbom <- phyloseq(OTU_f, TAX, current_samples)
carbom
sample_variables(carbom)
sample_sums(carbom)
carbom

?rarefy_even_depth

# rarefy without replacement
set.seed(1)
carbom.rarefied = rarefy_even_depth(carbom, rngseed=1, 
                                    sample.size=0.9*min(sample_sums(carbom)), replace=F)


?theme_light()

### alpha diversity
p = plot_richness(carbom, x="Kits_type", color = "Kits_type", shape="Sample.type",
              measures=c("Shannon")) + scale_shape_manual(values=c(18, 16, 17)) +
  geom_point(size=6, alpha=0.6) + 
  theme_minimal(base_size = 14) + theme(axis.text.x = element_text(angle = 65, hjust = 1))

newOrder = c("Stool", "Microbiome", "Fecal", "B&T", "PowerSoil", 
             "PureLink", "Monarch", "Soil")

#p$data$Kits_type <- soil_s$Kits_type
#print(soil_s$Kits_type)
#print(p$data$Kits_type)
p$data$Kits_type <- as.character(p$data$Kits_type)
print(p$data$Kits_type)
p$data$Kits_type <- factor(p$data$Kits_type, levels=newOrder)
p

print(soil_s$Kits_type)


plot_richness(carbom.rarefied, x="Extraction.kit", color = "Extraction_method",
              measures=c("Observed", "Shannon")) + 
  geom_boxplot() + theme_classic(base_size = 14)


### estimate richness and p-value
a_my_comparisons <- list( c("Blood&Tissue", "Microbiome"), c("Microbiome", "Saponin"), 
                          c("Saponin", "Stool"), c("Blood&Tissue","Stool"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                   symbols = c("****", "***", "**", "*", "ns"))

plot_richness(carbom.rarefied, x="Extraction_method", color = "Extraction_method",
              measures=c("Observed", "Shannon"))+
  geom_boxplot()+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+
  stat_compare_means(method = "wilcox.test", 
                     comparisons = a_my_comparisons, label = "p.signif", symnum.args = symnum.args) + 
  theme_classic(base_size = 14)



rich = estimate_richness(carbom.rarefied)
head(rich)
temp_stat <- pairwise.wilcox.test(rich$Observed, sample_data(carbom.rarefied)$Extraction_method)
df_stat <- data.frame(temp_stat$p.value)
write.table(otu_mat_contam, 'work/results_21_10_22/OTU_contam_either_05_realtive_values', col.names = TRUE,
            row.names = FALSE, sep = '\t')
write.table()

my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
             "#E16462FF", "#FCA636FF", "#F0F921FF")

ggballoonplot(df_stat, size.range = c(3,10), fill = 'value') + 
  scale_fill_gradientn(colors = my_cols)

df_stat <- data.frame(expand.grid(dimnames(temp_stat)),array(temp_stat))
df = NULL

### bars for control


########################
########################
########################

p = ggplot(mydf, aes(x = Type_of_kit, y = Shannon, color=Type_of_kit)) + 
  scale_shape_manual(values=c(18, 16, 17)) + scale_color_manual(c()) +
  geom_point(size=6, alpha=0.7) + facet_wrap(~ Sample.type, ncol = 4) + 
  theme_light(base_size = 14) + theme(axis.text.x = element_text(angle = 65, hjust = 1))

newOrder = c("Stool", "Microbiome", "Fecal", "B&T", "PwrSoil", 
             "PureLink", "Monarch", "Soil")
p$data$Type_of_kit = as.character(p$data$Type_of_kit)
p$data$Type_of_kit <- factor(p$data$Type_of_kit, levels=newOrder)
p + xlab("Kits") + ylab("Shannon")



write.table(mydf, 
            paste0(output_dir,"/Shannon_data.tsv"), 
            sep = "\t", row.names = FALSE, col.names = TRUE)




D = diversity(carbom.rarefied, "shannon")


contamin = sample_sums(carbom)
contamin

control_s = sample_data(carbom)
control_s$Sample_ID
data = control_s$Sample_ID

?diversity


alpha = c()
for (i in data){
  current_control = subset(samples, samples$Sample_ID == as.character(i))
  carbom_test <- phyloseq(OTU, TAX, current_control)
  otu_samples = data.frame(otu_table(carbom_test))
  freq_data = as.matrix(otu_samples[,1])
  alpha = c(alpha, diversity(freq_data[,1], "shannon"))
  print(i)
}

print(alpha)

alpha_stool = alpha
alpha_microb = alpha
alpha_fecal = alpha
alpha_pwrsoil = alpha

anim_alpha_stool = alpha
anim_alpha_fecal = alpha
anim_alpha_pwrsoil = alpha

length(alpha)
med_alpha = c()
for (i in seq(0, length(alpha)/3 - 1)){
  j = i * 3
  group = c(alpha[j+1], alpha[j+2], alpha[j+3])
  med_alpha = c(med_alpha, median(group))
}
alpha
med_alpha

Sample.type = rep(c('Soil', 'Water', 'Organism'), each=24)
print(Sample.type)


carbom_cont = phyloseq(OTU_c, TAX, current_samples)
contamin = sample_sums(carbom_cont)

i = 0
print(alpha[i+2+3+4])

alpha_mean = rep(0, each = 72)
for (j in rep(0:23)){
  i = j * 3
  start = i+1
  print(start)
  end = i+3
  alpha_mean[i+1] = mean(alpha[start:end])
  alpha_mean[i+2] = mean(alpha[start:end])
  alpha_mean[i+3] = mean(alpha[start:end])
}


mydf <-  data.frame(Sample_ID = data, Shannon = alpha, Shannon_mean = alpha_mean, Contamination = contamin, 
                    Type_of_kit = c("Stool", "Stool", "Stool", "Microbiome", "Microbiome", "Microbiome",
                                    "Fecal", "Fecal", "Fecal", "B&T", "B&T", "B&T",
                                    "PwrSoil", "PwrSoil", "PwrSoil", "PureLink", "PureLink", "PureLink",
                                    "Monarch", "Monarch", "Monarch", "Soil", "Soil", "Soil"),
                    Sample.type = rep(c('Soil', 'Water', 'Organism'), each=24))
p = ggplot(mydf, aes(x = Contamination, y = Shennon, color=Type_of_kit, shape = Type_of_kit, alpha=Type_of_kit)) + 
  geom_point(size =6) + 
  scale_fill_manual("Type_of_kit", labels = c("Stool", "Microbiome", "Fecal", "B&T", "PwrSoil", 
                                             "PureLink", "Monarch", "Soil"),
                    values = c(0.3, 0.3, 1, 1, 1, 0.25, 0.25, 0.25)) +
  scale_alpha_manual(name = "Type_of_kit", values = c(0.3, 0.3, 1, 1, 1, 0.25, 0.25, 0.25))
newOrder = c("Stool", "Microbiome", "Fecal", "B&T", "PwrSoil", 
             "PureLink", "Monarch", "Soil")
p$data$Type_of_kit = as.character(p$data$Type_of_kit)
p$data$Type_of_kit <- factor(p$data$Type_of_kit, levels=newOrder)
p + xlab("Contamination, %") + ylab("Shannon") +  
  scale_shape_manual(values=c(15, 16, 17, 18, 4, 10, 8, 9, 13)) + theme_classic(base_size = 16)


### for separated kits
contamin = sample_sums(carbom)

mydf <- animal_pwrsoil
mydf$Shennon = anim_alpha_pwrsoil
mydf$Contamination = contamin

p = ggplot(mydf, aes(x = Contamination, y = Shennon, color=Comments, 
                     shape = Comments, alpha=Comments)) + geom_point(size = 6, alpha = 0.6)
p + xlab("Уровень контаминации, %") + ylab("Индекс Шеннона") +  
  scale_shape_manual(values=c(15, 16, 17, 18, 4, 10, 8, 9, 13)) + theme_classic(base_size = 16)

###  scale_fill_manual("Type_of_kit", labels = mydf$Comments,
#                    values = c(0.7)) +
#  scale_alpha_manual(name = "Comments", values = c(0.3, 0.3, 1, 1, 1, 0.25, 0.25, 0.25))

p$data$Type_of_kit = as.character(p$data$Type_of_kit)
p$data$Type_of_kit <- factor(p$data$Type_of_kit, levels=newOrder)





### radar plot
contamin = sample_sums(carbom)


count = 0
median_contamin = c()
for (j in contamin){
  count = count + 1
  if (count %% 3 == 0){
    median_contamin = c(median_contamin, contamin[count - 1])
  }
}



water_df <- data.frame(Shannon = median_alpha, Type_of_kit = c("Stool", "Microbiome",
                                                               "Fecal", "B&T",
                                                               "PwrSoil", "PureLink",
                                                               "Monarch", "LSBio"))
water_df$Type_of_sample = "Water"

mydf <-  data.frame(Shannon = median_alpha, 
                    Type_of_kit = c("Stool", "Microbiome",
                                    "Fecal", "B&T",
                                    "PwrSoil", "PureLink",
                                    "Monarch", "LSBio"))
mydf$Type_of_samples = "Soil"

colnames(mydf) = c('Shannon', "Type_of_kit", "Type_of_sample")

animal_df <- data.frame(Shannon = median_alpha, Type_of_kit = c("Stool", "Microbiome",
                                                                "Fecal", "B&T",
                                                                "PwrSoil", "PureLink",
                                                                "Monarch", "LSBio"))
animal_df$Type_of_sample = "Animal"

df <- rbind(mydf, water_df, animal_df)
df$Contamination <- contamin


library(ggradar)
library(tidyverse)
library(tidyquant)
library(scales)
library(corrr)

df <- read.csv2('work/results_21_10_22/kit_features.csv', header=TRUE,
             sep='\t')

for (i in seq(2,8)){
  print(i)
  df[,i] = as.numeric(df[,i]) 
}
sapply(df, class)

?rescale

df %>% mutate_each(funs(rescale), -Type) -> df_radar
df_radar %>% ggradar(base.size = 10, axis.label.size = 4, legend.text.size = 12,
                     legend.position = "right", group.line.width = 1, group.point.size = 3, 
                     fill = T, fill.alpha = 0.25)

colnames(df_radar) = c("Type", "α S", "α W", "α A",
                       "P S", "P W", "P A", "Obj")

df_radar$Type <- factor(df_radar$Type,
                        levels = c(
                          "PowerSoil", "Fecal", "B&T", "PureLink",
                          "Soil",  "Monarch", "Stool", "Microbiome"
                        ))

r <- ggradar(df_radar, base.size = 10, axis.label.size = 4, legend.text.size = 12,
        legend.position = "right", group.line.width = 1, group.point.size = 3, 
        fill = T, fill.alpha = 0.25)
r + facet_wrap(~ Type, ncol = 4) + theme_void() + 
  theme(
    strip.text = element_text(
      size = 12,
      colour = "Black",
      margin = margin(t = 5, b = 5)
    ), strip.background = element_rect(fill = "white"),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )
