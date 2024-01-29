#by Viktor Mamontov, 2022
library("mia"); packageVersion("mia");
library("ANCOMBC"); packageVersion("ANCOMBC");

set.seed(1)

current_samples <- soil_s

carbom <- phyloseq(OTU_f, TAX, current_samples)
carbom

sample_data(carbom)$Sample.type[1]
print(sample_data(carbom)$Kits_type)

newOrder <- c("Median","Stool", "Microbiome", "Fecal", "B&T", 
              "PowerSoil", "PureLink", "Monarch", "Soil")

sample_data(carbom)$Kits_type <- factor(sample_data(carbom)$Kits_type, levels = newOrder)
type(sample_data(carbom)$Kits_type)
head(sample_data(carbom)$Kits_type, 10)


exp_data <- makeTreeSummarizedExperimentFromPhyloseq(carbom)

print(exp_data)
exp_data[["Kits_type"]]


?ancombc2
?SummarizedExperiment::assay

out <- ancombc2(data = exp_data, assay_name = "counts", tax_level = "Order",
                fix_formula = "Kits_type", rand_formula = NULL,
                p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                alpha = 0.05, n_cl = 4, verbose = TRUE,
                global = FALSE, pairwise = FALSE, 
                dunnet = FALSE, trend = FALSE,
                iter_control = list(tol = 1e-5, max_iter = 20, 
                                    verbose = FALSE),
                em_control = list(tol = 1e-5, max_iter = 100),
                lme_control = NULL, mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                trend_control = NULL)

res = out$res
res_global = out$res_global

write.table(res,
            paste0("work/results_21_10_22/tables/ANCOMBC_water_order.tsv"), 
            sep = "\t", row.names = FALSE, col.names = TRUE)


res = read.csv2("work/results_21_10_22/tables/ANCOMBC_animal_order.tsv", header = T,
                sep = "\t")

lfc <- res[,c(1,3:10)]

columns <- colnames(lfc[,2:9])
for (i in columns){
  print(i)
  lfc[,as.character(i)] = as.numeric(lfc[,as.character(i)])
}

p_val <- res[, c(1,39:46)]
columns <- colnames(p_val[,2:9])
for (i in columns){
  print(i)
  p_val[,as.character(i)] = as.numeric(p_val[,as.character(i)])
}


library(tidyverse)
long_lfc <- gather(lfc, key = "Kit", value = "LFC", -taxon)

long_p_val <- gather(p_val, key = "Kit_2", value = "p-value", -taxon)
long_p_val$Kit_2 = NULL

merged_data <- cbind(long_lfc, long_p_val$`p-value`)
colnames(merged_data) <- c("Taxon", "Kit", "LFC", "p-value")

soil_df <- checkTaxon(merged_data, pval = 0.05)

water_df <- checkTaxon(merged_data, pval = 0.05)

animal_df <- checkTaxon(merged_data, pval = 0.05)

checkTaxon <- function(merged_data, pval){
  check = FALSE
  for (i in seq(merged_data$"p-value")){
    if (merged_data[i, 4] < pval){
      if (check) {
        df = rbind(df, merged_data[i,])
      }
      else {
        df <- merged_data[i,]
        check = TRUE
      }
    }
  }
  return(df)
}



inter_s_w <- intersect(soil_df$Taxon, water_df$Taxon)
inter_s_w
inter_s_a <- intersect(soil_df$Taxon, animal_df$Taxon)
inter_s_a
inter_a_w <- intersect(animal_df$Taxon, water_df$Taxon)
inter_a_w


new_soil_df <- getIntersect(soil_df, inter_s_w, inter_s_a)
new_water_df <- getIntersect(water_df, inter_s_w, inter_a_w)
new_animal_df <- getIntersect(animal_df, inter_s_a, inter_a_w)


getIntersect <- function(df, inter_1, inter_2){
  check = FALSE
  for (i in seq(soil_df$"p-value")){
    if (df[i, 1] %in% inter_1 | df[i, 1] %in% inter_1){
      if (check) {
        new_df = rbind(new_df, df[i,])
      }
      else {
        new_df <- df[i,]
        check = TRUE
      }
    }
  }
  return(new_df)
}

library(dplyr)

?inner_join()

soil_water_df <- dplyr::inner_join(soil_df, water_df, by = c("Taxon", "Kit"),
                                   suffix = c(".soil", ".water"))
soil_animal_df <- dplyr::inner_join(soil_df, animal_df, by = c("Taxon", "Kit"),
                                    suffix = c(".soil", ".animal"))
water_animal_df <- dplyr::inner_join(water_df, animal_df, by = c("Taxon", "Kit"),
                                     suffix = c(".water", ".animal"))


write.table(water_animal_df,
            paste0("work/results_21_10_22/DA_tables/water_animal_class_intersection.tsv"), 
            sep = "\t", row.names = FALSE, col.names = TRUE)


soil_water_animal_df <- dplyr::inner_join(soil_water_df, water_animal_df, by = c("Taxon", "Kit"),
                                          suffix = c(".water", ".animal"))

soil_water_animal_df$LFC.water.water = NULL
soil_water_animal_df$"p-value.water.water" = NULL
colnames(soil_water_animal_df) <- c("Taxon", "Kit", "LFC.soil", "p-value.soil",
                                    "LFC.water", "p-value.water", "LFC.animal",
                                    "p-value.animal")


soil_list <- paste(c(soil_df$Taxon), c(soil_df$Kit))
water_list <- paste(c(water_df$Taxon), c(water_df$Kit))
animal_list <- paste(c(animal_df$Taxon), c(animal_df$Kit))

# install.packages("ggVennDiagram")
library(ggVennDiagram)
# install.packages("ggplot2")
library(ggplot2)

x <- list(Soil = soil_list, Water = water_list, Gut = animal_list)
ggVennDiagram(x) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + theme(legend.position = "none")

library("ggvenn")
ggvenn(x, fill_color = c("#957777","#6C9BCF", "#9CA777"), stroke_size = 1)
