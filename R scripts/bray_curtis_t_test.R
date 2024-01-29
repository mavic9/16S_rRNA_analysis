

pairwise_t_test <- function(distance_matrix, grouping_factor) {
  unique_groups <- unique(grouping_factor)
  n_groups <- length(unique_groups)
  print(unique_groups)
  print(n_groups)
  # Initialize a matrix to store pairwise p-values
  p_values <- matrix(NA, nrow = n_groups, ncol = n_groups)
  rownames(p_values) <- colnames(p_values) <- unique_groups
  
  for (i in 1:(n_groups - 1)) {
    for (j in (i + 1):n_groups) {
      physeq <- carbom
      factors <- c(unique_groups[i], unique_groups[j])
      print(factors)
      factor_1  <- as.logical(sample_data(physeq)$Kits_type == factors[1])
      factor_2  <- as.logical(sample_data(physeq)$Kits_type == factors[2])
      
      group_1 = c(as.vector(distance_matrix[factor_1, factor_1]),
                  as.vector(distance_matrix[factor_2, factor_2]))
      group_1 = group_1[group_1 != 0]
      group_1 = unique(group_1)
      print('Group 1: ')
      print(length(group_1))
      group_2 = c(as.vector(distance_matrix[factor_1, factor_2]),
                  as.vector(distance_matrix[factor_2, factor_1]))
      group_2 = unique(group_2)
      print('Group 2: ')
      print(length(group_2))
      t_test_result <- t.test(group_1, group_2)
      
      p_values[i, j] <- t_test_result$p.value
      print(t_test_result$p.value)
    }
  }
  
  return(p_values)
}


# Replace 'your_phyloseq_object' with the name of your phyloseq object
current_samples <- animal_s

carbom <- phyloseq(OTU_f, TAX, current_samples)
carbom
sample_variables(carbom)
sample_sums(carbom)
carbom

physeq <- carbom
physeq
sample_data_df <- as.data.frame(sample_data(physeq))
otu_table_matrix <- as.matrix(as.data.frame(otu_table(physeq)))
grouping_factor <- sample_data_df$Kits_type
grouping_factor
distance_matrix <- as.matrix(phyloseq::distance(physeq, method = "bray"))

pairwise_p_values <- pairwise_t_test(distance_matrix, grouping_factor)
#pairwise_p_values_abund <- pairwise_p_values 


write.table(pairwise_p_values, 
            paste0("work/results_21_10_22/water_bray_curtis_t_test.tsv"), 
            sep = "\t", row.names = TRUE, col.names = TRUE)


# database
pval_df = read.csv2('work/results_21_10_22/soil_bray_curtis_t_test.tsv',
                        header=TRUE, check.names = T, sep='\t')

pval_df <- pairwise_p_values
diag(pval_df) <- 1
pval_df


pval_df
pval_df[is.na(pval_df)] <- 0
pval_df

columns <- colnames(pval_df)
for (i in columns){
  print(i)
  pval_df[,as.character(i)] = as.numeric(pval_df[,as.character(i)])
}


pval_df

pval_df <- as.data.frame(sapply(pval_df, as.numeric))
pval_df

pval_df <- pval_df + t(pval_df)
pval_df

pval_df = replace(pval_df, pval_df == 0, 1)
pval_df

mat_rounded <- format(round(pval_df, 10), scientific = TRUE, digits = 10)
mat_rounded

df_long <- melt(pval_df, na.rm = TRUE)


Round_val = c()
for (i in df_long$value){
  if (i > 0.05){
    new_val = round(i, 3)
    Round_val = c(Round_val, new_val)
  } else if (i > 0.001 & i < 0.05) {
    new_val = round(i, 3)
    Round_val = c(Round_val, new_val)
  } else {
   new_val = "***"
   Round_val = c(Round_val, new_val)
  }
}
Round_val

df_long["New_val"] = Round_val


ggheatmap <-ggplot(df_long, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#d7191c", high = "#2c7bb6", mid = "#ffffbf", 
                       midpoint = 0.05, limit = c(0,1), space = "Lab", 
                       name="P-value") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.text.y = element_text(size = 10)) + coord_fixed()

ggheatmap + geom_text(aes(Var2, Var1, label = New_val), color = "black", size = 3) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


ggheatmap <-ggplot(df_long, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colors = custom_shades,
                       values = scales::rescale(c(0, 0.05, 0.5, 0.75, 1)), space = "Lab", 
                       name="P-value") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.text.y = element_text(size = 10)) + coord_fixed()


custom_shades <- c("#FAAB78", "#EEE9DA", "#BDCDD6", "#93BFCF", "#6096B4")

# Red shades
red_shades <- c("#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", "#de2d26")

# Blue shades
blue_shades <- c("#deebf7", "#b3cde3", "#8c96c6", "#8856a7", "#810f7c")

# Green shades
green_shades <- c("#edf8e9", "#bae4b3", "#74c476", "#31a354", "#006d2c")


?scale_fill_gradientn


newOrder = c("Stool", "Microbiome", "Fecal", "B&T", "PwrSoil", 
             "PureLink", "Monarch", "Soil")




# Load required libraries
library(ggplot2)
library(reshape2)

pval_df <- data.frame(Kit_1 = rownames(pval_df), pval_df, row.names = NULL)
pval_df
colnames(pval_df) <- c("Kit_1", newOrder)


df_long <- melt(pval_df, id.vars = "Kit_1", variable.name = "Kit_2", value.name = "P_value")
df_long$Kit_2 = factor(df_long$Kit_2, levels=newOrder)
type(df_long$Kit_2)


ggheatmap <- ggplot(df_long, aes(Kit_1, Kit_2, fill = P_value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#d7191c", high = "#2c7bb6", mid = "#ffffbf",
                       midpoint = 0, limit = c(0, 1), space = "Lab",
                       name = "Pearson\nCorrelation") +
  scale_x_discrete(limits = newOrder) + # Set the x-axis order based on the factor levels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1)) +
  coord_fixed()
ggheatmap$data$Kit_2 <- factor(ggheatmap$data$Kit_2, levels=newOrder)
ggheatmap


ggheatmap <- ggplot(df_long, aes(Kit_1, Kit_2, fill = P_value)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("red", "orange", "yellow", "green", "blue"),
                       values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1))) +
  scale_x_discrete(limits = newOrder) + # Set the x-axis order based on the factor levels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1)) +
  coord_fixed()
ggheatmap





ggplot(df_long, aes(Kit_1, Kit_2, fill = P_value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1)) +
  coord_fixed()


# Plot the heatmap
ggplot(df_long, aes(x = Kit_1, y = Kit_2, fill = P-value)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "blue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(fill = "P-value")



pval_colors <- colorRampPalette(c("blue", "white", "red"))(100)
pval_colors
map_pval_color <- function(p) {
  pval_colors[findInterval(p, seq(0, 1, length.out = 100))]
}

library(vegan)

carbom
set.seed(1)

current_samples <- animal_s

carbom <- phyloseq(OTU_f, TAX, current_samples)
carbom
sample_variables(carbom)
sample_sums(carbom)
carbom



distance_matrix <- phyloseq::distance(carbom, method = "bray")
sample_data_df <- as.data.frame(sample_data(carbom))
#grouping_factor <- sample_data_df$Kits_type
#grouping_factor
?adonis2

my_data <- data.frame(Sample_ID = sample_data_df$Sample_ID,
                      Kits_type = as.factor(sample_data_df$Kits_type))

permanova_result <- adonis2(distance_matrix ~ Kits_type, my_data)
permanova_result

p_values = c()
Df_values = c()
R2_values = c()
F_values = c()
SumOfSqs_values = c()

p_values = c(p_values, permanova_result$"Pr(>F)"[1])


Df_values = c(Df_values, permanova_result$Df[1])

R2_values = c(R2_values, permanova_result$R2[1])
SumOfSqs_values = c(SumOfSqs_values, permanova_result$SumOfSqs[1])

F_values = c(F_values, permanova_result$"F"[1])

p_values
SumOfSqs_values
Df_values
R2_values
F_values

permanova_df <- data.frame(Sample_type = c("Soil", "Water", "Gut"))
permanova_df$Df = Df_values
permanova_df$"SumOfSqs" = SumOfSqs_values
permanova_df$R2 = R2_values
permanova_df$"F" = F_values
permanova_df$"p_value" = p_values


write.table(permanova_df, 
            paste0("work/results_21_10_22/PERMANOVA_stat.tsv"), 
            sep = "\t", row.names = FALSE, col.names = TRUE)















