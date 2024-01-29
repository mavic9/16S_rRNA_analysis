#by Viktor Mamontov 2022
pairwise_permanova <- function(distance_matrix, grouping_factor, method = "holm") {
  unique_groups <- unique(grouping_factor)
  n_groups <- length(unique_groups)
  print(n_groups)
  
  # Initialize a matrix to store pairwise p-values
  p_values <- matrix(NA, nrow = n_groups, ncol = n_groups)
  rownames(p_values) <- colnames(p_values) <- unique_groups
  
  # Perform pairwise PERMANOVA tests
  for (i in 1:(n_groups - 1)) {
    for (j in (i + 1):n_groups) {
      subset <- grouping_factor %in% c(unique_groups[i], unique_groups[j])
      print(distance_matrix[subset, subset])
      #pairwise_result <- adonis2(distance_matrix[subset, subset] ~ grouping_factor[subset])
      pairwise_result <- anosim(distance_matrix[subset, subset], grouping_factor[subset])
      p_values[i, j] <- pairwise_result$'Pr(>F)'[1]
    }
  }
  
  # Adjust p-values for multiple comparisons
  #p_values <- p.adjust(p_values, method = method)
  
  # Make the matrix symmetric
  #p_values <- p_values + t(p_values)
  
  return(p_values)
}


physeq = carbom
sample_sums(physeq)


dist_matrix <- as.matrix(phyloseq::distance(physeq, method = "bray"))

grouping_factor <- as.factor(sample_data(physeq)$Kits_type)

pairwise_p_values <- pairwise_permanova(dist_matrix, grouping_factor)

pairwise_p_values




library(reshape2)

# Convert the OTU table to a matrix
otu_table_matrix <- as.matrix(otu_table(physeq))

# Calculate the distance matrix (Bray-Curtis dissimilarity)
dist_matrix <- as.matrix(phyloseq::distance(physeq, method = "bray"))

# Convert the distance matrix into a data frame
dist_df <- as.data.frame(dist_matrix)

# Add the grouping factor as a new column
dist_df$Group <- sample_data(physeq)$Kits_type


# Melt the distance data frame
dist_melted <- melt(dist_df, id.vars = "Group")

# Boxplot
ggplot(dist_melted, aes(x = Group, y = value)) +
  geom_boxplot() +
  xlab("Group") +
  ylab("Beta diversity") +
  theme_minimal()

kruskal_test <- kruskal.test(value ~ Group, data = dist_melted)



pathotype.anosim <- anosim(dist_matrix, grouping_factor)
print(pathotype.anosim)




color_scale <- colorRampPalette(c("red", "white", "blue"))(50)

# Create a heatmap using the pheatmap function
pheatmap(pairwise_p_values,
         color = color_scale,
         display_numbers = TRUE,
         number_format = "%.2f",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Pairwise PERMANOVA P-values")


