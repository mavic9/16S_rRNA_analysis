#by Viktor Mamontov, 2022
library("microbiomeMarker")

mm_welch <- run_test_two_groups(
  carbom,
  group = "Kits_type",
  method = "welch.test"
)

mm_welch

mm_welch <- function(otu_table_matrix, grouping_factor) {
  unique_groups <- unique(grouping_factor)
  n_groups <- length(unique_groups)
  print(unique_groups)
  print(n_groups)
  # Initialize a matrix to store pairwise p-values
  p_values <- matrix(NA, nrow = n_groups, ncol = n_groups)
  rownames(p_values) <- colnames(p_values) <- unique_groups
  
  for (i in 1:(n_groups - 1)) {
    for (j in (i + 1):n_groups) {
      factors <- c(unique_groups[i], unique_groups[j])
      set_samples = subset(samples, samples$Kits_type == factors[1] |
                             samples$Kits_type == factors[2])
      set_OTU <- otu_table(otu_table_matrix, taxa_are_rows = TRUE)
      physeq <- phyloseq(set_OTU, TAX, set_samples)
      mm_welch <- run_test_two_groups(
        physeq,
        group = "Kits_type",
        method = "welch.test"
      )
      
      p_values[i, j] <- kruskal_result$p.value
    }
  }
  
  return(p_values)
}

factors <- c(unique_groups[i], unique_groups[j])
factors
set_samples = subset(samples, samples$Kits_type == factors[1] |
                       samples$Kits_type == factors[2])
set_OTU <- otu_table(otu_table_matrix, taxa_are_rows = TRUE)
physeq <- phyloseq(set_OTU, TAX, set_samples)
physeq

mm_welch <- run_test_two_groups(
  physeq,
  group = "Kits_type",
  method = "welch.test",
  norm = "TSS"
)
mm_welch


carbom
metadata <- as.data.frame(sample_data(carbom))
dist_mat <- as.matrix(phyloseq::distance(carbom, method = 'bray'))
clustering_rows <- hclust(as.dist(dist_mat), method = "average")
clustering_cols <- hclust(as.dist(dist_mat), method = "average")

pheatmap(dist_mat,
         clustering_distance_rows = clustering_rows,
         clustering_distance_cols = clustering_cols,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),
         border_color = NA)





stat_res <- anosim(dist_mat, metadata$Kits_type, permutations = 1000)





library(amplicon)

?beta_pcoa_stat()

stat_res <- beta_pcoa_stat(
  dis_mat = dist_mat,
  metadata = metadata,
  groupID = "Kits_type",
  result = "beta_pcoa_stat.txt",
  pairwise = T,
  pairwise_list = "vignettes/compare.txt"
)

design2 = subset(metadata, group %in% c(sampleA, sampleB))

adonis_table = adonis2(sub_dis_table~group, data=design2, permutations = 10000)







pht <- run_posthoc_test(carbom, group = "Kits_type")

mg_anova <- run_test_multiple_groups(
  carbom,
  group = "Kits_type",
  method = "anova"
)
mg_anova

plot_ef_dot(pht)










