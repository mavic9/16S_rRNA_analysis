#by Viktor Mamontov, 2022
library(phyloseq)
library(ape)
library(vegan)
library(ade4)

WT = function(dm, f){
  f = as.factor(f)
  if(nlevels(f) != 2) return(NULL)
  lev = levels(f)
  ns = summary(f)
  N = sum(ns)
  dd = as.matrix(dm)
  dd[upper.tri(dd)]=0 ##
  SST = sum(dd^2)/N
  SSW1 = sum(dd[f==lev[1],f==lev[1]]^2)/ns[1] 
  SSW2 = sum(dd[f==lev[2],f==lev[2]]^2)/ns[2]
  SSW = SSW1 + SSW2
  
  s1 = SSW1/(ns[1]-1)
  s2 = SSW2/(ns[2]-1)
  
  t2.stat = (((ns[1]+ns[2])/(ns[1]*ns[2])))*(SST-SSW)/(s1/ns[1] + s2/ns[2])
  t2.stat
}


?sample()


wt_pairwise_test <- function(physeq, grouping_factor, nrep=999){
  unique_groups <- unique(grouping_factor)
  n_groups <- length(unique_groups)
  print(unique_groups)
  print(n_groups)
  mdf = as.data.frame(sample_data(physeq))
  # Initialize a matrix to store pairwise p-values
  p_values <- matrix(NA, nrow = n_groups, ncol = n_groups)
  rownames(p_values) <- colnames(p_values) <- unique_groups
  for (i in 1:(n_groups - 1)) {
    for (j in (i + 1):n_groups) {
      print(paste(unique_groups[i], unique_groups[j]))
      tested_data = subset(mdf, mdf$Kits_type == unique_groups[i] | 
                             mdf$Kits_type == unique_groups[j])
      tested_data = sample_data(tested_data)
      tested_physeq = phyloseq(OTU_f, TAX, tested_data)
      
      dm = phyloseq::distance(tested_physeq, method="bray")
      f = sample_data(tested_physeq)$Kits_type
      stats = c(WT(dm, f), replicate(nrep, WT(dm, f[sample(length(f))])))
      p_values[i, j] = sum(stats>=stats[1])/(nrep+1)
      print(p_values[i, j])
    }
  }
  return(p_values)
}



carbom
sample_data(carbom)$Sample.type[1]
grouping_factor <- sample_data(carbom)$Kits_type
wt_pairwise_table <- wt_pairwise_test(carbom, grouping_factor, nrep=999)

write.table(wt_pairwise_table, 
            paste0("work/results_21_10_22/animal_bray_curtis_welch_test.tsv"), 
            sep = "\t", row.names = TRUE, col.names = TRUE)




