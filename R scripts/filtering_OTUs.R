#by Viktor Mamontov, 2022
contam_table = OTU_name
current_samples = animal_all
conc_val = as.numeric(current_samples$DNA.concentration..ng.mkl)
print(conc_val)
carbom <- phyloseq(OTU, TAX, current_samples)
sample_data(carbom)$is.neg <- sample_data(carbom)$Sample_type == "control";
contamdf.prev <- isContaminant(carbom, conc = conc_val, method="either", neg="is.neg", 
                               threshold=0.5, normalize = T)
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

carbom.pos <- prune_samples(sample_data(carbom)$Sample_type == "sample", carbom)
carbom.noncontam <- prune_taxa(!contamdf.prev$contaminant, carbom.pos)
carbom.noncontam

sample_ID = sample_data(carbom.noncontam)$Sample_ID
for (i in sample_ID){
  contam_table[as.character(i)] = contamdf.prev$contaminant
}


write.table(contam_table, "work/results_21_10_22/tables/contam_table_05_conc.tsv", 
            sep = "\t", row.names = TRUE, col.names = TRUE)

rownames(contam_table) = contam_table$OTU
contam_table$OTU = NULL

current_samples <- all_samples
carbom <- phyloseq(OTU, TAX, current_samples)
carbom

?paste

otu_frq <- data.frame(otu_table(carbom), check.names = F)
decontam_otu <- otu_frq
contam_otu <- otu_frq
print(length(otu_frq))
print(length(otu_frq[,1]))

for (i in seq(1, length(otu_frq[,1]))){
  for (j in seq(1, length(otu_frq))){
    if (contam_table[i,j]){
      decontam_otu[i,j] = 0
    }
    else {
      contam_otu[i,j] = 0
    }
  }
  if (i %% 10000 == 0){
    print(paste(i, "/67045", sep=""))
  }
}

save_table <- function(df, name){
  write.table(df, paste("work/results_21_10_22/tables/", name, sep=""),
              sep = "\t", row.names = TRUE, col.names = TRUE)
}

save_table(decontam_otu, "decontam_otu_05_conc_all.tsv")
save_table(contam_otu, "contam_otu_05_conc_all.tsv")

contam_otu_relative = contam_otu
for (i in seq(1, length(otu_frq))){
  total = sum(otu_frq[i])
  contam_otu_relative[i] = round(contam_otu[i]/total, 4) * 100
  print(i)
}

save_table(contam_otu_relative, "contam_otu_rel_05_conc_all.tsv")


