#by Viktor Mamontov, 2022
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")

contam_table = OTU_name

cur_samples = as.character(soil_all$Sample_ID)

cur_control
cur_samples

for (i in seq(0, 7)){
  n = i + 1
  j = i * 3
  print(n)
  test_s = subset(samples, 
                  samples$Sample_ID == cur_control[n] | 
                  samples$Sample_ID == cur_samples[j+1] |
                  samples$Sample_ID == cur_samples[j+2] | 
                  samples$Sample_ID == cur_samples[j+3])
  conc_val = as.numeric(test_s$DNA.concentration..ng.mkl)
  print(conc_val)
  current_samples <- test_s
  carbom <- phyloseq(OTU, TAX, current_samples)
  print(sample_sums(carbom))
  sample_data(carbom)$is.neg <- sample_data(carbom)$Sample_type == "control";
  contamdf.prev <- isContaminant(carbom, conc = conc_val, method="either", neg="is.neg", 
                                 threshold=0.5, normalize = T)
  print(table(contamdf.prev$contaminant))
  contam_table[cur_samples[j+1]] = contamdf.prev$contaminant
  contam_table[cur_samples[j+2]] = contamdf.prev$contaminant
  contam_table[cur_samples[j+3]] = contamdf.prev$contaminant
}

cur_samples = as.character(water_s$Sample_ID)
cur_control = as.character(water_control$Sample_ID)

for (i in seq(0, 7)){
  n = i + 1
  j = i * 3
  print(n)
  test_s = subset(samples, 
                  samples$Sample_ID == cur_control[n] | 
                    samples$Sample_ID == cur_samples[j+1] |
                    samples$Sample_ID == cur_samples[j+2] | 
                    samples$Sample_ID == cur_samples[j+3])
  conc_val = as.numeric(test_s$DNA.concentration..ng.mkl)
  print(conc_val)
  current_samples <- test_s
  carbom <- phyloseq(OTU, TAX, current_samples)
  print(sample_sums(carbom))
  sample_data(carbom)$is.neg <- sample_data(carbom)$Sample_type == "control";
  contamdf.prev <- isContaminant(carbom, conc = conc_val, method="either", neg="is.neg", 
                                 threshold=0.5, normalize = T)
  print(table(contamdf.prev$contaminant))
  contam_table[cur_samples[j+1]] = contamdf.prev$contaminant
  contam_table[cur_samples[j+2]] = contamdf.prev$contaminant
  contam_table[cur_samples[j+3]] = contamdf.prev$contaminant
}

cur_samples = as.character(animal_s$Sample_ID)
cur_control = as.character(animal_control$Sample_ID)

for (i in seq(0, 7)){
  n = i + 1
  j = i * 3
  print(n)
  test_s = subset(samples, 
                  samples$Sample_ID == cur_control[n] | 
                    samples$Sample_ID == cur_samples[j+1] |
                    samples$Sample_ID == cur_samples[j+2] | 
                    samples$Sample_ID == cur_samples[j+3])
  conc_val = as.numeric(test_s$DNA.concentration..ng.mkl)
  print(conc_val)
  current_samples <- test_s
  carbom <- phyloseq(OTU, TAX, current_samples)
  print(sample_sums(carbom))
  sample_data(carbom)$is.neg <- sample_data(carbom)$Sample_type == "control";
  contamdf.prev <- isContaminant(carbom, conc = conc_val, method="either", neg="is.neg", 
                                 threshold=0.5, normalize = T)
  print(table(contamdf.prev$contaminant))
  contam_table[cur_samples[j+1]] = contamdf.prev$contaminant
  contam_table[cur_samples[j+2]] = contamdf.prev$contaminant
  contam_table[cur_samples[j+3]] = contamdf.prev$contaminant
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

save_table(decontam_otu, "decontam_otu_05_conc.tsv")
save_table(contam_otu, "contam_otu_05_conc.tsv")

contam_otu_relative = contam_otu
for (i in seq(1, length(otu_frq))){
  total = sum(otu_frq[i])
  contam_otu_relative[i] = round(contam_otu[i]/total, 4) * 100
  print(i)
}

save_table(contam_otu_relative, "contam_otu_rel_05_conc.tsv")
