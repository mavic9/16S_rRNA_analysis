#by Viktor Mamontov, 2022
otu_abundance <- as.matrix(as.data.frame(otu_table(carbom)))
otu_median <- as.matrix(apply(otu_abundance, 1, median))

otu_median = cbind(otu_median, otu_median, otu_median)
colnames(otu_median) <- c("7", "8", "9")

otu_abundance <- cbind(otu_abundance, otu_median)
otu_abundance <- otu_table(otu_abundance, taxa_are_rows = T)


# database
samples_all = read.csv2('work/results_21_10_22/Metadata_test.tsv',
                        header=TRUE, check.names = F, sep='\t', row.names = 1)

samples = sample_data(samples_all)

soil_s = subset(samples, (samples$Expedition.station.number == "76" & samples$Expedition.type == "Marine"
                          | samples$Sample.type == "for_soil"))
soil_s$Source.sample.ID = as.character(soil_s$Source.sample.ID)


### subset for water
water_s = subset(samples, (samples$Expedition.station.number == "Technopark" & samples$Sample_type == "sample" 
                           | samples$Sample.type == "for_water"))
water_s$Source.sample.ID = as.character(water_s$Source.sample.ID)


### subset for animals
animal_s = subset(samples, (samples$Expedition.station.number == "M2" & samples$Sample_type == "sample" 
                            | samples$Sample.type == "for_animal"))
animal_s$Source.sample.ID = as.character(animal_s$Source.sample.ID)




current_samples <- animal_s

carbom <- phyloseq(otu_abundance, TAX, current_samples)
carbom
sample_variables(carbom)
sample_sums(carbom)
carbom

