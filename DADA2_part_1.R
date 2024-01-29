#by Viktor Mamontov, 2022
library(dada2); packageVersion("dada2")

args <- commandArgs(trailingOnly = TRUE)

print('Start 16S_script.R')
directory <- args[1] ### The first input after start pipeline
path <- paste0(directory, '/trimmed')
list.files(path)
fnFs <- sort(list.files(path, pattern="fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), ".fastq"), `[`, 1)
filtFs <- file.path(directory, "filtered", paste0(sample.names, "_filt.fastq.gz"))
names(filtFs) <- sample.names
print('Do output of samples')
out <- filterAndTrim(fnFs, filtFs, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
rdataFile <- "./metagenomes.RData"
save.image(rdataFile)
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
save.image(rdataFile)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaFs[[1]]
save.image(rdataFile)
seqtab <- makeSequenceTable(dadaFs)

#Filtering seqtab table

print('######################')
print('Reads table dimension:')
dim(seqtab)
thresh <- args[2] ### The second input after start pipeline
seqtab = seqtab[,colSums(seqtab)>thresh]
print('')
print('######################')
print(paste0('Reads table dimension after filtration: ', thresh))
dim(seqtab)

table(nchar(getSequences(seqtab)))
save.image(rdataFile)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
save.image(rdataFile)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
save.image(rdataFile)

print('Save information:frequency ASV from seqtab.nochim')
#save information:frequency ASV from seqtab.nochim
seqtab.nochim.export <- t(seqtab.nochim)
rownames(seqtab.nochim.export) <- paste("ASV_", seq(nrow(seqtab.nochim.export)), sep = '')
head(seqtab.nochim.export)
write.table(seqtab.nochim.export, file = paste(directory, "/seqtabnochim_n.tsv", sep = ""), col.names = NA, sep = '\t')

#save information: make fasta file with ASV sequences
library("ape")
fastaMx <- as.matrix(colnames(seqtab.nochim))
rownames(fastaMx) <- rownames(seqtab.nochim.export)
write.dna(fastaMx, file = paste(directory, "/seqtabnochim.fasta", sep = ""), 
          format="fasta", nbcol = -1, colw = 700)
