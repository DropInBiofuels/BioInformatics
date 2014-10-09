## This script creates a distance matrix from an multiple sequence alignment and finds the 
## Proteins with the greatest Distance to each other. This works quite fast (approx. 10 min on an i5)
## for very large alignments (>10K Sequences), but requires some memory (approx. 5GB) for each calculation.
## Thus please run the calculation of the two matrices seperatly.
## J. Kabisch 2013
## kabisch@uni-greifswald.de

# load library and sequence alignment
library(seqinr)
setwd("/home/jfk/Arbeit/Projekte/Alberto/")
myseqs <- read.alignment("Full3dm_orn.fasta", format = "fasta")

# calculate distance based on similarity and output max distance and protein pair
mat2simi <- as.matrix(dist.alignment(myseqs, matrix = "similarity" )) # calculate distances
maxPosition <- arrayInd(which.max(mat2simi), dim(mat2simi)) # find max. distance in matrix
MaxSimiProteinName<- rownames(mat2simi)[maxPosition[1,]] # extract name of protein pair with max. distance
MaxSimiValue <- mat2simi[maxPosition] # extract distance value

# calculate distance based on similarity and output max distance and protein pair
mat2ident <- as.matrix(dist.alignment(myseqs, matrix = "identity" )) # calculate distances
maxPosition <- arrayInd(which.max(mat2ident), dim(mat2ident)) # find max. distance in matrix
MaxIdentProteinPair<- rownames(mat2ident)[maxPosition[1,]] # extract name of protein pair with max. distance
MaxIdentValue <- mat2ident[maxPosition] # extract distance value

