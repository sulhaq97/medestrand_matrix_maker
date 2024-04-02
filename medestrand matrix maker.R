setwd("C:/Users/Sami/OneDrive - University of Toronto/Masters/Code/Reference Scripts/Basics_Bioinformatics/")

library(MeDEStrand)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Repitools)

BSgenome <- "BSgenome.Hsapiens.UCSC.hg19"

# Only one read is kept if multiple reads are mapped to the same position
uniq = 1

# Extend the reads to the estimated fragment length
extend = 200
shift = 0

# window size
ws = 300 

chr.select = paste0('chr', c(22:22) )

bam.files.names <- Sys.glob("*bam")

# creates 300bp tiled genome
genome_300bp <- data.frame(genomeBlocks(BSgenome.Hsapiens.UCSC.hg19, chrs=seqnames(BSgenome.Hsapiens.UCSC.hg19)[1:22], width=300))
genome_300bp$name_format <- paste(genome_300bp$seqnames, genome_300bp$start, genome_300bp$end, sep=".")

# this matrix will hold all medestrand beta values
matrix.of.pbls <- matrix(nrow=length(genome_300bp$name_format))
# the rownames correspond to windows
rownames(matrix.of.pbls) <- genome_300bp$name_format

# each bam file is examined and the beta-value is infered using MeDEStrand and then added as a column to the holder matrix
for(each.bam.file in bam.files.names) {
  MeDIP_seq = MeDEStrand.createSet(file=each.bam.file, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, window_size=ws, chr.select=chr.select, paired = T)
  
  #  count CpG pattern in the bins
  CS = MeDEStrand.countCG(pattern="CG", refObj=MeDIP_seq)
  result.methylation = MeDEStrand.binMethyl(MSetInput = MeDIP_seq, CSet = CS, Granges = FALSE)
  
  # each sample is added
  matrix.of.pbls <- cbind(matrix.of.pbls, result.methylation)
}

# the first column of the matrix contains NAs
matrix.of.pbls <- matrix.of.pbls[,-c(1)]

# the colnames are the same as bam file names 
colnames(matrix.of.pbls) <- substr(bam.files.names, 1, 9)
