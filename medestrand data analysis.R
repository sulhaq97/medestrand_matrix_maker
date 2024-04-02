############################################
# This script finds hypermethylated windows
# in PBLs.
#
# Written by Sami Ul Haq
# Sept 15, 2021
############################################

setwd("C:/Users/Sami/OneDrive - University of Toronto/Masters/Code/Reference Scripts/MeDEStrand Stuff")

# loads the matrix containing MeDEStrand converted beta-values
load("medestrand.300bp.matrix.of.pbls.RData")

# this examines the median beta values per windwo
median.beta.vals <- apply(matrix.of.pbls, MARGIN=1, FUN = median)
names(median.beta.vals) <- rownames(matrix.of.pbls)

# Filter out blacklist windows
load("hg19.encode.blacklist.windows.RData")
# matrix is removed for blacklist windows
median.beta.vals <- median.beta.vals[ !(names(median.beta.vals) %in% hg19.encode.blacklist.windows) ]


# selects windows with beta values greater than 0.7
hyper.meth.windows <- median.beta.vals[which(median.beta.vals > 0.7)]

save(hyper.meth.windows, file="hypermethylated.windows.PBLs.RData")

