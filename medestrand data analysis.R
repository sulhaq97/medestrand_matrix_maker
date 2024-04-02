#######################################################
# This script converts MeDIP data to MeDEStrand data 
# and finds windows of interest with specific beta-val
# cut-offs
#
# Written by Sami Ul Haq
# Sept 15, 2021
############################################

# loads the matrix containing MeDEStrand converted beta-values
load("saved_MeDEStrand_matrix_of_interest.RData")

# this examines the median beta values per windwo
median.beta.vals <- apply(matrix.of.medestrand, MARGIN=1, FUN = median)
names(median.beta.vals) <- rownames(matrix.of.medestrand)

# Filter out blacklist windows
load("hg19.encode.blacklist.windows.RData")
# matrix is removed for blacklist windows
median.beta.vals <- median.beta.vals[ !(names(median.beta.vals) %in% hg19.encode.blacklist.windows) ]


# selects windows with beta values greater than 0.7
hyper.meth.windows <- median.beta.vals[which(median.beta.vals > 0.7)]

save(hyper.meth.windows, file="hypermethylated.windows.RData")

