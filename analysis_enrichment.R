# Code for serialized gene set enrichment analysis

# Load expression data ----------------------------------------------------
expr.fc
# assuming that the fold change data has been calculated between two groups
# assuming that the fold change data is in the form like this:
# gene_id     | Time1 | Time2 | Time3  ...
# g1          | 
# g2          |
# g3          |
# g4          |
# g5          |

# Load test set/s ---------------------------------------------------------
# assuming that there are groups of genes you want to test
# assuming the gene ids are the same as expression data, load them into a (R) character vector

test.set1
test.set2


# Calculate the test value for all genes ------------------------------------
# uses the absolute adjacent time point difference in fold change per gene
# the calculation labels time intervals by alphabetical letters

get_adj_diff <- function(geneexpr, timeunit) {
  #where geneexpr is fold change in chronological order by col
  #where timeunit is the no. of time units between samples
  adjacent.differences <- geneexpr/timeunit
  names(adjacent.differences) <- LETTERS[1:length(adjacent.differences)]
  return(abs(adjacent.differences))
}

# to use you need:
# time point differences are standardized into the unit of time used in the data set, eg hours, days
time.units <- c(24,24,24) #eg there are 24 hours between each sample
#time.units <- c(1,1,4,7) #eg there are 1, 1, 4, 7 hours between each corresponding sample

# use like so:
all.adj <- as.data.frame(t(as.data.frame(apply(expr.fc[,-1], 1, function(x) get_adj_diff(x, time.units)))))
row.names(all.adj) <- expr.fc[,1]


# Apply test on test.set1 ----------------------------------------------
# Calculate the probability of seeing at least q success from size trials at p rate
# Success = test value is higher than the Q% quantile of center of data
# Center of data = densist cluster of genes with high average test value
# q = no. of genes with test value a Q% quantile distance away from center of data
# size = total number of genes in test.set1
# p = probability of a gene with test value further than the Q% quantile of center of data, from full data set

get_quantile_threshold_interval <- function(values, quant) {
  # where values is the test value values in a matrix where time intervals are columns and rows are genes
  # where quant is the quantile of the values you want to calculate, in decimal like 0.5 for median
  values.dist <- as.matrix(dist(values)) #distance matrix between genes
  average.dist <- colMeans(values.dist) #average distance per gene over all times
  min.dist <- average.dist[average.dist == min(average.dist)] #get minimum average distance
  target.gene <- values[rownames(values) %in% names(min.dist),] #gene data with the minimum average distance
  target.gene <- target.gene[order(target.gene$A, target.gene$B, target.gene$C, target.gene$D, decreasing = TRUE),] #sort largest to smallest from left to right col
  target.gene <- target.gene[1,] #gene with max value with min average distance; gene that is closest to all other genes
  target.threshold <- target.gene
  return(target.threshold)
}                                            

get_successes_interval <- function(value, data) { #number data with value or higher
  # where value is the threshold of which to count in data
  # where data is the SE values in a matrix where time intervals are columns and rows are genes
  #return(sum(rowSums(t(apply(data, 1, function(x) x >= value)))==length(value))) #more conservative
  return(sum(rowSums(t(apply(data, 1, function(x) x >= value)))>0)) #more liberal
}

get_probability_interval <- function(value, data) { #ratio of data with value or higher
  # where value is the threshold of which to count in data
  # where data is the SE values in a matrix where time intervals are columns and rows are genes
  #subset_data_by_value <- sum(rowSums(t(apply(data, 1, function(x) x >= value)))==length(value))
  #at least one column of data has to be higher than corresponding value; more liberal
  subset_data_by_value <- sum(rowSums(t(apply(data, 1, function(x) x >= value)))>0)
  return(subset_data_by_value/nrow(data))
}
                                            
# to use you need:
# decide on a quantile distance for calculating the test threshold 
quantile.dist <- 0.25
# the test threshold
test.set1.sethreshold <- get_quantile_threshold_interval(all.adj[which(rownames(all.adj) %in% test.set1),], quantile.dist)

# ecdf of SE from full data set
se.ecdf <- ecdf(all.adj)
# the parameters for the test
se.q <- get_successes_interval(test.set1.sethreshold, all.adj[which(rownames(all.adj) %in% test.set1),])-1
se.size <- length(test.set1)
se.p <- get_probability_interval(binom3.sim.strongdeg, all.adj)

# use like so:
pbinom(q = se.q, size = se.size, prob = se.p, lower.tail = FALSE)

