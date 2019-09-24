# Code for enrichment analysis for average variation enrichment (AE) and sequential variation enrichment (SE)

# Load expression data ----------------------------------------------------
expr.val
# assuming that the expression data has been normalized by DESeq2, edgeR, etc
# assuming that the expression data has the form like this:
# gene_id     | Group1.Time1 | Group2.Time1 | Group1.Time2 | Group2.Time2 ...
# g1          | 
# g2          |
# g3          |
# g4          |
# g5          |
# or like this:
# assuming that the expression data has the form like this:
# gene_id     | Group1.Time1 | Group1.Time2 | Group1.Time3 | Group2.Time1 ...
# g1          | 
# g2          |
# g3          |
# g4          |
# g5          |
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

# Calculate the AE value for all genes ------------------------------------
# AE uses the absolute mean adjacent time point difference in fold change per gene

get_avg_adj_diff <- function(geneexpr, timeunit) {
  #where geneexpr is fold change in chronological order by col
  #where timeunit is the no. of time units between samples
  adjacent.differences <- tail(geneexpr,-1) - head(geneexpr,-1)
  adjacent.differences <- adjacent.differences/timeunit
  return(mean(abs(adjacent.differences)))
}

# to use you need:
# time point differences are standardized into the unit of time used in the data set, eg hours, days
time.units <- c(24,24,24) #eg there are 24 hours between each sample
#time.units <- c(1,1,4,7) #eg there are 1, 1, 4, 7 hours between each corresponding sample

# use like so:
average.adj <- apply(expr.fc[,-1], 1, function(x) get_avg_adj_diff(x, time.units))
names(average.adj) <- expr.fc[,1]


# Calculate the SE value for all genes ------------------------------------
# SE uses the absoluate adjacent time point difference in fold change per gene
# the calculation labels time intervals by alphabetical letters

get_adj_diff <- function(geneexpr, timeunit) {
  #where geneexpr is fold change in chronological order by col
  #where timeunit is the no. of time units between samples
  adjacent.differences <- tail(geneexpr,-1) - head(geneexpr,-1)
  adjacent.differences <- adjacent.differences/timeunit
  names(adjacent.differences) <- LETTERS[1:length(adjacent.differences)]
  return(abs(adjacent.differences))
}

# to use you need:
# time point differences are standardized into the unit of time used in the data set, eg hours, days
time.units <- c(24,24,24) #eg there are 24 hours between each sample
#time.units <- c(1,1,4,7) #eg there are 1, 1, 4, 7 hours between each corresponding sample

# use like so:
all.adj <- as.data.frame(t(as.data.frame(apply(expr.fc.fc[,-1], 1, function(x) get_adj_diff(x, time.units)))))
row.names(all.adj) <- expr.fc[,1]


# Apply AE test on test.set1 ---------------------------------------------
# Calculate the probability of seeing at least q success from size trials at p rate
# Success = AE is higher than the Q% quantile of center of data
# Center of data = densist cluster of genes with highest AE
# q = no. of genes with AE higher than the Q% quantile of center of data
# size = total number of genes in test.set1
# p = probability of a gene with AE higher than the Q% quantile of center of data, from full data set

get_quantile_threshold <- function(values, quantile) {
  # where values are AE values
  # where quantile is the quantile of the values you want to calculate, in decimal like 0.5 for median
  values.dist <- as.matrix(dist(values))
  average.dist <- colMeans(values.dist)
  min.dist <- average.dist[average.dist == min(average.dist)] #find center of data by looking for gene with minimum average distance to all other genes
  target.gene <- values[names(values) %in% names(min.dist)]
  target.gene <- target.gene[target.gene == max(target.gene)] #max stat value with min average distance
  target.threshold <- target.gene - quantile(values.dist[,colnames(values.dist) %in% names(target.gene)], quantile)
  return(target.threshold)
}

# to use you need:
# decide on a quantile distance for calculating the test threshold 
quantile.dist <- 0.25
# the test threshold
test.set1.aethreshold <- get_quantile_threshold(average.adj[which(names(average.adj) %in% test.set1)], quantile.dist)
# ecdf of AE from full data set
ae.ecdf <- ecdf(average.adj)
# the parameters for the test
ae.q <- sum(average.adj[which(names(average.adj) %in% test.set1)] >= test.set1.aethreshold)-1
ae.size <- length(test.set1)
ae.p <- 1-ae.ecdf(test.set1.aethreshold)

# use like so:
pbinom(q = ae.q, size = ae.size, prob = ae.p, lower.tail = FALSE)

# Apply SE test on test.set1 ----------------------------------------------
# Calculate the probability of seeing at least q success from size trials at p rate
# Success = SE is higher than the Q% quantile of center of data
# Center of data = densist cluster of genes with high average SE
# q = no. of genes with SE a Q% quantile distance away from center of data
# size = total number of genes in test.set1
# p = probability of a gene with SE further than the Q% quantile of center of data, from full data set

get_quantile_threshold_interval <- function(values, quant) {
  # where values is the SE values in a matrix where time intervals are columns and rows are genes
  # where quant is the quantile of the values you want to calculate, in decimal like 0.5 for median
  values.dist <- as.matrix(dist(values)) #distance matrix between genes
  average.dist <- colMeans(values.dist) #average distance per gene over all times
  min.dist <- average.dist[average.dist == min(average.dist)] #get minimum average distance
  target.gene <- values[rownames(values) %in% names(min.dist),] #gene data with the minimum average distance
  target.gene <- target.gene[do.call(order, as.list(target.gene), decreasing = TRUE),] #sort largest to smallest from left to right col
  target.gene <- target.gene[1,] #gene with max value with min average distance; gene that is closest to all other genes
  interval.matrix <- do.call(rbind, apply(values, 1, function(x) x - target.gene))
  interval.matrix.col <- which(rownames(values) %in% rownames(target.gene)) #need index to remove diagonal
  interval.matrix.fixed <- interval.matrix[-interval.matrix.col,]
  target.threshold.lhs <- do.call(c, lapply(interval.matrix.fixed, function(x) quantile(x, quant)))
  target.threshold <- target.gene + target.threshold.lhs
  return(target.threshold)
}

get_successes_interval <- function(value, data) { #number data with value or higher
  # where value is the threshold of which to count in data
  # where data is the SE values in a matrix where time intervals are columns and rows are genes
  return(sum(rowSums(t(apply(data, 1, function(x) x >= value)))==length(value)))
}

get_probability_interval <- function(value, data) { #ratio of data with value or higher
  # where value is the threshold of which to count in data
  # where data is the SE values in a matrix where time intervals are columns and rows are genes
  subset_data_by_value <- sum(rowSums(t(apply(data, 1, function(x) x >= value)))==length(value))
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

