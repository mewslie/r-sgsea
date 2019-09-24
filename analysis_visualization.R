# Code for PCA visualization of gene expression time series of a group of genes

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

# Identify the cluster to center the test.set1 PCA on ---------------------
# Calculate the distance between time points per gene
# Cluster genes by time distance
# Assign cluster designation to small tight clusters 

make_clusters_by_time <- function(data.norm, type="bygroup") {
  # where data.norm is a scaled and centered expression data
  # where type describes the type of format of the expression data; assumes two groups only!
  # type = "bygroup" means the columns are arranged like group1(T0,1,2,3),group2(T0,1,2,3) etc
  # type = "bycontrol" means the columns are control(rep1,rep2),T1(rep1,rep2),T2(rep1,rep2) etc

  subdiag <- function(dist_obj, d) { #return the dth diagonal of a distance matrix
    n <- attr(dist_obj,"Size")
    j_1 <- c(0, seq.int(from = n - 1, by = -1, length = n - d - 1))
    subdiag_ind <- d + cumsum(j_1)
    dist_obj[subdiag_ind]
  }
  
  if (type=="bygroup") {
    #from:
    # ctrl0 ctrl1 ctrl2 treat0 treat1 treat2
    #to:
    # ctrl0 treat0
    # ctrl1 treat1
    # ctrl2 treat2
    data.dists <- apply(data.norm, 1, function(x) {
      points.dat <- matrix(x, ncol=2, byrow = F) #each point = row
      points.dist <- dist(points.dat) #eucl distance
      subdiag(points.dist,2) #get 2nd diagonal for adjacent time point distance
    })
    data.dists <- t(data.dists)
  } else if (type=="bycontrol") {
    #from:
    # ctrl r1t1 r1t2 r1t3 r2t1 r2t2 r2t3
    #to:
    # t1 t1
    # t2 t2
    # t3 t3
    rep.list <- data.frame(rep1 = seq(3,ncol(data.norm), by=2), rep2 = seq(4,ncol(data.norm), by=2)) #cols of replicates
    data.dists <- apply(data.norm, 1, function(x) {
      rm.nfe.data <- matrix(x, ncol = 2, byrow = T) #row = same time, col = replicate, without NFE data
      points.dist <- dist(rm.nfe.data) #eucl distance
      subdiag(points.dist,2)
    })
    data.dists <- t(data.dists)
  } else { #bytime
    print("Unknown data format!")
    return(NA)
  }
  return(hclust(dist(data.dists), method = "ward.D2")) #Ward
}

get_cluster_assignments <- function(hclust.obj) {
  # where hclust.obj is an object returned from hclust call
  hclust.merging <- hclust.obj$merge
  hclust.merging <- sign(hclust.merging) #neg = leaf, pos = cluster of leaves
  hclust.merging <- cbind(hclust.merging, as.numeric(apply(hclust.merging, 1, function(x) sum(x==-1)))) #add no. of neg no. per row
  #go through each row to flag each merge that only merges leaves
  clust.designation <- rep(0, nrow(hclust.merging)) #"unclustered" merge is 0
  cluster.name <- 0 #keep track of cluster number
  for (i in 1:nrow(hclust.merging)) {
    if (hclust.merging[i,3] == 2 & (i+1 < nrow(hclust.merging))) {
      #hclust.merging[i,3] == 2; start of new cluster designation
      #i+1 < nrow(hclust.merging); can look at next row
      if (hclust.merging[i+1,3] %in% c(1,2)) {
        #hclust.merging[i+1,3] == 1; the next row is a leaf merge
        cluster.name <- cluster.name + 1
        clust.designation[i] <- cluster.name
      }
    } else if (i-1 > 0) {
      #(i-1 > 0); there was a previous row
      if (hclust.merging[i,3] == 1 & hclust.merging[i-1,3] %in% c(1,2) & clust.designation[i-1] != 0) {
        clust.designation[i] <- cluster.name
      }
    }
  }
  #reformat clust.designation for output; get names and cluster designation
  cluster.df <- cbind(hclust.obj$merge, clust.designation)
  cluster.names <- as.data.frame(subset(cluster.df, clust.designation > 0)) #only genes in selected clusters
  cluster.names$one <- NA
  cluster.names$two <- NA
  cluster.names$one[which(cluster.names[,1]<0)] <- hclust.obj$labels[abs(cluster.names[cluster.names[,1]<0,1])]
  cluster.names$two[which(cluster.names[,2]<0)] <- hclust.obj$labels[abs(cluster.names[cluster.names[,2]<0,2])]
  #return vector with names and cluster designation starting from 1; last no. is "unclustered"
  hclust.centers <- as.numeric(c(cluster.names$clust.designation[!is.na(cluster.names$one)],cluster.names$clust.designation[!is.na(cluster.names$two)]))
  names(hclust.centers) <- as.character(c(cluster.names$one[!is.na(cluster.names$one)],cluster.names$two[!is.na(cluster.names$two)]))
  hclust.unclustered <- hclust.obj$labels[!(hclust.obj$labels %in% names(hclust.centers))]
  hclust.centers <- c(hclust.centers, rep(cluster.name+1, length(hclust.unclustered)))
  names(hclust.centers)[nchar(names(hclust.centers))==0] <- hclust.unclustered
  hclust.centers <- hclust.centers[order(match(names(hclust.centers), hclust.obj$labels))]
  return(hclust.centers)
}

# to use you need:
library("data.table")
# scaled and centered data
expr.scale <- as.data.frame(scale(expr.val, center = T, scale = T))
# data of only test.set1
test.set1.expr <- subset(expr.scale, row.names(expr.scale) %in% test.set1)

# use like so:
set1.hclust <- make_clusters_by_time(test.set1.expr, "bygroup") #calc distance and hclust by ward
set1.clusters <- get_cluster_assignments(set1.hclust)


# Identify PCA projection to use ------------------------------------------
# Using previously identified clusters
# Select clusters that are tight using a threshold
# Use PCA on selected clusters only and calculate PCA that gives biggest distance between groups

get_centers <- function(data.hclust, data.clusts, qthreshold) {
  # where data.hclust is the hclust object of make_clusters_by_time()
  # where data.clusts is the cluster assignment from get_cluster_assignments()
  # where qthreshold is the threshold to define a tight cluster
  center.threshold <- quantile(data.hclust$height, qthreshold) #cluster merges above this is "too far"
  data.with.height <- data.frame(clust = data.clusts)
  data.with.height$ind <- match(data.hclust$labels, rownames(data.with.height)) #get clustering labels in same order as df rows
  max.merge.height <- do.call("rbind", as.list(by(data.with.height, list(clust=data.with.height$clust), function(x) {
    hclust.heights <- cbind(as.data.frame(data.hclust$merge), data.hclust$height)
    this.height <- subset(hclust.heights, V1 %in% (x$ind*-1) | V2 %in% (x$ind*-1))
    return(this.height[nrow(this.height),3]) #return highest height for each cluster
  })))
  data.with.height$height <- max.merge.height[match(data.with.height$clust, rownames(max.merge.height)),]
  return(rownames(data.with.height)[data.with.height$height < center.threshold]) #return gene names that are in "central" clusters
}

get_PCA <- function(central.data, pca.runs) {
  # where central.data is the expression of the data centers only which the PCA will be calculated on
  # where pca.runs is a data frame with columns that specify which central.data column to run PCA on and the last row is the column with cluster identity
  all.pcas <- lapply(pca.runs, function(x) {
    PCA(central.data[,x], graph = F, quali.sup = 3:length(x))
  })
  all.pcas.dist <- lapply(all.pcas, function(x) {
    if (length(unique(central.data$group))==1) {
      #there's only 1 group so there's no distance between groups so use the individals (row of central.data)
      return(mean(dist(x$ind$coord)))
    } else {
      return(mean(dist(x$quali.sup$coord)))
    }
  })
  best.pca <- which(unlist(all.pcas.dist) == max(unlist(all.pcas.dist)))
  return(all.pcas[[best.pca]])
}

# to use you need:
library("FactoMineR")
# a quantile threshold for distinguishing central clusters of cluster merge heights
center.quantile <- 0.4
# expression data of the test set with the cluster assignments per gene
test.set1.expr.df <- test.set1.expr
test.set1.expr$group <- factor(set1.clusters)
# dataframe defines the columns in test.set1.expr.df of the same time, and the last column is the one with cluster identity
test.set1.col.def <- data.frame(T0=c(1,5,9),
                                T24=c(2,6,9),
                                T48=c(3,7,9),
                                T72=c(4,8,9))

# use like so:
test.set1.centers <- get_centers(set1.hclust, set1.clusters, center.quantile)
test.set1.centers.df <- test.set1.expr.df[which(rownames(test.set1.expr.df) %in% test.set1.centers),]
test.set1.pca <- get_PCA(test.set1.centers.df, test.set1.col.def)


# Use gganimate and plot test.set1 ----------------------------------------

get_animate_df <- function(this.expr, the.pca, time.txt, time.cols) {
  # where this.expr is the full data to be put into animation df
  # where the.pca is the pca used for prediction
  # where time.txt is the times used for the time column
  # where time.cols is the columns in this.expr for each time
  #initialize the gganimate dataframe for rbind later on
  data.animate <- data.frame(gene=character(), Dim.1=numeric(), Dim.2=numeric(), clust=numeric(), time=numeric(), stringsAsFactors = F)
  #calculate remaining predictions
  make_prediction <- function(this.data, predict.data, header, time, factors) {
    #this.data is the data to be used in prediction, predict.data is the data used to make prediction, header is new colnames, time is a number used to make time column with T#, factors is additional factor variable to include for this.data
    colnames(this.data) <- rownames(predict.data$var$coord) #same var name to this.data as predict.data
    this.prediction <- predict(predict.data, this.data) #do prediction
    this.prediction.df <- as.data.frame(this.prediction$coord) #take new coordinates
    this.prediction.df <- cbind(gene=rownames(this.prediction.df), this.prediction.df) #add gene names
    rownames(this.prediction.df) <- c() #remove rownames
    this.prediction.df <- cbind(this.prediction.df, factors[which(rownames(factors) %in% this.prediction.df$gene),,drop=F]) #add factors
    rownames(this.prediction.df) <- c()
    this.prediction.df$time <- rep(time, nrow(this.prediction.df)) #add time
    colnames(this.prediction.df) <- header
    this.prediction.df
  }
  #calculate remaining predictions in time.cols
  predict.list <- lapply(1:ncol(time.cols), function(predict.loop) {
    isolate.factors <- this.expr[,tail(time.cols[,predict.loop,drop=T],-2),drop=F] #factors are this.expr but without the initial dimension data
    predict <- make_prediction(this.expr[,time.cols[,predict.loop,drop=T]], the.pca, colnames(data.animate), time.txt[predict.loop], isolate.factors)
  })
  predict.list.to.df <- do.call("rbind", predict.list)
  data.animate <- predict.list.to.df
  colnames(data.animate) <- c("gene","Dim.1","Dim.2","clust","time")
  return(data.animate)
}

# to use you need:
library("ggplot2")
library("gganimate")
library("gifski")
library("png")
# the times of the time column
time.text <- c(0,24,48,72)
time.column <- data.frame(T0=c(1,5,9),
                          T24=c(2,6,9),
                          T48=c(3,7,9),
                          T72=c(4,8,9))

# use like so:
test.set1.animate <- get_animate_df(test.set1.expr.df, test.set1.pca, time.text, time.column)
ggplot(test.set1.animate, aes(x = Dim.1, y=Dim.2, colour = clust)) +
  geom_point(show.legend = FALSE, alpha = 0.7) +
  labs(x = "PC1", y = "PC2") +
  transition_time(time)

