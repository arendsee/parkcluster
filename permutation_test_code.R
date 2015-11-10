## R Code for permuation based test for determining the significance of
## clusters as described in:

## Park, P.J., Manjourides, J., Bonetti, M., and Pagano, M.
## A permutation test for determining significance of clusters with
## applications to spatial and gene expression data.
## Journal of Computational Statistics and Data Analysis. 53:4290-4300, 2009.

## Copyright 2001,2009
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


## -------------------------------------------------------------------
## This code was refactored by Zebulun Arendsee <arendsee@iastate.edu>
## -------------------------------------------------------------------



#======================================================================
# FUNCTIONS FOR CALCULATE P-VALUES
#======================================================================

#----------------------------------------------------------------------
# Given a tree ($merge) and index (ith splitting from the bottom) and 
# branch (left or right), get all the elements in that branch.  
# $merge contains a negative number if it is at the end; a positive 
# number points to a subsequent branch.  This subroutine updates
# until all the items are negative.
find.leaves <- function (tree,i,side) {

    v <- tree[i,side]
    j <- 1
    l <- 1

    while (j <= l) {
        if (v[j] < 0) { # branch ends; move on
            j <- j+1
        } 
        else { # add a lower branch
            if (j < l) { 
                v[(j+2):(l+1)] <- v[(j+1):l]  # move over to make room 
            }
            v[j:(j+1)] <- tree[v[j],]  # put the new branch
        }
        l <- length(v)
    }
    return (-v)  # negative number should now be inverted
}


#----------------------------------------------------------------------
find.nodes <- function(tree,i,pv.cutoff) {
    v <- tree[i,]
    j <- 1
    l <- 2
    while (j<=l){
        if ((v[j] > 0) && (pv[v[j]] > pv.cutoff)) {
            v <- c(v,tree[v[j],])
        }
        j <- j+1
        l <- length(v)
    }
    return (v[v>0])
}


#----------------------------------------------------------------------
get.var <- function (y1) {
    mu <- if(!is.matrix(y1)) y1 else rowMeans(y1)
    y1.var <- (y1 - mu) %*% t(y1 - mu)
    y1.var
}


#----------------------------------------------------------------------
# Choose a metric for comparing within- and between-cluster variances.
get.metric <- function(method){
    switch(
        method,
        function(V) { # 1) Trace
            sum(diag(V))
        },
        function(V) { # 2) Product of singular values
            sv <- svd(V)$d
            sum(log(sv[which(sv > sv[1] * 1e-8)]))
        }
    )
}


#----------------------------------------------------------------------
permutation.function <- function(m, i, y, metric, ntrials=ntrials){
    left.branch  <- find.leaves(m,i,1)
    right.branch <- find.leaves(m,i,2)
    branch       <- c(left.branch,right.branch) 
    V0           <- get.var(y[, right.branch]) + get.var(y[, left.branch])
    V.r0         <- metric(V0)
    V.r <- replicate(ntrials, { 
        v1 <- sample(branch, length(right.branch))
        V0 <- get.var(y[,v1]) + get.var(y[, setdiff(branch,v1)])
        metric(V0)
    })
    sum(V.r <= V.r0) / length(V.r)
}


#----------------------------------------------------------------------
get.hclust <- function(y, method="average") {
    hclust(dist(t(y)), method=method)
}


#----------------------------------------------------------------------
calculate.cluster.pvalues <- function(y, hcl,  method=1, ntrials=100){
    n <- ncol(y)
    v <- hcl$order
    m <- hcl$merge
    pv <- within.var <- rep(0,(n-1))
    total.var <- get.var(y[, v])
    metric <- get.metric(method)
    for (i in 1:(n-1)) {
        if ( m[i, 1] < 0 && m[i, 2] < 0 ) {
            within.var[i] <- 0
            total.var[i] <- 1  # want the ratio to be 0
        } else {
            pv[i] <- permutation.function(m, i, y, metric, ntrials)
        }
    }
    pv
}



#======================================================================
# PLOTTING FUNCTIONS
#======================================================================

#----------------------------------------------------------------------
plot_cluster <- function(y){
    plot(y[1, ], y[2, ], xlab="x", ylab="y", type="n")
    text(y[1, ], y[2, ], labels=as.character(1:ncol(y)))
}


#----------------------------------------------------------------------
# plot the dendrogram and p-values
plot_confidence_tree <- function(hcl, pv){
    # (n-1)x2 matrix. Row i describes the merging of clusters at step i of the clustering
    m <- hcl$merge

    # a vector giving the permutation of the original observations suitable for plotting.
    v <- hcl$order	

    n <- length(v) # each column represents one patient

    # Find x-coordinates of labels. 
    x.c <- rep(0,n-1)
    for (i in 1:(n-1)){
        t1 <- m[i,1]
        t2 <- m[i,2]
        x1 <- if (t1 < 0) which(v == -t1) else x1 <- x.c[t1]
        x2 <- if (t2 < 0) which(v == -t2) else x2 <- x.c[t2]
        x.c[i] <- (x1 + x2) / 2
    }

    plot(hcl, axes=TRUE, xlab="", ylab=NULL, ann=FALSE)
    h.inc <- max(hcl$height) / n
    for (i in 1:(n-1)) {
        if (m[i,1]>0 || m[i,2]>0) {
            text(x.c[i] + 1.2, hcl$h[i] + h.inc, pv[i], cex=.7)
        }
    }
}


#----------------------------------------------------------------------
# adjust the branch lengths based on the p-values
plot_adjusted_confidence_tree <- function(hcl, pv){
    hcl1 <- hcl
    n <- length(hcl$order)
    pv.cutoff <- 0.05 / (n-1)
    for (i in 1:length(hcl1$height)) {
        if (pv[i] > pv.cutoff) {
            v <- c(i, find.nodes(hcl1$merge, i, pv.cutoff))
            for (j in 1:length(v)) {
                hcl1$height[v[j]] <- min(hcl1$height[v])
            }
        }
    }
    plot(hcl1, axes=TRUE, xlab="", ylab=NULL, ann=FALSE)
}



#======================================================================
# MISCELLANEOUS
#======================================================================
build.sample.dataset <- function(){
    require(MASS)
    sample.size <- 20
    sigma <- matrix(c(1, .9, .9, 1), 2, 2)
    a <- mvrnorm(n=sample.size, mu=c(0,0), Sigma=sigma)
    b <- mvrnorm(n=sample.size, mu=c(0,4), Sigma=sigma)
    y <- t(rbind(a, b))
    y
}
