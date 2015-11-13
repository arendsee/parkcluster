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



# =====================================================================
# Exported Functions
# =====================================================================

phclust_pvalues <- function(x, hcl=NULL, method='svd', ntrials=100){
    if(is.null(hcl)){
        hcl <- phclust(x)
    }
    n <- nrow(x)
    v <- hcl$order
    m <- hcl$merge
    pv <- within.var <- rep(0,(n-1))
    total.var <- get_var( t(x[v, ]) )
    metric <- get_metric(method)
    for (i in 1:(n-1)) {
        if ( m[i, 1] < 0 && m[i, 2] < 0 ) {
            within.var[i] <- 0
            total.var[i] <- 1  # want the ratio to be 0
        } else {
            pv[i] <- permutation_function(m, i, t(x), metric, ntrials)
        }
    }
    pv
}


phclust_plot <- function(x, pv=NULL, group=NULL, adjusted=FALSE, show.pval=!adjusted, ...){

    if(class(x) == 'hclust'){
        hcl <- x
        if(is.null(pv)){
            stop('Must provide pv vector if x is an hclust object')
        }
    } else {
        hcl <- phclust(x)
    }

    if(is.null(pv)){
        pv <- phclust_pvalues(x, hcl=hcl, ...)
    }

    # (n-1)x2 matrix. Row i describes the merging of clusters at step i of the clustering
    m <- hcl$merge

    # a vector giving the permutation of the original observations suitable for plotting.
    v <- hcl$order	

    # number of observations
    n <- length(v) 

    # Find x-coordinates of labels. 
    x.c <- rep(0,n-1)
    for (i in 1:(n-1)){
        t1 <- m[i,1]
        t2 <- m[i,2]
        x1 <- if (t1 < 0) which(v == -t1) else x1 <- x.c[t1]
        x2 <- if (t2 < 0) which(v == -t2) else x2 <- x.c[t2]
        x.c[i] <- (x1 + x2) / 2
    }

    # Adjust branch length according to p-value
    if(adjusted){
        pv.cutoff <- 0.05 / (n-1)
        for (i in 1:length(hcl$height)) {
            if (pv[i] > pv.cutoff) {
                v <- c(i, find_nodes(hcl$merge, pv, i, pv.cutoff))
                for (j in 1:length(v)) {
                    hcl$height[v[j]] <- min(hcl$height[v])
                }
            }
        }
    }

    if(!is.null(group)){
        group <- factor(group)
        col <- rainbow(nlevels(group))[group[v]]
    } else {
        col <- 'black'
    }
    plot(hcl, axes=TRUE, xlab="", ylab=NULL, ann=FALSE, col=col)

    # Show p-value next to each node
    if(show.pval){
        h.inc <- max(hcl$height) / n
        for (i in 1:(n-1)) {
            if (m[i,1]>0 || m[i,2]>0) {
                text(x.c[i] + 1.2, hcl$h[i] + h.inc, pv[i], cex=.7)
            }
        }
    }
}



# =====================================================================
# Internal Functions
# =====================================================================

phclust <- function(x, method="average") {
    hclust(dist(x), method=method)
}

#----------------------------------------------------------------------
# Given a tree ($merge) and index (ith splitting from the bottom) and 
# branch (left or right), get all the elements in that branch.  
# $merge contains a negative number if it is at the end; a positive 
# number points to a subsequent branch.  This subroutine updates
# until all the items are negative.
find_leaves <- function (tree,i,side) {

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
find_nodes <- function(tree, i, pv, pv.cutoff) {
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
get_var <- function (x1) {
    mu <- if(!is.matrix(x1)) x1 else rowMeans(x1)
    x1.var <- (x1 - mu) %*% t(x1 - mu)
    x1.var
}


#----------------------------------------------------------------------
# Choose a metric for comparing within- and between-cluster variances.
get_metric <- function(method){
    switch(
        method,
        trace = function(V) { # 1) Trace
            sum(diag(V))
        },
        svd = function(V) { # 2) Product of singular values
            sv <- svd(V)$d
            sum(log(sv[which(sv > sv[1] * 1e-8)]))
        }
    )
}


#----------------------------------------------------------------------
permutation_function <- function(m, i, x, metric, ntrials=ntrials){
    left.branch  <- find_leaves(m,i,1)
    right.branch <- find_leaves(m,i,2)
    branch       <- c(left.branch,right.branch) 
    V0           <- get_var(x[, right.branch]) + get_var(x[, left.branch])
    V.r0         <- metric(V0)
    V.r <- replicate(ntrials, { 
        v1 <- sample(branch, length(right.branch))
        V0 <- get_var(x[,v1]) + get_var(x[, setdiff(branch,v1)])
        metric(V0)
    })
    sum(V.r <= V.r0) / length(V.r)
}
