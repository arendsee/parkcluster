## Copyright 2001,2009
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

# =====================================================================
# Exported Functions
# =====================================================================

#' Calculate p-values for each branch in the tree
#' 
#' @export
#' @param x A data.frame or matrix with one row for each individual
#' @param hcl A "hclust" distance object, it will be created if missing.
#' @param method Either 'trace' or 'svd', see Park 2009
#' @param nperm The number of permutations to use for predicting branch confidence
#' @return A vector of p-values
#' @examples
#' phclust_pvalues(iris[1:3], nperm=10)
phclust_pvalues <- function(x, hcl=phclust(x), method='svd', nperm=100){
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
            pv[i] <- permutation_function(m, i, t(x), metric, nperm)
        }
    }
    pv
}


#' Plot the hierarchical cluster significance tree
#'
#' @export
#' @param x A data.frame, matrix, or hclust object
#' @param pv A vector of p-values, an output of phclust_pvalues
#' @param group A factor defining a hypothesized group, these groups will be colored on the tree, but will not affect the main algorithm
#' @param cutoff Minimum p-value required to resolved a branch
#' @param show.pval Whether to show the branch p-values on the tree
#' @param ... Arguments that will be passed to phclust_pvalues (ignored if pv is given)
#' @return A tree plot
#' @examples
#' phclust_plot(cars)
#' phclust_plot(iris[1:3], group=iris$Species, cutoff=0.1, nperm=10)
phclust_plot <- function(x, pv=NULL, group=NULL, cutoff=1, show.pval=TRUE, ...){
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

    # Adjust branch length according to p-value
    for (i in which(pv > cutoff)) {
        nodes <- c(i, find_nodes(m, i, pv, cutoff))
        str(nodes)
        hcl$height[nodes] <- min(hcl$height[nodes])
    }

    col <- if(is.null(group)) 1 else as.numeric(group) 
    tre <- hcl %>%
        as.dendrogram %>%
        dendextend::set('labels_col', col) %>%
        plot

    # Show p-value next to each node
    if(show.pval){
        # Find x-coordinates of labels. 
        x.c <- rep(0,n-1)
        for (i in 1:(n-1)){
            t1 <- m[i,1]
            t2 <- m[i,2]
            x1 <- if (t1 < 0) which(v == -t1) else x.c[t1]
            x2 <- if (t2 < 0) which(v == -t2) else x.c[t2]
            x.c[i] <- (x1 + x2) / 2
        }
        # Add p-values to significant branches
        i <- (m[, 1] > 0 | m[, 2] > 0) & pv <= cutoff
        xpos <- x.c[i] + 1.2
        ypos <- hcl$height[i] + max(hcl$height) / n
        text(xpos, ypos, pv[i], cex=.7)
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
find_nodes <- function(tree, i, pv, cutoff) {
    v <- tree[i,]
    j <- 1
    l <- 2
    while (j<=l){
        if ((v[j] > 0) && (pv[v[j]] > cutoff)) {
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
permutation_function <- function(m, i, x, metric, nperm=nperm){
    left.branch  <- find_leaves(m,i,1)
    right.branch <- find_leaves(m,i,2)
    branch       <- c(left.branch,right.branch) 
    V0           <- get_var(x[, right.branch]) + get_var(x[, left.branch])
    V.r0         <- metric(V0)
    V.r <- replicate(nperm, { 
        v1 <- sample(branch, length(right.branch))
        V0 <- get_var(x[,v1]) + get_var(x[, setdiff(branch,v1)])
        metric(V0)
    })
    sum(V.r <= V.r0) / length(V.r)
}
