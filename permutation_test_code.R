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


source("subroutines.R")
require(MASS)


##Must choose the proper method (Trace or Singular Value Decomposition), as described in the paper.
method <- 2 # 1: trace; 2: product of singular values


##Example Data Set
random.samples <- 100

  sample.size <- 20
  sigma <- matrix(c(1, .9, .9, 1), 2, 2)
  a <- mvrnorm(n=sample.size, mu=c(0,0), Sigma=sigma)
  b <- mvrnorm(n=sample.size, mu=c(0,4), Sigma=sigma)
  y <- t(rbind(a, b))
  plot(y[1,],y[2,],xlab="x",ylab="y",type="n")
  text(y[1,],y[2,],labels=as.character(1:(dim(y)[2])))

  plot(hclust(dist(t(y)), method="ave"))
#----------------------------------------------------------------------

n <- ncol(y) # each column represents one patient
d <- matrix(0,n,n)


hcl <- hclust(dist(t(y)), method="average")

# (n-1)x2 matrix. Row i describes the merging of clusters at step i of the clustering
m <- hcl$merge

# a vector giving the permutation of the original observations suitable for plotting.
v <- hcl$order	

#----------------------------------------------------------------------
# Find x-coordinates of labels.  Start building this information from
# the bottom to minimize work.

x.c <- get.coord(m,v,rep(0,n-1),n)
pv <- within.var <- rep(0,(n-1))
total.var <- get.var(y[,v])


#----------------------------------------------------------------------
#          M  A  I  N     L  O  O  P
#----------------------------------------------------------------------

# Choose a metric for comparing within- and between-cluster variances.
W.metric <- switch(
    method,
    function(V) { # 1) Trace
        sum(diag(V))
    },
    function(V) { # 2) Product of singular values
        sv <- svd(V)$d
        sum(log(sv[which(sv > sv[1] * 1e-8)]))
    }
)

permutation_function <- function(m, i, y, ntrials=100){
    left.branch  <- find.leaves(m,i,1)
    right.branch <- find.leaves(m,i,2)
    branch       <- c(left.branch,right.branch) 
    V0           <- get.var(y[, right.branch]) + get.var(y[, left.branch])
    V.r0         <- W.metric(V0)
    V.r <- replicate(ntrials, { 
        v1 <- sample(branch, length(right.branch))
        V0 <- get.var(y[,v1]) + get.var(y[, setdiff(branch,v1)])
        W.metric(V0)
    })
    sum(V.r <= V.r0) / length(V.r)
}

for (i in 1:(n-1)) {
    if ( m[i, 1] < 0 && m[i, 2] < 0 ) {
        within.var[i] <- 0
        total.var[i] <- 1  # want the ratio to be 0
    } else {
        pv[i] <- permutation_function(m, i, y, ntrials=random.samples)
        cat("i = ", i, " p.values = ", pv[i], "\n")
    }
}


#----------------------------------------------------------------------
# plot the dendrogram and p-values
#----------------------------------------------------------------------
plot(hcl, axes=TRUE, xlab="", ylab=NULL, ann=FALSE)
h.inc <- max(hcl$height) / n
for (i in 1:(n-1)) {
  if (m[i,1]>0 || m[i,2]>0) {
    text(x.c[i] + 1.2, hcl$h[i] + h.inc, pv[i], cex=.7)
  }
}


#----------------------------------------------------------------------
# adjust the branch lengths based on the p-values
#----------------------------------------------------------------------

hcl1 <- hcl
pv.cutoff <- 0.05 / (n-1)
#pv.cutoff <- .005
for (i in 1:length(hcl1$height)) {
  if (pv[i] > pv.cutoff) {
    v <- c(i, find.nodes(hcl1$merge, i, pv.cutoff))
    for (j in 1:length(v)) {
      hcl1$height[v[j]] <- min(hcl1$height[v])
    }
  }
}

plot(hcl1, axes=TRUE, xlab="", ylab=NULL, ann=FALSE)
