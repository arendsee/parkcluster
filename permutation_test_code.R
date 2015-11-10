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


source("subroutines.R")

##Example Data Set
random.samples <- 100

  sample.size <- 20
  sigma <- matrix(c(1,.9,.9,1),2,2)
  a <- rmultnorm(sample.size,c(0,0),sigma)
  b <- rmultnorm(sample.size,c(0,4),sigma)
  y <- t(rbind(a,b))
  plot(y[1,],y[2,],xlab="x",ylab="y",type="n")
  text(y[1,],y[2,],labels=as.character(1:(dim(y)[2])))

  plclust(hclust(dist(t(y)),method="ave"))
#----------------------------------------------------------------------

num.genes <- dim(y)[1]      # each row is for a gene
n <- dim(y)[2]              # each column is for a patient
d <- matrix(0,n,n)


hcl <- hclust(dist(t(y)), method="average")
m <- hcl$merge	# (n-1)x2 matrix. Row i describes the merging of clusters at step i of the clustering

v <- hcl$order	# a vector giving the permuation of the original observations suitable for plotting.

#----------------------------------------------------------------------
# Find x-coordinates of labels.  Start building this information from
# the bottom to minimize work.

x.c <- get.coord(m,v,rep(0,n-1),n)
pv <- within.var <- rep(0,(n-1))
total.var <- get.var(y[,v])




#----------------------------------------------------------------------
#          M  A  I  N     L  O  O  P
#----------------------------------------------------------------------

##Must choose the proper method (Trace or Singular Value Decomposition), as described in the paper.

method <- 2 # 1: trace; 2: product of singular values

for (i in 1:(n-1)) {

  if ( m[i,1]<0 && m[i,2]<0 ) {
    within.var[i] <- 0
    total.var[i] <- 1  # want the ratio to be 0
  }
  else {
    left.branch <- find.leaves(m,i,1)
    right.branch <- find.leaves(m,i,2)
    branch <- c(left.branch,right.branch) 
    V0  <- get.var(y[,right.branch]) + get.var(y[,left.branch])
    if (method==1) {
      V.r0 <- sum(diag(V0)) 
    }
    else {
      sv <- svd(V0)$d
      V.r0 <- sum(log(sv[which(sv>sv[1]*1e-8)]))
    }
    V.r <- rep(0,random.samples)
    lr <- length(right.branch)

    for (j in 1:random.samples) { 
      v1 <- sample(branch,lr)
      V0 <- get.var(y[,v1]) + get.var(y[,setdiff(branch,v1)])

      if (method==1) {
        V.r[j] <- sum(diag(V0))
      }
      else {
        sv <- svd(V0)$d
        V.r[j] <- sum(log(sv[which(sv>sv[1]*1e-8)]))
      }
    }
    pv[i] <- sum(V.r <= V.r0) /length(V.r)
    cat("i = ",i," p.values = ",pv[i],"\n")
  }
}


#----------------------------------------------------------------------
# plot the dendrogram and p-values
#----------------------------------------------------------------------
plclust(hcl,axes=F,xlab="",ylab=NULL,ann=F,labels=FALSE)
h.inc <- max(hcl$height)/n
for (i in 1:(n-1)) {
  if (m[i,1]>0 || m[i,2]>0) {
    text(x.c[i]+1.2,hcl$h[i]+h.inc,pv[i],cex=.7)
  }
}


#----------------------------------------------------------------------
# adjust the branch lengths based on the p-values
#----------------------------------------------------------------------

hcl1 <- hcl
pv.cutoff <- .05/(n-1)
#pv.cutoff <- .005
for (i in 1:length(hcl1$height)) {
  if (pv[i]>pv.cutoff) {
    v <- c(i,find.nodes(hcl1$merge,i,pv.cutoff))
    for (j in 1:length(v)) {
      hcl1$height[v[j]] <- min(hcl1$height[v])
    }
  }
}

plot(hcl1,axes=T,xlab="",ylab=NULL,ann=F)




