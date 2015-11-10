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


#----------------------------------------------------------------------
get.coord <- function(m,v,x.coord,n) {
  for (i in 1:(n-1))  {
    t1 <- m[i,1]
    t2 <- m[i,2]
    x1 <- if (t1 < 0) which(v == -t1) else x1 <- x.coord[t1]
    x2 <- if (t2 < 0) which(v == -t2) else x2 <- x.coord[t2]
    x.coord[i] <- (x1 + x2) / 2
  }
  x.coord
}


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
