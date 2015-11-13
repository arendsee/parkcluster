# Cluster Permutation Test

This clustering algorithm was created by Park et. al.

Park, P. J., Manjourides, J., Bonetti, M., & Pagano, M. (2009). A permutation test for determining significance of clusters with applications to spatial and gene expression data. Computational statistics & data analysis, 53(12), 4290-4300.  

I have simply refactered it.

# Installation

```R
library(devtools)
install_github('zbwrnz/parkcluster')
```

# Example
```R
library(parkcluster)

?phcluster_plot
?phcluster_pvalues

phclust_plot(cars)

d <- iris[5*(1:30), c('Sepal.Length', 'Sepal.Width', 'Species')]

plot(d[[1]], d[[2]], col=as.numeric(d$Species))

phclust_plot(d[1:2], group=d$Species, cutoff=0.1, ntrials=1000)
```
