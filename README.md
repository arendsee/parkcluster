# Cluster Permutation Test

This clustering algorithm was created by Park et. al.

Park, P. J., Manjourides, J., Bonetti, M., & Pagano, M. (2009). A permutation test for determining significance of clusters with applications to spatial and gene expression data. Computational statistics & data analysis, 53(12), 4290-4300.  

I have simply refactered it.

# Example
```
source('permutation_test_code.R')

y <- build.sample.dataset()
plot_cluster(y)

hcl <- get.hclust(y)
plot(hcl)

pvalues <- calculate.cluster.pvalues(y, hcl, method=2)
plot_confidence_tree(hcl, pvalues)
```
