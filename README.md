# GeLaToLab
A collection of visualization tools and algorithms for data analysis and machine learning focused on source code analysis.

## Installation

```r
# Install devtools from CRAN
install.packages("devtools")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("amirms/GeLaToLab")
```

## Multi-view Learning

### Multi-view Clustering
Function `perform.clustering` in file `clusterValidate.R` performs evaluation of single and multi-view clustering for each project. It finds the best kernel functions + parameters for each view, and uses the same choices to perform multi-view clustering for different methods.

For clustering, `complete` hierarchical cluster analysis is used, and path difference (PD) is used to measure the performance.

### Multi-view CF-based Recommendation System
Function `perform.prediction` in file `recommenderValidate.R` perform single-view and multi-view evaluation of CF-based recommendations for each project using a nested k-fold cross validation. The default setup is a 10-fold nested cross-validation.

It first finds the best kernel parameters for each view, followed by evaluating three multi-view learning approaches to perform CF-based recommendation.

The scores used to measure the performance are: `ROC AUC`, `PR AUC`, and `max F1` scores.

### Uni- and Cross-modal search

In `crossModalRetrieval.R` file, function `findCrossModal` performs the experiment for uni-modal and cross-modal retrieval based on the best kernel parameters for each view.






