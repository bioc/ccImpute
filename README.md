# ccImpute 
<a><img src="https://www.bioconductor.org/shields/years-in-bioc/ccImpute.svg" alt="Years in BioConductor badge" title="How long since the package was first in a released Bioconductor version (or is it in devel only)."/> </a>
<a href="https://bioconductor.org/checkResults/release/bioc-LATEST/ccImpute/"><img src="https://www.bioconductor.org/shields/build/release/bioc/ccImpute.svg" alt="Build results badge" title="build results; click for full report"/></a>
<a href='https://www.gnu.org/licenses/gpl-3.0'><img src="https://img.shields.io/badge/License-GPLv3-blue.svg"/></a>
<a href='https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04814-8'><img src="https://img.shields.io/badge/Manuscript-gray"/></a>
<a href='https://www.bioconductor.org/packages/release/bioc/vignettes/ccImpute/inst/doc/ccImpute.html'><img src="https://img.shields.io/badge/Package-Manual-blue"/></a>
<a href='https://www.bioconductor.org/packages/release/bioc/manuals/ccImpute/man/ccImpute.pdf'><img src="https://img.shields.io/badge/Reference-Manual-blue"/></a>


<a href='https://www.bioconductor.org/packages/release/bioc/html/ccImpute.html'><img src='man/figures/logo.png' align="right" height="180"/></a>


In single-cell RNA sequencing (scRNA-seq) data, dropout events often mask lowly expressed genes, making them indistinguishable from true zero expression and obscuring subtle differences between cell types. This hinders accurate downstream analysis. ccImpute is an innovative imputation algorithm designed to address this issue. It leverages consensus clustering to identify similar cells and imputes the most probable dropout events. The rigorous evaluation demonstrates that ccImpute outperforms existing methods while minimizing the introduction of noise, as evidenced by clustering performance on datasets with known cell identities. Please see the [manuscript](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04814-8) for more details.

## Installation
You can install the release version of
*[ccImpute](https://www.bioconductor.org/packages/release/bioc/html/ccImpute.html)*
from BioConductor:
``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ccImpute")
```
Alternatively, you can install the development version of the package from [GitHub](https://github.com/khazum/ccImpute) with:
``` r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("khazum/ccImpute")
```
## Quick Start
Here's a simple example that loads scRNA-seq data and imputes dropout events, creating a new assay named 'imputed' within the 'SingleCellExperiment' object. Please refer to the package manual for detailed instructions and strategies to optimize the performance. [package manual](https://www.bioconductor.org/packages/release/bioc/vignettes/ccImpute/inst/doc/ccImpute.html).
``` r
library(scRNAseq)
library(scater)
library(BiocParallel)

# Load Usoskin brain cells dataset:
sce <- UsoskinBrainData()

X <- cpm(sce)
labels <- colData(sce)$"Level 1"

# Do pre-processing of the data:
#Filter bad cells
filt <- !grepl("Empty well", labels) &
        !grepl("NF outlier", labels) &
        !grepl("TH outlier", labels) &
        !grepl("NoN outlier", labels) &
        !grepl("NoN", labels) &
        !grepl("Central, unsolved", labels) &
        !grepl(">1 cell", labels) &
        !grepl("Medium", labels)

labels <-labels[filt]
X <- as.matrix(X[,filt])

#Remove genes that are not expressed in any cells:
X <- X[rowSums(X)>0,]

#Recreate the SingleCellExperiment and add log-transformed data:
ann <- data.frame(cell_id = labels)
sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(X)), 
                               colData = ann)
logcounts(sce) <- log(normcounts(sce) + 1)

cores <- 2
BPPARAM = MulticoreParam(cores)
sce <- ccImpute(sce, BPPARAM=BPPARAM)
summary(assay(sce, "imputed"))
```
## Release Notes
### Version 1.6.1
- Performance Optimizations:
    - Significantly enhanced calculation speed for Pearson and Spearman 
      correlation matrices, including weighted versions.
    - Leveraged the Irlba package for efficient truncated Singular Value 
      Decomposition (SVD) computation.
    - Optimized imputation by limiting the number of singular components while 
      maintaining the accuracy of downstream analysis, with adjustable maximum 
      limits based on dataset size.
    - Optimized the identification of dropout events.
    - Introduced a fast dropout calculation method based on non-zero expression
      value means, preserving imputation performance and greatly improving 
      runtime speed.
    - Replaced SIMLR with Tracy-Widom Bound for estimating k when not provided,
      resulting in faster calculations and improved empirical performance.
- Expanded Functionality:
    - Added support for sparse matrices in dgCmatrix format, allowing increased memory 
      efficiency.
- Documentation Enhancements:
    - Expanded the package manual with detailed guidance and practical examples for
      maximizing the package's value and computational speed.
    - Included comparative benchmarking against previous release in the
      package manual, demonstrating the performance improvements.
- Overall Impact:
    - The ccImpute package is now substantially faster and more efficient.
    - Users can expect a smoother experience with improved documentation and
      expanded functionality.
### Version 0.99.x:
- The initial version was submitted to Bioconductor.

## Issues
Please use [this page](https://github.com/khazum/ccImpute/issues) to report bugs, comments and suggestions.
