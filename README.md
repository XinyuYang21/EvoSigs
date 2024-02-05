<<<<<<< HEAD
# EvoSig

EvoSig is an R packge for quantifying evolutionary dynamics signatures from bulk whole genome/exome data. This repository also contains all code and documentation necessary to reproduce the analysis undertaken in the manuscript [Pan-cancer evolution signatures link clonal expansion to dynamic changes in the tumour immune microenvironment](https://www.biorxiv.org/content/10.1101/2023.10.12.560630v1.abstract).

## Installation
=======
# EvoSigs
EvoSig is an R package for quantifying evolutionary dynamics signatures from bulk whole genome/exome data (TCGA/PCAWG). This repository also contains all code and documentation necessary to reproduce the analysis undertaken in the manuscript 

### Data and result files
### How to estimate existing evolutionary dynamics signatures in your own samples
Once you have formatted your mutation data (VCF/MAF etc.,) into CCF distribution matrices, run the following functions: 
-
>>>>>>> 18b6d22c1311852450904ffdf2075f4d62587564

```{r}
require(devtools)
devtools::install_github("XinyuYang21/EvoSigs")
```

<<<<<<< HEAD
## Getting Started

### How to estimate existing evolutionary dynamics signatures in your own samples

Follow Rmarkdown manuscript `Extract_signature.Rmd`

```{r}
manuscript_RMD/Extract_signature.Rmd
```

### How to generate consensus sigantures of evolutionary dynamics

Follow Rmarkdown manuscript `Esigs_identification.Rmd` and rank estimate analysis

```{r}
manuscript_RMD/Esigs_identification.Rmd
manuscript_RMD/rank_estimate
=======
>>>>>>> 18b6d22c1311852450904ffdf2075f4d62587564
```
