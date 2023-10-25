# vcfGrapher
An R package that generates plots directly from a vcf file

# Installation
To use vcfGrapher locally you can either install the package (method 1) on your machine or source it (method 2)

Installation method 1
```
library(devtools)
install_github("tomstamour/vcfGrapher")
```
Installation method 2 (devtools not needed)
```
source("https://raw.githubusercontent.com/tomstamour/vcfGrapher/main/R/vcfGrapher.R")
```
# Usage
library("vcfGrapher")
read_depth_per_position(vcf_name = "Genotype_file.vcf")
![Alt text](/posts/path/to/Read_depth_per_genomic_position.png)
