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
#### Loading the package and running the functions
```
library("vcfGrapher")

data <- "Your_genotypes_file.vcf"

read_depth_per_position(vcf_name = data)
\
SNP_heterosigousity_per_position(vcf_name = data)

SNP_heterosigousity_histogram(vcf_name = data)
Samples_heterosigousity_histogram(vcf_name = data)
Samples_heterosigousity_boxplot(vcf_name = data)
Missing_data_per_position(vcf_name = data)
Samples_proportion_missing_histogram(vcf_name = data)
Samples_proportion_missing_boxplot(vcf_name = data)

```


<br>

#### An example of the output produced by the read_depth_per_position function

![Alt text](/images/Read_depth_per_genomic_position.png)

