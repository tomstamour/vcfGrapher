#################
# vcfGrapher
#################

# setwd("Labo_Belzile/Autre/scripts/vcfGrapher/input_files/")
# vcf_name <- "FastGBS_platypus_noScaf_onlyPASS_rmIndels_minMaxAll2_mac4_max_missing02.vcf"

#-----------------------------------------------------------------

#                 READ DEPTH per genomic position

#-----------------------------------------------------------------
# Plotting the average read depth as a function of genomic position

#' Read depth per genomic position
#'
#' Plotting read depth per genomic position
#' @param vcf_name The name of your vcf file
#' @return A plot for each chromosome with the position of the SNP in base pair on the x axis and the average number of reads on the y axis.
#' @export
read_depth_per_position <- function(vcf_name) {

  # Loading the necessary packages
  library(data.table)

  # Reading the VCF file using readLines
  vcf_lines <- readLines(vcf_name)

  # Filtering out header lines (lines that start with "##")
  header_lines <- vcf_lines[grepl("^##", vcf_lines)]

  # Extracting the column names from the header line starting with "#CHROM"
  column_names <- unlist(strsplit(vcf_lines[grepl("^#CHROM", vcf_lines)], "\t"))

  # Reading the VCF data using data.table, skipping the header lines
  vcf_data <- fread(vcf_name, skip = length(header_lines) + 1, sep = "\t", header = FALSE, fill = TRUE)
  colnames(vcf_data)[1] <- "CHROM"

  # Seting the column names
  colnames(vcf_data) <- column_names
  colnames(vcf_data)[1] <- "CHROM"

  # outputing the information in the vcf file ("FORMAT" column)
  format_infos <- vcf_data[1, "FORMAT"]

  # creating a matrix
  vcf_matrix <- as.matrix(vcf_data[ , 10:ncol(vcf_data)])

  # Extracting the data from the matrix
  sub_mat <- strsplit(vcf_matrix, ":")
  infos <- unlist(strsplit(as.character(format_infos), split = ":"))
  data_pos <- grep("NR|DP", infos)

  if (!identical(data_pos, integer(0))) {

    vector <- sapply(sub_mat, "[", data_pos)
    data_matrix <- matrix(vector, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))
  }

  # Dealing with triallelic sites which have read depth coded as x,x
  sub_mat2 <- strsplit(data_matrix, ",")
  vector2 <- sapply(sub_mat2, "[", 1)
  data_matrix2 <- matrix(vector2, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))

  # converting the character matrix to a numeric matrix
  data_matrix2 <- matrix(as.numeric(data_matrix2), ncol = ncol(data_matrix2))

  # Calculating the average read depth per position
  average_data_per_position <- data.frame(Chromosome = vcf_data$CHROM, Position = vcf_data$POS)
  average_data_per_position$average_data <- rowMeans(data_matrix2)

  # Plotting the read depth as a function of genomic position
  library(ggplot2)
  ggplot(data = average_data_per_position, aes(x = Position, y = average_data)) +
    geom_point() +
    facet_wrap(~Chromosome) +
    labs(y= "Average read depth", x = "Position (bp)") +
    ggtitle("Average read depth as a function of genomic position") +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), plot.title = element_text(size = 25)) +
    labs(caption = paste("Input file:", vcf_name))
}

# To run the function
# read_depth_per_position(vcf_name = "FastGBS_platypus_noScaf_onlyPASS_rmIndels_minMaxAll2_mac4_max_missing02.vcf")



#-----------------------------------------------------------------

#                 HETEROZYGOUSITY per SNP per genomic position

#-----------------------------------------------------------------
# Plotting the proportion of heterosigousity per SNP (or any variant) as a function of the position on the genome

#' Heterozigousity per genomic position
#'
#' Plotting the proportion of heterosigousity per SNP on the y axis as a function of the genomic position of the SNP
#' @param vcf_name The name of your vcf file
#' @return A plot of heterosigousity ratio per SNP as a function of the position on the genome
#' @export
SNP_heterosigousity_per_position <- function(vcf_name) {

  # Loading the necessary packages
  library(data.table)
  library(ggplot2)

  # Reading the VCF file using readLines
  vcf_lines <- readLines(vcf_name)

  # Filtering out header lines (lines that start with "##")
  header_lines <- vcf_lines[grepl("^##", vcf_lines)]

  # Extracting the column names from the header line starting with "#CHROM"
  column_names <- unlist(strsplit(vcf_lines[grepl("^#CHROM", vcf_lines)], "\t"))

  # Reading the VCF data using data.table, skipping the header lines
  vcf_data <- fread(vcf_name, skip = length(header_lines) + 1, sep = "\t", header = FALSE, fill = TRUE)
  colnames(vcf_data)[1] <- "CHROM"

  # Seting the column names
  colnames(vcf_data) <- column_names
  colnames(vcf_data)[1] <- "CHROM"

  # outputing the information in the vcf file ("FORMAT" column)
  format_infos <- vcf_data[1, "FORMAT"]

  # creating a matrix
  vcf_matrix <- as.matrix(vcf_data[ , 10:ncol(vcf_data)])

  # Extracting the data from the matrix
  sub_mat <- strsplit(vcf_matrix, ":")
  infos <- unlist(strsplit(as.character(format_infos), split = ":"))
  data_pos <- grep("GT", infos)

  if (!identical(data_pos, integer(0))) {

    vector <- sapply(sub_mat, "[", data_pos)
    data_matrix <- matrix(vector, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))
  }

  # Dealing with triallelic sites which have read depth coded as x,x
  sub_mat2 <- strsplit(data_matrix, ",")
  vector2 <- sapply(sub_mat2, "[", 1)
  data_matrix2 <- matrix(vector2, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))

  # Creating a data frame to store the ratio of heterosigousity per SNP (variant)
  data_genotypes_perSNPs <- data.frame(CHROM = vcf_data$CHROM, Position = vcf_data$POS)
  data_genotypes_perSNPs$Proportion_of_Hetero <-  rowSums(data_matrix2 == "1/0" | data_matrix2 == "0/1", na.rm=TRUE)/(ncol(data_matrix2) - rowSums(data_matrix2 == "./.", na.rm=TRUE))

  # Plotting the proportion of heterozigousity per SNP (variant)
  ggplot(data = data_genotypes_perSNPs, aes(x = Position, y = Proportion_of_Hetero)) +
    geom_point() +
    facet_wrap(~CHROM)+
    labs(y= "Proportion heterosigousity", x = "Position (bp)") +
    ggtitle("Proportion heterosigousity as a function of genomic position") +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), plot.title = element_text(size = 25)) +
    labs(caption = paste("Input file:", vcf_name))
}

# To run the function:
# SNP_heterosigousity_per_position(vcf_name = "FastGBS_platypus_noScaf_onlyPASS_rmIndels_minMaxAll2_mac4_max_missing02.vcf")


#-----------------------------------------------------------------

#                 HETEROZYGOUSITY per SNP histogram

#-----------------------------------------------------------------
# Plotting an histogram of the proportion of heterosigousity per SNP (or any variant)

#' Histogram of the proportion of heterozigousity per SNP
#'
#' Plotting an histogram for the proportion of heterosigousity per SNP
#' @param vcf_name The name of your vcf file
#' @return An histogram of heterosigousity ratio per SNP
#' @export
SNP_heterosigousity_histogram <- function(vcf_name) {

  # Loading the necessary packages
  library(data.table)
  library(ggplot2)

  # Reading the VCF file using readLines
  vcf_lines <- readLines(vcf_name)

  # Filtering out header lines (lines that start with "##")
  header_lines <- vcf_lines[grepl("^##", vcf_lines)]

  # Extracting the column names from the header line starting with "#CHROM"
  column_names <- unlist(strsplit(vcf_lines[grepl("^#CHROM", vcf_lines)], "\t"))

  # Reading the VCF data using data.table, skipping the header lines
  vcf_data <- fread(vcf_name, skip = length(header_lines) + 1, sep = "\t", header = FALSE, fill = TRUE)
  colnames(vcf_data)[1] <- "CHROM"

  # Seting the column names
  colnames(vcf_data) <- column_names
  colnames(vcf_data)[1] <- "CHROM"

  # outputing the information in the vcf file ("FORMAT" column)
  format_infos <- vcf_data[1, "FORMAT"]

  # creating a matrix
  vcf_matrix <- as.matrix(vcf_data[ , 10:ncol(vcf_data)])

  # Extracting the data from the matrix
  sub_mat <- strsplit(vcf_matrix, ":")
  infos <- unlist(strsplit(as.character(format_infos), split = ":"))
  data_pos <- grep("GT", infos)

  if (!identical(data_pos, integer(0))) {

    vector <- sapply(sub_mat, "[", data_pos)
    data_matrix <- matrix(vector, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))
  }

  # Dealing with triallelic sites which have read depth coded as x,x
  sub_mat2 <- strsplit(data_matrix, ",")
  vector2 <- sapply(sub_mat2, "[", 1)
  data_matrix2 <- matrix(vector2, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))

  # Creating a data frame to store the ratio of heterosigousity per SNP (variant)
  data_genotypes_perSNPs <- data.frame(CHROM = vcf_data$CHROM, Position = vcf_data$POS)
  data_genotypes_perSNPs$Proportion_of_HomoRef <-  rowSums(data_matrix2 == "0/0", na.rm=TRUE)/ncol(data_matrix2)
  data_genotypes_perSNPs$Proportion_of_HomoAlt <-  rowSums(data_matrix2 == "1/1", na.rm=TRUE)/ncol(data_matrix2)
  data_genotypes_perSNPs$Proportion_of_Hetero <-  rowSums(data_matrix2 == "1/0" | data_matrix2 == "0/1", na.rm=TRUE)/(ncol(data_matrix2) - rowSums(data_matrix2 == "./.", na.rm=TRUE))
  data_genotypes_perSNPs$Proportion_of_Missing <-  rowSums(data_matrix2 == "./.", na.rm=TRUE)/ncol(data_matrix2)

  # Plotting an histogram of the proportion of heterozigousity per SNP (variant)
  ggplot(data = data_genotypes_perSNPs, aes(x = Proportion_of_Hetero)) +
    geom_histogram() +
    labs(y= "Count", x = "Proportion of heterosigousity") +
    ggtitle("Distribution of proportion heterosigousity per SNP") +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), plot.title = element_text(size = 25)) +
    labs(caption = paste("Input file:", vcf_name))
}

# To run the function:
# SNP_heterosigousity_histogram(vcf_name = "FastGBS_platypus_noScaf_onlyPASS_rmIndels_minMaxAll2_mac4_max_missing02.vcf")


#-----------------------------------------------------------------

#               HETEROZYGOUSITY per sample histogram

#-----------------------------------------------------------------
# Plotting an histogram of the proportion of heterosigousity per sample
Samples_heterosigousity_histogram <- function(vcf_name, Pattern_to_highlight_in_sample_name) {

  # Loading the necessary packages
  library(data.table)
  library(ggplot2)

  # Reading the VCF file using readLines
  vcf_lines <- readLines(vcf_name)

  # Filtering out header lines (lines that start with "##")
  header_lines <- vcf_lines[grepl("^##", vcf_lines)]

  # Extracting the column names from the header line starting with "#CHROM"
  column_names <- unlist(strsplit(vcf_lines[grepl("^#CHROM", vcf_lines)], "\t"))

  # Reading the VCF data using data.table, skipping the header lines
  vcf_data <- fread(vcf_name, skip = length(header_lines) + 1, sep = "\t", header = FALSE, fill = TRUE)
  colnames(vcf_data)[1] <- "CHROM"

  # Seting the column names
  colnames(vcf_data) <- column_names
  colnames(vcf_data)[1] <- "CHROM"

  # outputing the information in the vcf file ("FORMAT" column)
  format_infos <- vcf_data[1, "FORMAT"]

  # creating a matrix
  vcf_matrix <- as.matrix(vcf_data[ , 10:ncol(vcf_data)])

  # Extracting the data from the matrix
  sub_mat <- strsplit(vcf_matrix, ":")
  infos <- unlist(strsplit(as.character(format_infos), split = ":"))
  data_pos <- grep("GT", infos)

  if (!identical(data_pos, integer(0))) {

    vector <- sapply(sub_mat, "[", data_pos)
    data_matrix <- matrix(vector, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))
  }

  # Dealing with triallelic sites which have read depth coded as x,x
  sub_mat2 <- strsplit(data_matrix, ",")
  vector2 <- sapply(sub_mat2, "[", 1)
  data_matrix2 <- matrix(vector2, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))

  # creating a vector of heterosigousity ratio for all samples
  by_samples_metrics <- data.frame(Sample_name = colnames(vcf_data)[10:ncol(vcf_data)])
  by_samples_metrics$Proportion_heterosigousity <- colSums(data_matrix2 == "1/0" | data_matrix2 == "0/1", na.rm = FALSE) / (nrow(data_matrix2) - colSums(data_matrix2 == "./.", na.rm = FALSE))

  # Plotting an histogram of the proportion of heterozigousity per SNP (variant)
  ggplot(data = by_samples_metrics, aes(x = Proportion_heterosigousity)) +
    geom_histogram() +
    labs(y= "Count", x = "Proportion of heterosigousity") +
    ggtitle("Distribution of proportion heterosigousity per sample") +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), plot.title = element_text(size = 25)) +
    labs(caption = paste("Input file:", vcf_name))
}

# To run the function:
# Samples_heterosigousity_histogram(vcf_name = "FastGBS_platypus_noScaf_onlyPASS_rmIndels_minMaxAll2_mac4_max_missing02.vcf")



#-----------------------------------------------------------------

#               Box plot HETEROZYGOUSITY per sample

#-----------------------------------------------------------------
# Plotting a boxplot of the proportion of heterosigousity per sample
# Pattern_to_highlight_in_sample_name : will color the dots on the boxplot regarding the patern has been found in the name or not
Samples_heterosigousity_boxplot <- function(vcf_name, Pattern_to_highlight_in_sample_name) {

  # Loading the necessary packages
  library(data.table)
  library(ggplot2)

  # Reading the VCF file using readLines
  vcf_lines <- readLines(vcf_name)

  # Filtering out header lines (lines that start with "##")
  header_lines <- vcf_lines[grepl("^##", vcf_lines)]

  # Extracting the column names from the header line starting with "#CHROM"
  column_names <- unlist(strsplit(vcf_lines[grepl("^#CHROM", vcf_lines)], "\t"))

  # Reading the VCF data using data.table, skipping the header lines
  vcf_data <- fread(vcf_name, skip = length(header_lines) + 1, sep = "\t", header = FALSE, fill = TRUE)
  colnames(vcf_data)[1] <- "CHROM"

  # Setting the column names
  colnames(vcf_data) <- column_names
  colnames(vcf_data)[1] <- "CHROM"

  # outputing the information in the vcf file ("FORMAT" column)
  format_infos <- vcf_data[1, "FORMAT"]

  # creating a matrix
  vcf_matrix <- as.matrix(vcf_data[ , 10:ncol(vcf_data)])

  # Extracting the data from the matrix
  sub_mat <- strsplit(vcf_matrix, ":")
  infos <- unlist(strsplit(as.character(format_infos), split = ":"))
  data_pos <- grep("GT", infos)

  if (!identical(data_pos, integer(0))) {

    vector <- sapply(sub_mat, "[", data_pos)
    data_matrix <- matrix(vector, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))
  }

  # Dealing with triallelic sites which have read depth coded as x,x
  sub_mat2 <- strsplit(data_matrix, ",")
  vector2 <- sapply(sub_mat2, "[", 1)
  data_matrix2 <- matrix(vector2, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))

  # creating a vector of heterosigousity ratio for all samples
  by_samples_metrics <- data.frame(Sample_name = colnames(vcf_data)[10:ncol(vcf_data)])
  by_samples_metrics$Proportion_heterosigousity <- colSums(data_matrix2 == "1/0" | data_matrix2 == "0/1", na.rm = FALSE) / (nrow(data_matrix2) - colSums(data_matrix2 == "./.", na.rm = FALSE))

  # Tagging rows for the presence or absence of the desired pattern in the sample name
  rows <- grep(Pattern_to_highlight_in_sample_name, by_samples_metrics$Sample_name)
  by_samples_metrics$pattern_presence <- rep("NO", nrow(by_samples_metrics))
  by_samples_metrics[rows, "pattern_presence"] <- "YES"

  # Calculatin quartile metrics that will be used in the plot
  Q3 <- summary(by_samples_metrics[["Proportion_heterosigousity"]])[5]
  Q1 <- summary(by_samples_metrics[["Proportion_heterosigousity"]])[2]
  IQR <- Q3-Q1
  treshold <- Q3+(1.5*IQR)

  library(ggplot2)
  ggplot(data = by_samples_metrics, aes(y = Proportion_heterosigousity, x = "")) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = pattern_presence), width = 0.4, shape = 1, size = 5) +
    geom_hline(aes(yintercept= treshold, linetype = "Q3 + 1.5*IQR"), colour= "blue") +
    ggtitle("Proportion heterosigousity per sample") +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    labs(x = "") +
    scale_fill_discrete(name = "New Legend Title") +
    labs(caption = paste("Input file:", vcf_name))
}

# To run the function:
# Samples_heterosigousity_boxplot(vcf_name = "FastGBS_platypus_noScaf_onlyPASS_rmIndels_minMaxAll2_mac4_max_missing02.vcf", Pattern_to_highlight_in_sample_name = "01A")



#-----------------------------------------------------------------

#                 MISSING DATA per SNP per position

#-----------------------------------------------------------------
# Plotting the proportion of missing data per SNP (or any variant) as a function of the position on the genome
Missing_data_per_position <- function(vcf_name) {

  # Loading the necessary packages
  library(data.table)
  library(ggplot2)

  # Reading the VCF file using readLines
  vcf_lines <- readLines(vcf_name)

  # Filtering out header lines (lines that start with "##")
  header_lines <- vcf_lines[grepl("^##", vcf_lines)]

  # Extracting the column names from the header line starting with "#CHROM"
  column_names <- unlist(strsplit(vcf_lines[grepl("^#CHROM", vcf_lines)], "\t"))

  # Reading the VCF data using data.table, skipping the header lines
  vcf_data <- fread(vcf_name, skip = length(header_lines) + 1, sep = "\t", header = FALSE, fill = TRUE)
  colnames(vcf_data)[1] <- "CHROM"

  # Seting the column names
  colnames(vcf_data) <- column_names
  colnames(vcf_data)[1] <- "CHROM"

  # outputing the information in the vcf file ("FORMAT" column)
  format_infos <- vcf_data[1, "FORMAT"]

  # creating a matrix
  vcf_matrix <- as.matrix(vcf_data[ , 10:ncol(vcf_data)])

  # Extracting the data from the matrix
  sub_mat <- strsplit(vcf_matrix, ":")
  infos <- unlist(strsplit(as.character(format_infos), split = ":"))
  data_pos <- grep("GT", infos)

  if (!identical(data_pos, integer(0))) {

    vector <- sapply(sub_mat, "[", data_pos)
    data_matrix <- matrix(vector, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))
  }

  # Dealing with triallelic sites which have read depth coded as x,x
  sub_mat2 <- strsplit(data_matrix, ",")
  vector2 <- sapply(sub_mat2, "[", 1)
  data_matrix2 <- matrix(vector2, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))

  # Creating a data frame to store the ratio of heterosigousity per SNP (variant)
  data_genotypes_perSNPs <- data.frame(CHROM = vcf_data$CHROM, Position = vcf_data$POS)
  data_genotypes_perSNPs$Proportion_of_missing <-  rowSums(data_matrix2 == "./.", na.rm=TRUE)/(ncol(data_matrix2))

  # Plotting the proportion of heterozigousity per SNP (variant)
  ggplot(data = data_genotypes_perSNPs, aes(x = Position, y = Proportion_of_missing)) +
    geom_point() +
    facet_wrap(~CHROM)+
    labs(y= "Proportion missing data", x = "Position (bp)") +
    ggtitle("Proportion of missing data as a function of genomic position") +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), plot.title = element_text(size = 25)) +
    labs(caption = paste("Input file:", vcf_name))
}

# To run the function:
# Missing_data_per_position(vcf_name = "FastGBS_platypus_noScaf_onlyPASS_rmIndels_minMaxAll2_mac4_max_missing02.vcf")


#-----------------------------------------------------------------

#               histogram MISSING DATA per sample

#-----------------------------------------------------------------
# Plotting an histogram of the proportion of missing data per sample
Samples_proportion_missing_histogram <- function(vcf_name) {

  # Loading the necessary packages
  library(data.table)
  library(ggplot2)

  # Reading the VCF file using readLines
  vcf_lines <- readLines(vcf_name)

  # Filtering out header lines (lines that start with "##")
  header_lines <- vcf_lines[grepl("^##", vcf_lines)]

  # Extracting the column names from the header line starting with "#CHROM"
  column_names <- unlist(strsplit(vcf_lines[grepl("^#CHROM", vcf_lines)], "\t"))

  # Reading the VCF data using data.table, skipping the header lines
  vcf_data <- fread(vcf_name, skip = length(header_lines) + 1, sep = "\t", header = FALSE, fill = TRUE)
  colnames(vcf_data)[1] <- "CHROM"

  # Seting the column names
  colnames(vcf_data) <- column_names
  colnames(vcf_data)[1] <- "CHROM"

  # outputing the information in the vcf file ("FORMAT" column)
  format_infos <- vcf_data[1, "FORMAT"]

  # creating a matrix
  vcf_matrix <- as.matrix(vcf_data[ , 10:ncol(vcf_data)])

  # Extracting the data from the matrix
  sub_mat <- strsplit(vcf_matrix, ":")
  infos <- unlist(strsplit(as.character(format_infos), split = ":"))
  data_pos <- grep("GT", infos)

  if (!identical(data_pos, integer(0))) {

    vector <- sapply(sub_mat, "[", data_pos)
    data_matrix <- matrix(vector, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))
  }

  # Dealing with triallelic sites which have read depth coded as x,x
  sub_mat2 <- strsplit(data_matrix, ",")
  vector2 <- sapply(sub_mat2, "[", 1)
  data_matrix2 <- matrix(vector2, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))

  # creating a vector of heterosigousity ratio for all samples
  by_samples_metrics <- data.frame(Sample_name = colnames(vcf_data)[10:ncol(vcf_data)])
  by_samples_metrics$Proportion_missing <- colSums(data_matrix2 == "./.", na.rm = FALSE) / nrow(data_matrix2)

  # Plotting an histogram of the proportion of heterozigousity per SNP (variant)
  ggplot(data = by_samples_metrics, aes(x = Proportion_missing)) +
    geom_histogram() +
    labs(y= "Count", x = "Proportion of missing data") +
    ggtitle("Distribution of proportion missing data per sample") +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), plot.title = element_text(size = 25)) +
    labs(caption = paste("Input file:", vcf_name))
}

# To run the function
# Samples_proportion_missing_histogram(vcf_name = "FastGBS_platypus_noScaf_onlyPASS_rmIndels_minMaxAll2_mac4_max_missing02.vcf")

#-----------------------------------------------------------------

#               Box plot MISSING DATA per sample

#-----------------------------------------------------------------
# Plotting a boxplot of the proportion of missing data per sample
# Pattern_to_highlight_in_sample_name : will color the dots on the boxplot regarding the patern has been found in the name or not
Samples_proportion_missing_boxplot <- function(vcf_name, Pattern_to_highlight_in_sample_name) {

  # Loading the necessary packages
  library(data.table)
  library(ggplot2)

  # Reading the VCF file using readLines
  vcf_lines <- readLines(vcf_name)

  # Filtering out header lines (lines that start with "##")
  header_lines <- vcf_lines[grepl("^##", vcf_lines)]

  # Extracting the column names from the header line starting with "#CHROM"
  column_names <- unlist(strsplit(vcf_lines[grepl("^#CHROM", vcf_lines)], "\t"))

  # Reading the VCF data using data.table, skipping the header lines
  vcf_data <- fread(vcf_name, skip = length(header_lines) + 1, sep = "\t", header = FALSE, fill = TRUE)
  colnames(vcf_data)[1] <- "CHROM"

  # Setting the column names
  colnames(vcf_data) <- column_names
  colnames(vcf_data)[1] <- "CHROM"

  # outputing the information in the vcf file ("FORMAT" column)
  format_infos <- vcf_data[1, "FORMAT"]

  # creating a matrix
  vcf_matrix <- as.matrix(vcf_data[ , 10:ncol(vcf_data)])

  # Extracting the data from the matrix
  sub_mat <- strsplit(vcf_matrix, ":")
  infos <- unlist(strsplit(as.character(format_infos), split = ":"))
  data_pos <- grep("GT", infos)

  if (!identical(data_pos, integer(0))) {

    vector <- sapply(sub_mat, "[", data_pos)
    data_matrix <- matrix(vector, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))
  }

  # Dealing with triallelic sites which have read depth coded as x,x
  sub_mat2 <- strsplit(data_matrix, ",")
  vector2 <- sapply(sub_mat2, "[", 1)
  data_matrix2 <- matrix(vector2, nrow = nrow(vcf_matrix), ncol = ncol(vcf_matrix))

  # creating a vector of heterosigousity ratio for all samples
  by_samples_metrics <- data.frame(Sample_name = colnames(vcf_data)[10:ncol(vcf_data)])
  by_samples_metrics$Proportion_missing <- colSums(data_matrix2 == "./.", na.rm = FALSE) / nrow(data_matrix2)

  # Tagging rows for the presence or absence of the desired pattern in the sample name
  rows <- grep(Pattern_to_highlight_in_sample_name, by_samples_metrics$Sample_name)
  by_samples_metrics$pattern_presence <- rep("NO", nrow(by_samples_metrics))
  by_samples_metrics[rows, "pattern_presence"] <- "YES"

  # Plotting the boxplot and jitter
  library(ggplot2)
  ggplot(data = by_samples_metrics, aes(y = Proportion_missing, x = "")) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = pattern_presence), width = 0.4, shape = 1, size = 5) +
    ggtitle("Proportion missing data per sample") +
    theme_classic() +
    theme(text = element_text(size = 20)) +
    labs(x = "") +
    scale_fill_discrete(name = "New Legend Title") +
    labs(caption = paste("Input file:", vcf_name))
}

# To run the function:
# Samples_proportion_missing_boxplot(vcf_name = "FastGBS_platypus_noScaf_onlyPASS_rmIndels_minMaxAll2_mac4_max_missing02.vcf", Pattern_to_highlight_in_sample_name = "01A")

