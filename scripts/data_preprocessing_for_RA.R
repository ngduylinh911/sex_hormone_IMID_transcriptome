rm(list=ls()) # clear all variables in the environment
getwd()

#set working directory
setwd("C:/Users/ngduy/Documents/Project/Final_code")

# --- Configuration Section ---
# Define the GEO Series Accession ID for the dataset
geo_accession_id <- "GSE97810"

# Define the Bioconductor annotation package name for your microarray platform (GPL ID)
annotation_package_name <- "hgu133plus2.db"
# --- End Configuration Section ---


# # 1. Install and Load necessary packages
# # Check if BiocManager is installed, if not, install it.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# # Install required Bioconductor packages
# required_bioc_packages <- c("GEOquery", "limma", "Biobase", "AnnotationDbi", annotation_package_name)
# for (pkg in required_bioc_packages) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     message(paste0("Installing Bioconductor package: ", pkg))
#     BiocManager::install(pkg)
#   }
# }

# Load the installed packages
library(GEOquery) # For downloading GEO datasets
library(dplyr) # For data manipulation
library(limma) # For linear modeling and differential expression analysis
library(Biobase) # For handling ExpressionSet objects
library(AnnotationDbi) # For annotation database interface
library(stringr) # For string manipulation
# Dynamically load the specific annotation package
library(package = annotation_package_name, character.only = TRUE)

message(paste0("Attempting to download GEO dataset: ", geo_accession_id))

# 2. Download the GEO dataset
gse <- getGEO(geo_accession_id, GSEMatrix = TRUE, AnnotGPL = FALSE)

# Handle cases where getGEO might return a list of ExpressionSet objects (e.g., multiple platforms)
if (length(gse) > 1) {
  message(paste0("Warning: GEO series ", geo_accession_id, " contains multiple platforms."))
  message("Using the first platform's ExpressionSet object. Please inspect 'gse' object if this is not desired.")
  eset <- gse[[1]] # Take the first ExpressionSet object
} else {
  eset <- gse[[1]] # Directly access the ExpressionSet object
}

message(paste0("Successfully downloaded and parsed GEO dataset: ", geo_accession_id))
message(paste0("Platform used: ", annotation(eset)))

# 3. Extract Expression Data
expression_data <- exprs(eset)

# Get the probe IDs from the row names of your expression data
probe_ids <- rownames(expression_data)

message(paste0("Raw expression matrix dimensions: ", nrow(expression_data), " probes x ", ncol(expression_data), " samples"))
message("First 5 rows and columns of raw expression data:")
print(head(expression_data[, 1:min(5, ncol(expression_data))]))

# 4. Extract and Clean Sample Metadata (Phenotype Data)
metadata <- pData(eset)

#--- function to extract new column name from cell value ---
get_name_from_cell <- function(column_vector) {
  # Find the first non-NA character string
  first_val <- na.omit(as.character(column_vector))[1]
  
  if (is.na(first_val) || !str_detect(first_val, ": ")) {
    return(NULL) # Return NULL if no valid pattern found
  }
  # Extract the part before the first ": " and trim whitespace
  new_name <- str_extract(first_val, "^[^:]+") %>% str_trim()
  return(new_name)
}
#--- End of function ---

metadata.subset <- metadata %>%
  # Select columns by their numerical index
  dplyr::select(1, 10, c(11:58), c(1407:1428))

# --- Dynamically Rename Other Columns Based on Their Cell Values ---

# Get current names of columns to consider for dynamic renaming
# Exclude 'id' and 'tissue' because they are already handled.
cols_to_dynamically_rename <- names(metadata.subset)

# Create a mapping for renaming
new_name_mapping <- list()
for (col_name in cols_to_dynamically_rename) {
  current_col_data <- metadata.subset[[col_name]]
  if (is.character(current_col_data)) { # Only process character columns
    derived_name <- get_name_from_cell(current_col_data)
    if (!is.null(derived_name) && derived_name != col_name) {
      new_name_mapping[[col_name]] <- derived_name
    }
  }
}

# Apply the dynamic renaming
colnames(metadata.subset) <- sapply(
  colnames(metadata.subset),
  function(col) {
    if (col %in% names(new_name_mapping)) {
      return(new_name_mapping[[col]])
    } else {
      return(col) # Keep original name if no mapping exists
    }
  }
)

# --- Filter and save data for female sample ---
female_metadata <- metadata.subset %>%
  dplyr::mutate(across(where(is.character), ~ stringr::str_replace(.x, ".*: ", ""))) %>%
  dplyr::filter(tissue == "PBMCs_for_RNA") %>%
  dplyr::rename(age = 51, sex = 54) %>%
  dplyr::mutate(age = as.numeric(age)) %>%  # Convert 'age' to a numeric type
  dplyr::filter(between(age, 20, 36), sex == "Female") # filter for female in the age from 20 to 36

# Ensure sample order matches between expression data and metadata
# This is crucial for correct alignment in Python
female_expression_data <- expression_data[, rownames(female_metadata)]

# Save the Processed Data to CSV Files
output_female_gene_expression_file <- "./result/RA/female_RA_expression_data.csv"
output_female_metadata_file <- "./result/RA/female_RA_metadata.csv"

write.csv(female_expression_data, output_female_gene_expression_file, row.names = TRUE)
write.csv(female_metadata, output_female_metadata_file, row.names = TRUE)

message(paste0("\n--- Data Preparation Complete ---"))
message(paste0("Gene symbol-mapped expression data saved to '", output_female_gene_expression_file, "'"))
message(paste0("Cleaned metadata saved to '", output_female_metadata_file, "'"))
#--- End of filtering female samples ---

# --- Filter and save data for male sample ---
male_metadata <- metadata.subset %>%
  dplyr::mutate(across(where(is.character), ~ stringr::str_replace(.x, ".*: ", ""))) %>%
  dplyr::filter(tissue == "PBMCs_for_RNA") %>%
  dplyr::rename(age = 51, sex = 54) %>%
  dplyr::mutate(age = as.numeric(age)) %>%  # Convert 'age' to a numeric type
  dplyr::filter(between(age, 18, 45), sex == "Male") # filter for male in the age from 20 to 36

# Ensure sample order matches between expression data and metadata
# This is crucial for correct alignment in Python
male_expression_data <- expression_data[, rownames(male_metadata)]

# Save the Processed Data to CSV Files
output_male_gene_expression_file <- "./result/RA/male_RA_expression_data.csv"
output_male_metadata_file <- "./result/RA/male_RA_metadata.csv"

write.csv(male_expression_data, output_male_gene_expression_file, row.names = TRUE)
write.csv(male_metadata, output_male_metadata_file, row.names = TRUE)

message(paste0("\n--- Data Preparation Complete ---"))
message(paste0("Gene symbol-mapped expression data saved to '", output_male_gene_expression_file, "'"))
message(paste0("Cleaned metadata saved to '", output_male_metadata_file, "'"))
#--- End of filtering male samples ---

# Optional: Clean up memory
rm(list = ls())
gc()
