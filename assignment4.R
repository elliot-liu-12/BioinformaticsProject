# Installing CRAN packages (if not already installed)
#install.packages(c("tidymodels", "e1071", "randomForest", "kernlab", "recipes","ranger"))

# Load necessary libraries
library(tidymodels)
library(recipes)   # Add this line
library(e1071)  # For SVM
library(randomForest)  # For Random Forest
library(kernlab)
library(ranger) # For Random Forest

setwd("C:/Users/donov/Downloads/Bioinformatics dataset/SRP192714")

#getting needed data from assignment 2
#read files as dataframes
metadata <- readr::read_tsv("metadata_SRP192714.tsv") 
data <- readr::read_tsv("SRP192714.tsv") %>% 
  tibble::column_to_rownames("Gene")

data <- data %>%
  tibble::rownames_to_column("Gene")

# Filter the metadata to include only 'late acute' and 'early acute' samples
filtered_metadata <- metadata %>%
  filter(refinebio_time %in% c('late acute', 'early acute'))

# Join with the expression data to get the corresponding samples
filtered_data <- data %>%
  select(Gene, all_of(filtered_metadata$refinebio_accession_code))

rm(metadata, data)  # Remove these objects after filtering
gc()

# Create a vector of gene subset sizes
gene_subsets <- c(5000, 10, 100, 1000, 10000)


for (genes_n in gene_subsets) {
  # 1. Subset the data to the top most variable genes
  var_genes <- filtered_data %>%
    rowwise() %>%
    mutate(variance = var(c_across(-Gene))) %>%
    arrange(-variance) %>%
    head(genes_n)
  
  # Pivot var_genes to long format
  long_var_genes <- var_genes %>%
    tidyr::pivot_longer(
      cols = -Gene,
      names_to = "refinebio_accession_code",
      values_to = "expression_value"
    )
  
  # Join with filtered_metadata
  var_genes_with_labels <- left_join(
    long_var_genes, 
    filtered_metadata, 
    by = "refinebio_accession_code"
  )
  
  rm(long_var_genes,var_genes)  # Remove this object after joining
  gc()
  
  # Splitting the data into train and test sets
  set.seed(123)  # for reproducibility
  train_indices <- sample(1:nrow(var_genes_with_labels), size = 0.8 * nrow(var_genes_with_labels))
  train_data_vg <- var_genes_with_labels[train_indices, ]
  test_data <- var_genes_with_labels[-train_indices, ]
  
  rm(var_genes_with_labels)  # If you don't need the complete dataset anymore.
  gc()
  
  #Check for duplicates:
  train_data_vg <- train_data_vg[!duplicated(train_data_vg), ]
  
  #Drop columns that have only NA values:
  train_data_vg <- train_data_vg[, colSums(is.na(train_data_vg)) != nrow(train_data_vg)]
  
  #Convert character columns with few unique values to factors (let's say, columns with less than 100 unique values for this example):
  char_cols <- sapply(train_data_vg, is.character)
  for(col in names(train_data_vg)[char_cols]) {
    if(length(unique(train_data_vg[[col]])) < 100) {
      train_data_vg[[col]] <- as.factor(train_data_vg[[col]])
    }
  }
  
  #Handle missing values. For simplicity, let's remove rows with NA. Depending on your analysis, imputation might be a better choice.
  train_data_vg <- na.omit(train_data_vg)
  
  # 2. Supervised Analysis
  
  #a.) SVM (Support Vector Machine)
  
  # Model specification
  svm_spec <- svm_poly(degree = 1, cost = 1) %>%
    set_engine("kernlab") %>%
    set_mode("classification")
  
  # Recipe
  # Define a vector of potentially problematic columns
  columns_to_remove <- c("Gene", "refinebio_cell_line", "refinebio_compound", "refinebio_disease", 
                         "refinebio_disease_stage", "refinebio_genetic_information", 
                         "refinebio_race", "refinebio_source_archive_url", 
                         "refinebio_specimen_part", "refinebio_treatment")
  
  # Create an updated recipe with only the columns that are present in train_data_vg
  existing_columns_in_train <- columns_to_remove[columns_to_remove %in% names(train_data_vg)]
  
  single_level_factors <- sapply(train_data_vg, function(col) {
    is.factor(col) && length(levels(col)) == 1
  })
  
  train_data_vg <- train_data_vg[, !single_level_factors]
  
  svm_recipe <- recipe(refinebio_time ~ ., data = train_data_vg) %>%
    step_rm(all_of(existing_columns_in_train)) %>%
    step_nzv(all_predictors())
  
  svm_wf <- workflow() %>%
    add_recipe(svm_recipe) %>%
    add_model(svm_spec)
  
  single_level_factors_updated <- sapply(train_data_vg, function(col) {
    is.factor(col) && length(levels(col)) == 1
  })
  
  # Fit the updated workflow
  set.seed(123) # for reproducibility
  train_data_vg_subset <- train_data_vg[sample(1:nrow(train_data_vg), size = genes_n), ]
  svm_fit <- svm_wf %>% fit(data = train_data_vg_subset)
  
  # Train SVM
  svm_res <- svm_wf %>%
    fit(data = train_data_vg_subset)
  
  # Predict and calculate accuracy
  svm_preds <- predict(svm_res, new_data = test_data)
  
  trimmed_test_data <- test_data[1:length(svm_preds), ]
  
  results <- data.frame(
    Prediction = svm_preds,
    Truth = trimmed_test_data$refinebio_time
  )
  
  # Convert 'Truth' and '.pred_class' to factor
  results$Truth <- factor(results$Truth, levels = c("early acute", "late acute"))
  results$.pred_class <- as.factor(results$.pred_class)
  
  # Compute accuracy
  acc_result <- accuracy(results, Truth, .pred_class)
  
  # Print the accuracy for SVM
  cat("For", genes_n, "genes, SVM Accuracy:", toString(acc_result), "\n")
  
  rm(svm_fit, svm_res, svm_preds)
  gc()
  
  #b.) Logistic Regression
  
  # Model specification
  log_reg_spec <- logistic_reg() %>%
    set_engine("glm") %>%
    set_mode("classification")
  
  set.seed(123)  # Setting a seed for reproducibility
  
  # Separate out early acute and late acute data
  early_acute_data <- train_data_vg_subset[train_data_vg_subset$refinebio_time == "early acute", ]
  late_acute_data  <- train_data_vg_subset[train_data_vg_subset$refinebio_time == "late acute", ]
  
  # Sample from the "late acute" category to match the size of the "early acute" category
  late_acute_sampled <- late_acute_data[sample(1:nrow(late_acute_data), nrow(early_acute_data)), ]
  
  # Combine the downsampled "late acute" data with the "early acute" data
  balanced_data <- rbind(early_acute_data, late_acute_sampled)
  
  # Check the balance
  table(balanced_data$refinebio_time)
  
  # Update the recipe to exclude the non-numeric predictors
  svm_recipe_standardized <- recipe(refinebio_time ~ ., data = balanced_data) %>%
    step_rm(Gene, refinebio_accession_code, refinebio_processed, refinebio_title, refinebio_processor_id) %>%
    step_center(all_predictors(), -all_nominal()) %>%
    step_scale(all_predictors(), -all_nominal())
  
  # Workflow
  log_reg_wf <- workflow() %>%
    add_recipe(svm_recipe_standardized) %>%
    add_model(log_reg_spec)
  
  # Train logistic regression on balanced_data
  log_reg_res <- log_reg_wf %>%
    fit(data = balanced_data)
  
  # Predict on test data
  log_reg_preds <- predict(log_reg_res, new_data = test_data) %>%
    bind_cols(test_data)  # to combine predictions with true values
  
  # Convert the 'refinebio_time' and '.pred_class' to factor type
  log_reg_preds$refinebio_time <- as.factor(log_reg_preds$refinebio_time)
  log_reg_preds$.pred_class <- as.factor(log_reg_preds$.pred_class)
  
  # Compute accuracy
  acc_result <- accuracy(log_reg_preds, refinebio_time, .pred_class)
  
  # Print the accuracy for Logistic Regression
  cat("For", genes_n, "genes, Logistic Regression Accuracy:", toString(acc_result), "\n")
  
  rm(log_reg_res, log_reg_preds)
  gc()
  
  # c.) Random Forest
  
  # Convert problematic columns to ordered factors
  balanced_data$refinebio_accession_code <- factor(balanced_data$refinebio_accession_code, ordered = TRUE)
  balanced_data$refinebio_title <- factor(balanced_data$refinebio_title, ordered = TRUE)
  balanced_data$Gene <- factor(balanced_data$Gene, ordered = TRUE)
  
  # Simple random forest specification
  rf_spec <- rand_forest(trees = 1800) %>%
    set_engine("ranger", respect.unordered.factors = "partition") %>%
    set_mode("classification")
  
  # Updated recipe to remove additional problematic columns
  simple_recipe <- recipe(refinebio_time ~ ., data = balanced_data) %>%
    step_rm(refinebio_accession_code, refinebio_title, Gene)
  
  # Workflow
  rf_wf <- workflow() %>%
    add_recipe(simple_recipe) %>%
    add_model(rf_spec)
  
  # Training on the small subset
  rf_res <- rf_wf %>%
    fit(data = balanced_data)
  
  # Prediction (using the same subset for simplicity)
  rf_preds <- predict(rf_res, new_data = balanced_data)
  
  correct_predictions <- sum(rf_preds$.pred_class == balanced_data$refinebio_time)
  total_predictions <- nrow(balanced_data)
  accuracy <- correct_predictions / total_predictions
  
  # Print the accuracy for Random Forest
  cat("For", genes_n, "genes, Random Forest Accuracy:", accuracy, "\n")
  
  rm(rf_res, rf_preds)
  gc()
}

# Load necessary libraries
library(readr)
library(dplyr)
library(tibble)
library(pheatmap)

# Load data (adjust paths as needed)
metadata <- read_tsv("metadata_SRP192714.tsv")
data <- read_tsv("SRP192714.tsv") %>%
  column_to_rownames("Gene")

# Filter metadata and data for 'late acute' and 'early acute' samples
filtered_metadata <- metadata %>%
  filter(refinebio_time %in% c('late acute', 'early acute'))
filtered_data <- data %>%
  select(all_of(filtered_metadata$refinebio_accession_code))

# Identify top 1000 most variable genes
var_genes <- filtered_data %>%
  rowwise() %>%
  mutate(variance = var(c_across(everything()))) %>%
  arrange(-variance) %>%
  head(1000)

# Prepare annotation data
annotation_data <- filtered_metadata %>%
  select(refinebio_accession_code, refinebio_time) %>%
  column_to_rownames("refinebio_accession_code")

# Create a heatmap
pheatmap(mat = as.matrix(var_genes), 
         annotation_col = annotation_data, 
         main = "Heatmap of Top 1000 Most Variable Genes",
         xlab = "Samples",
         ylab = "Genes",
         show_colnames = FALSE,
         annotation_colors = list(refinebio_time = c('late acute' = 'blue', 'early acute' = 'red')),
         annotation_legend = TRUE)