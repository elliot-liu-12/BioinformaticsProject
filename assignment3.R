#Install packages
install.packages("stats")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("ggfortify")
install.packages("matrixStats")
install.packages("ggalluvial")

library(cluster)
library(matrixStats)
library(ggplot2)
library(stats)
library(dplyr)
library(ggplot2)
library(ggfortify)

# Calculate the variance of each gene and select the top 5,000 genes
#numeric_columns <- sapply(log_scaled_table, is.numeric)
#variance_values <- apply(log_scaled_table[, numeric_columns], 2, var)  # Calculate variance for each gene
#top_genes <- names(sort(variance_values, decreasing = TRUE)[1:1021])  # Select top 1021 gene names

# Calculate the variance of each gene and select the top 5,000 genes
log_scaled_matrix <- data.matrix(log_scaled_table) # Calculate variance for each gene
log_scaled_table$row_var <- rowVars(log_scaled_matrix)
#sorts dataframe by variance
log_scaled_table <- log_scaled_table %>% arrange(desc(row_var))


# Subset the data to include only the top 1021 genes
subset_data <- log_scaled_table %>% slice_head(n = 5000)

numeric_columns <- setdiff(names(subset_data), c("Ensembl", "HUGO"))  # Get column names excluding "Ensembl" and "HUGO"



subset_data_cleaned <- subset_data[complete.cases(subset_data), ]
subset_data_cleaned <- subset_data[, colSums(is.na(subset_data)) == 0]  # Removes columns with missing values
numeric_columns <- sapply(subset_data, is.numeric)
subset_data_numeric <- subset_data[, numeric_columns]
subset_data_cleaned <- na.omit(subset_data_numeric)


#Kmeans Clustering
KM = kmeans(subset_data_cleaned,5);
autoplot(KM,subset_data_cleaned,frame=TRUE)

#Top 10 Genes
subset_data <- log_scaled_table %>% slice_head(n = 10)
numeric_columns <- setdiff(names(subset_data), c("Ensembl", "HUGO"))  # Get column names excluding "Ensembl" and "HUGO"
subset_data_cleaned <- subset_data[complete.cases(subset_data), ]
subset_data_cleaned <- subset_data[, colSums(is.na(subset_data)) == 0]  # Removes columns with missing values
numeric_columns <- sapply(subset_data, is.numeric)
subset_data_numeric <- subset_data[, numeric_columns]
subset_data_cleaned <- na.omit(subset_data_numeric)
KM10 = kmeans(subset_data_cleaned,5);
autoplot(KM10,subset_data_cleaned,frame=TRUE)

#Top 100 genes
subset_data <- log_scaled_table %>% slice_head(n = 100)
numeric_columns <- setdiff(names(subset_data), c("Ensembl", "HUGO"))  # Get column names excluding "Ensembl" and "HUGO"
subset_data_cleaned <- subset_data[complete.cases(subset_data), ]
subset_data_cleaned <- subset_data[, colSums(is.na(subset_data)) == 0]  # Removes columns with missing values
numeric_columns <- sapply(subset_data, is.numeric)
subset_data_numeric <- subset_data[, numeric_columns]
subset_data_cleaned <- na.omit(subset_data_numeric)
KM100 = kmeans(subset_data_cleaned,5);
autoplot(KM100,subset_data_cleaned,frame=TRUE)

#Top 1000 genes
subset_data <- log_scaled_table %>% slice_head(n = 1000)
numeric_columns <- setdiff(names(subset_data), c("Ensembl", "HUGO"))  # Get column names excluding "Ensembl" and "HUGO"
subset_data_cleaned <- subset_data[complete.cases(subset_data), ]
subset_data_cleaned <- subset_data[, colSums(is.na(subset_data)) == 0]  # Removes columns with missing values
numeric_columns <- sapply(subset_data, is.numeric)
subset_data_numeric <- subset_data[, numeric_columns]
subset_data_cleaned <- na.omit(subset_data_numeric)
KM1000 = kmeans(subset_data_cleaned,5);
autoplot(KM1000,subset_data_cleaned,frame=TRUE)

#Top 10000 genes
subset_data <- log_scaled_table %>% slice_head(n = 10000)
numeric_columns <- setdiff(names(subset_data), c("Ensembl", "HUGO"))  # Get column names excluding "Ensembl" and "HUGO"
subset_data_cleaned <- subset_data[complete.cases(subset_data), ]
subset_data_cleaned <- subset_data[, colSums(is.na(subset_data)) == 0]  # Removes columns with missing values
numeric_columns <- sapply(subset_data, is.numeric)
subset_data_numeric <- subset_data[, numeric_columns]
subset_data_cleaned <- na.omit(subset_data_numeric)
KM10000 = kmeans(subset_data_cleaned,5);
autoplot(KM10000,subset_data_cleaned,frame=TRUE)

#Alluvial diagram

library(ggalluvial)

KM_alluvial_df <- as.data.frame(KM10$cluster)

#KM_alluvial_df <- KM_alluvial_df[1:43363, ]
#metadata <- metadata[1:43363, ]
KM_alluvial_df <- cbind(KM_alluvial_df,
                        as.data.frame(KM100$cluster),
                        as.data.frame(KM1000$cluster),
                        metadata$refinebio_time)
colnames(KM_alluvial_df) <- c("ten_genes", "hundred_genes", "thousand_genes", "stage")

KM_alluvial_df <- KM_alluvial_df %>%
  group_by(ten_genes, hundred_genes, thousand_genes, stage) %>%
  summarise(freq = n())



# Plot alluvial diagram
ggplot(KM_alluvial_df,
       aes(y = freq, axis1 = ten_genes, axis2 = hundred_genes, axis3 = thousand_genes)) +
  geom_alluvium(aes(fill = stage), width = 1/12) +
  geom_stratum(width = 1/8, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("10 genes", "100 genes", "1000 genes"), expand = c(.02, .02)) +
  ylab("Frequency") +
  guides(fill = guide_legend(title = "Stage")) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("KM clustering alluvial diagram based on number of genes")

#Hierarchical
dist_matrix <- dist(subset_data_cleaned, method = "euclidean")
hclust_result <- hclust(dist_matrix, method = "complete")
autoplot(hclust_result,frame=TRUE)
