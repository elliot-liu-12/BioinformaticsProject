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


# Subset the data to include only the top 5000 genes
subset_data <- log_scaled_table %>% slice_head(n = 5000)

numeric_columns <- setdiff(names(subset_data), c("Ensembl", "HUGO"))  # Get column names excluding "Ensembl" and "HUGO"



subset_data_cleaned <- subset_data[complete.cases(subset_data), ]
subset_data_cleaned <- subset_data[, colSums(is.na(subset_data)) == 0]  # Removes columns with missing values
numeric_columns <- sapply(subset_data, is.numeric)
subset_data_numeric <- subset_data[, numeric_columns]
subset_data_cleaned <- na.omit(subset_data_numeric)


#Kmeans Clustering for 5000
KM = kmeans(subset_data_cleaned,5);
autoplot(KM,subset_data_cleaned,frame=TRUE)

# Hierarchical clustering for Top 5000 Genes
dist_matrix_5000 <- dist(subset_data_cleaned, method = "euclidean")
hclust_result_5000 <- hclust(dist_matrix_5000, method = "complete")
plot(hclust_result_5000, main = "Hierarchical Clustering for Top 5000 Genes", xlab = "Genes", ylab = "Distance")

#Top 10 Genes Kmeans
subset_data <- log_scaled_table %>% slice_head(n = 10)
numeric_columns <- setdiff(names(subset_data), c("Ensembl", "HUGO"))  # Get column names excluding "Ensembl" and "HUGO"
subset_data_cleaned <- subset_data[complete.cases(subset_data), ]
subset_data_cleaned <- subset_data[, colSums(is.na(subset_data)) == 0]  # Removes columns with missing values
numeric_columns <- sapply(subset_data, is.numeric)
subset_data_numeric <- subset_data[, numeric_columns]
subset_data_cleaned <- na.omit(subset_data_numeric)
KM10 = kmeans(subset_data_cleaned,5);
autoplot(KM10,subset_data_cleaned,frame=TRUE)

# Hierarchical clustering for Top 10 Genes
dist_matrix_10 <- dist(subset_data_cleaned, method = "euclidean")
hclust_result_10 <- hclust(dist_matrix_10, method = "complete")
plot(hclust_result_10, main = "Hierarchical Clustering for Top 10 Genes", xlab = "Genes", ylab = "Distance")

#Top 100 genes Kmeans
subset_data <- log_scaled_table %>% slice_head(n = 100)
numeric_columns <- setdiff(names(subset_data), c("Ensembl", "HUGO"))  # Get column names excluding "Ensembl" and "HUGO"
subset_data_cleaned <- subset_data[complete.cases(subset_data), ]
subset_data_cleaned <- subset_data[, colSums(is.na(subset_data)) == 0]  # Removes columns with missing values
numeric_columns <- sapply(subset_data, is.numeric)
subset_data_numeric <- subset_data[, numeric_columns]
subset_data_cleaned <- na.omit(subset_data_numeric)
KM100 = kmeans(subset_data_cleaned,5);
autoplot(KM100,subset_data_cleaned,frame=TRUE)

# Hierarchical clustering for Top 100 Genes
dist_matrix_100 <- dist(subset_data_cleaned, method = "euclidean")
hclust_result_100 <- hclust(dist_matrix_100, method = "complete")
plot(hclust_result_100, main = "Hierarchical Clustering for Top 100 Genes", xlab = "Genes", ylab = "Distance")




#Top 1000 genes Kmeans
subset_data <- log_scaled_table %>% slice_head(n = 1000)
numeric_columns <- setdiff(names(subset_data), c("Ensembl", "HUGO"))  # Get column names excluding "Ensembl" and "HUGO"
subset_data_cleaned <- subset_data[complete.cases(subset_data), ]
subset_data_cleaned <- subset_data[, colSums(is.na(subset_data)) == 0]  # Removes columns with missing values
numeric_columns <- sapply(subset_data, is.numeric)
subset_data_numeric <- subset_data[, numeric_columns]
subset_data_cleaned <- na.omit(subset_data_numeric)
KM1000 = kmeans(subset_data_cleaned,5);
autoplot(KM1000,subset_data_cleaned,frame=TRUE)


# Hierarchical clustering for Top 1000 Genes
dist_matrix_1000 <- dist(subset_data_cleaned, method = "euclidean")
hclust_result_1000 <- hclust(dist_matrix_1000, method = "complete")
plot(hclust_result_1000, main = "Hierarchical Clustering for Top 1000 Genes", xlab = "Genes", ylab = "Distance")


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


# Hierarchical clustering for Top 10,000 Genes
dist_matrix_10000 <- dist(subset_data_cleaned, method = "euclidean")
hclust_result_10000 <- hclust(dist_matrix_10000, method = "complete")
plot(hclust_result_10000, main = "Hierarchical Clustering for Top 10,000 Genes", xlab = "Genes", ylab = "Distance")

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



# Plot alluvial diagram for KMeans
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


#Consensu Cluster Plus
#ConsensusClusterPlus
#Clean data before analysis
subset_data_numeric <- sweep(subset_data_numeric, 1, apply(subset_data_numeric, 1, median, na.rm = T))
#convert cleaned data into matrix
matrix <- data.matrix(subset_data_numeric)
rownames(matrix) <- subset_data$Ensembl

#precompute distance matrix to increase speed
distance <- as.dist(1-cor(matrix, method="pearson"))
#run clustering
title = "Consensus Cluster"
results_top5000  = ConsensusClusterPlus(distance, maxK = 12 ,reps = 1000, pItem = 0.8, pFeature = 1, title = title, clusterAlg = "hc", distance = "pearson",
                                        seed = 1262118388.71279, plot = "png")
#top 10 genes
top_10 <- subset_data_numeric[0:10,]

#Clean data before analysis
top_10 <- sweep(top_10, 1, apply(top_10, 1, median, na.rm = T))
#convert cleaned data into matrix
matrix <- data.matrix(top_10)
#precompute distance matrix to increase speed
distance <- as.dist(1-cor(matrix, method="pearson"))
#run clustering
title = "Top 10 Consensus Cluster"
results_top10  = ConsensusClusterPlus(matrix, maxK = 3 ,reps = 1000, pItem = 0.8, pFeature = 1, title = title, clusterAlg = "hc", distance = "pearson", 
                                      seed = 1262118388.71279, plot = "png")

#top 100 genes
top_100 <- subset_data_numeric[0:100,]

#Clean data before analysis
top_100 <- sweep(top_100, 1, apply(top_100, 1, median, na.rm = T))
#convert cleaned data into matrix
matrix <- data.matrix(top_100)
#precompute distance matrix to increase speed
distance <- as.dist(1-cor(matrix, method="pearson"))
#run clustering
title = "Top 100 Consensus Cluster"
results_top100  = ConsensusClusterPlus(matrix, maxK = 12 ,reps = 1000, pItem = 0.8, pFeature = 1, title = title, clusterAlg = "hc", distance = "pearson", 
                                       seed = 1262118388.71279, plot = "png")

#top 1000 genes
top_1000 <- subset_data_numeric[0:1000,]

#Clean data before analysis
top_1000 <- sweep(top_1000, 1, apply(top_1000, 1, median, na.rm = T))
#convert cleaned data into matrix
matrix <- data.matrix(top_1000)
#precompute distance matrix to increase speed
distance <- as.dist(1-cor(matrix, method="pearson"))
#run clustering
title = "Top 1000 Consensus Cluster"
results_top1000  = ConsensusClusterPlus(matrix, maxK = 12 ,reps = 1000, pItem = 0.8, pFeature = 1, title = title, clusterAlg = "hc", distance = "pearson", 
                                        seed = 1262118388.71279, plot = "png")
#top 10000 genes


#Clean data before analysis
top_10000 <- sweep(top_10000, 1, apply(top_10000, 1, median, na.rm = T))
#convert cleaned data into matrix
matrix <- data.matrix(top_10000)
#precompute distance matrix to increase speed
distance <- as.dist(1-cor(matrix, method="pearson"))
#run clustering
title = "Top 10000 Consensus Cluster"
results_top10000  = ConsensusClusterPlus(matrix, maxK = 12 ,reps = 1000, pItem = 0.8, pFeature = 1, title = title, clusterAlg = "hc", distance = "pearson", 
                                         seed = 1262118388.71279, plot = "png")

#Do alluvial table - varying numbers of genes
consensus_combined <- as.data.frame(cbind(results_top10[[3]][["consensusClass"]], results_top100[[12]][["consensusClass"]], results_top1000[[12]][["consensusClass"]],results_top5000[[12]][["consensusClass"]],results_top10000[[12]][["consensusClass"]]))
freqs <- data.frame(Freq = rep(1, 1021))
for(i in 1:1021)
{
  for(j in 1:5)
  {
    consensus_combined[i,j] = paste0("cluster ", consensuscombined[i,j], "", j)
  }
}
consensus_combined <- cbind(consensus_combined, freqs)


ggplot(consensus_combined,
       aes(y = Freq, axis1 = V1, axis2 = V2, axis3 = V3, axis4 = V4, axis5 = V5, fill=V1)) +
  geom_alluvium(aes(fill = V2), width = 1/8)+
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("10 genes", "100 genes", "1000 genes", "5000 genes", "10000 genes"), expand = c(.05, .05)) +
  ylab("Frequency") +
  guides(fill = guide_legend(title = "Cluster number")) +
  scale_fill_manual(values = c("cluster 1_2" = "#06b8a6", "cluster 2_2" = "#ffc4bf", "cluster 7_2" = "#80e8ff", "cluster 3_2" = "#34a853", "cluster 4_2" = "#dc493a", "cluster 5_2" = "#fff49b", "cluster 6_2" = "#d7b2ff", "cluster 8_2" = "#b7ffbe", "cluster 9_2" = "#f0bc8c", "cluster 10_2" = "#4392f1", "cluster 11_2" = "#873e23", "cluster 12_2" = "#063970")) +
  ggtitle("Consensus clustering alluvial diagram based on number of genes")


#generate heatmap
#process hierarchical cluster results
hclust_vector <- cutree(hclust_result_5000, k=9)
matrix <- log_scaled_table[1:5000,]
matrix <- matrix[-c(1,2,1024)]
matrix <- data.matrix(matrix)
par(ask = TRUE)
png("clustering-heatmap.png", width = 8, height = 6, units = "in", res = 500)
sidebars <- ComplexHeatmap::rowAnnotation(kmeans = kmeans_clusters, hierarchical = hclust_vector, consensus = results[[9]][["consensusClass"]])
topbar <- ComplexHeatmap::HeatmapAnnotation(stage = metadata$refinebio_time)
heatmap <- ComplexHeatmap::Heatmap(matrix, row_title = "Expressed genes", column_title = "Sample data", top_annotation = topbar, right_annotation = sidebars, name="Heatmap")
draw(heatmap)
dev.off()

#run tests
sample_groups <- metadata[,c(1, 21)]
temp <- as.data.frame(results_top5000[[12]][["consensusClass"]])
contingency_table <- table(temp, sample_groups)


