#Necessary dependencies
#install.packages("BiocManager", repos = "https://cloud.r-project.org")
#BiocManager::install(c("org.Hs.eg.db"))
#BiocManager::install("AnnotationDbi")
#BiocManager::install("DESeq2")
#BiocManager::install("EnhancedVolcano", update = FALSE)
#BiocManager::install("apeglm", update = FALSE)
#BiocManager::install("ComplexHeatmap")
#install.packages("tidyverse")
#install.packages("umap")
#install.packages("clue")
#install.packages("cluster")
#install.packages("ms3")
#install.packages("magick")
#BiocManager::install("M3C")
#change working directory to correct folder - machine specific
setwd("C:/Users/donov/Downloads/Bioinformatics dataset/SRP192714")

#attach library:
library(org.Hs.eg.db)
library(magrittr)
library(AnnotationDbi)
library(tidyr)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(M3C)
library(EnhancedVolcano)
library(apeglm)
library(ComplexHeatmap)
library(cluster)
library(circlize)
library(magick)
#read files as dataframes
metadata <- readr::read_tsv("metadata_SRP192714.tsv") 
data <- readr::read_tsv("SRP192714.tsv") %>% 
  tibble::column_to_rownames("Gene")

data <- data %>%
  tibble::rownames_to_column("Gene")

#map ENSMBL IDs to HUGO ids
mappedList <- mapIds(
  org.Hs.eg.db,
  keys=data$Gene,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

#preview resulting list
#head(mappedList)

#turn list into a data frame
mapped_dataframe <- mappedList %>% 
  tibble::enframe(name="Ensembl", value="HUGO") %>%
  tidyr::unnest(cols = HUGO)
remove(mappedList)
#add HUGO to dataframe
hugo_table <- full_join(mapped_dataframe, data,
                        join_by("Ensembl" == "Gene"),
                        copy = FALSE,
                        suffix = c(".hugo", ".gene"),
                        keep = NULL,
                        na_matches = c("na", "never"),
                        relationship = NULL
)
remove(mapped_dataframe)
#make a backup of the HUGO table before log scaling
backup_table <- hugo_table
#log scale that HUGO table
log_scaled_table <- hugo_table %>% mutate(across(SRR8907879:SRR8908899, ~ log2(.x)))
#remove(hugo_table)
backup_table <- log_scaled_table
#find the median of each row and skips over non-numeric values and NA cells.
median_table <- log_scaled_table %>% rowwise() %>% mutate(row_median = median(c_across(where(is.numeric)), na.rm=TRUE))
remove(mapped_dataframe)
remove(log_scaled_table)
ggplot(median_table, aes(median_table$HUGO, median_table$row_median)) + geom_point() + geom_smooth()

#creates the column data table by extracting relevant columns from metadata
col_data <- metadata[c("refinebio_accession_code","refinebio_time")]
#confirm that all of the metadata samples are in the dataset
all(col_data$refinebio_accession_code %in% colnames(hugo_table))
#confirm that the metadata samples are in the right order
all(col_data$refinebio_accession_code == colnames(hugo_table[-c(1,2)]))

dds <- DESeqDataSetFromMatrix(countData = round(data[-1]), colData = col_data, design = ~ refinebio_time)
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup="refinebio_time")

#UMAP plot
tsne(data[,-1], labels=as.factor(col_data$refinebio_time))

#generate volcano plot
set.seed(32575)
#removes rows in the dataset that have a total read count of
#less than ten
rownames(data) <- data$Gene
filtered_data <- data[,-1] %>% 
  dplyr::filter(rowSums(.) >= 10)
#reminder: add back gene name column if necessary

#round the values
filtered_data <- round(filtered_data)
#use this data to create a new dataset
filtered_dds <- DESeqDataSetFromMatrix(countData = filtered_data, colData = col_data, design = ~ refinebio_time)
deseq_object <- DESeq(filtered_dds)
#put the results of DESeq in another object
deseq_results <- results(deseq_object)
head(deseq_results)

#process the data by adding significance column and sorting
deseq_table <- deseq_results %>%
  as.data.frame() %>% 
  tibble::rownames_to_column("Gene") %>%
  dplyr::mutate(threshold = padj < 0.001) %>%
  dplyr::arrange(dplyr::desc(log2FoldChange))

#plot 1 gene to check results
plotCounts(filtered_dds, gene="ENSG00000000003", intgroup="refinebio_time")

#save for further reference - change filepath depending on local machine
readr::write_tsv(deseq_table, file.path("C:", "Users", "donov", "Downloads", "Bioinformatics dataset", "Results", "deseqtable.tsv", fsep="\\"))

#generate volcano plot
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_table,
  lab = deseq_table$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01
)

#display it
volcano_plot

#filter out columns in the deseq table that do not meet the significance threshold
significant_genes <- deseq_table %>% filter(threshold == TRUE)
#now join main table data with it
significant_genes <- left_join(deseq_table, data, by = join_by(Gene), keep=FALSE)
#remove rows where the threshold is below that of significance
significant_genes <- significant_genes %>% filter(threshold == TRUE)
#remove unnecessary (first 8) columns and convert row names to gene names
rownames(significant_genes) <- significant_genes$HUGO
significant_genes <- significant_genes[,-c(1:8)]


rownames(significant_genes) <- NULL
colnames(significant_genes) <- NULL
matrix <- data.matrix(significant_genes)

is.numeric(matrix)

#generate heatmap
par(ask = TRUE)
png("heatmap.png", width = 8, height = 6, units = "in", res = 500)
topbar <- ComplexHeatmap::HeatmapAnnotation(stage = metadata$refinebio_time)
heatmap <- ComplexHeatmap::Heatmap(matrix, row_title = "Expressed genes", column_title = "Sample data", top_annotation = topbar, name="Heatmap")
draw(heatmap)
dev.off()


#TopGo
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("topGO")
library(topGO)
library(org.Hs.eg.db)


deseq_df <- file.path("C:", "Users", "donov", "Downloads", "Bioinformatics dataset", "Results", "deseqtable.tsv", fsep="\\")
deseq_data <- readr::read_tsv(deseq_df)

#named vector of p-values
all_genes <- setNames(deseq_data$pvalue, deseq_data$Gene)
#statistically sig genes
geneSelectionFunc <- function(pvalues) {
  return(pvalues < 0.01)
}


selected_genes <- sum(geneSelectionFunc(all_genes))
print(selected_genes)

library(org.Hs.eg.db)

# Pick a random gene from your list
keytypes(org.Hs.eg.db)
head(names(all_genes))

valid_keys <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
head(intersect(valid_keys, names(all_genes)))
length(valid_keys)

invalid_keys <- setdiff(names(all_genes), valid_keys)
head(invalid_keys)
length(invalid_keys)

all_genes_valid <- all_genes[names(all_genes) %in% valid_keys]


sample_gene <- sample(names(all_genes_valid), 1)
AnnotationDbi::select(org.Hs.eg.db, keys = sample_gene, columns = "GO", keytype = "ENSEMBL")


# Create topGOdata Object:
GOdata <- new("topGOdata",
              ontology = "BP", 
              allGenes = all_genes, 
              geneSel = geneSelectionFunc, 
              nodeSize = 10, 
              annot = annFUN.org,
              mapping = "org.Hs.eg.db",
              ID = "ENSEMBL")

GOdata

GOresults <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
significant_terms <- GenTable(GOdata, 
                              classicFisher = GOresults, 
                              topNodes = length(GOresults@score), 
                              numChar = 80)
head(significant_terms)
top_terms <- head(sort(GOresults@score, decreasing = FALSE), 10)
barplot(-log10(top_terms), las = 2, names.arg = names(top_terms), main = "TopGO Enrichment Results", ylab = "-log10(p-value)", cex.names = 0.7)



#clustProfiler
#BiocManager::install("clusterProfiler")
#BiocManager::install("DOSE")

# Load necessary libraries
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

# Read the DESeq2 results
deseq_table <- read.delim("C:/Users/donov/Downloads/Bioinformatics dataset/Results/deseqtable.tsv", header=TRUE, sep="\t")

# Filter DEGs based on criteria (p-value < 0.05 and abs(log2FoldChange) > 1)
DEGs <- as.character(subset(deseq_table, padj < 0.05 & abs(log2FoldChange) > 1)$Gene)

# Convert Ensembl gene IDs to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db, 
                     keys = DEGs, 
                     column = "ENTREZID", 
                     keytype = "ENSEMBL", 
                     multiVals = "first")

# Perform GO enrichment using the Entrez IDs
enrichGO_result <- enrichGO(entrez_ids, 
                            OrgDb='org.Hs.eg.db', 
                            keyType="ENTREZID", 
                            pAdjustMethod="BH", 
                            pvalueCutoff=0.05, 
                            qvalueCutoff=0.2)

# Print out the top GO enrichment results
head(enrichGO_result)

# Visualize the GO enrichment results
dotplot(enrichGO_result, showCategory=10)


#gProfiler2
#install.packages("gprofiler2")
library(gprofiler2)

deseq_table <- read.delim("C:/Users/donov/Downloads/Bioinformatics dataset/Results/deseqtable.tsv", header=TRUE, sep="\t")
DEGs <- as.character(subset(deseq_table, padj < 0.05 & abs(log2FoldChange) > 1)$Gene)

#Unlike clusterProfiler, gprofiler2 does not need the conversion of gene IDs to Entrez format. It can work with Ensembl IDs directly.
gp_results <- gost(DEGs, organism="hsapiens", sources=c("GO:BP"))
concise_results <- gp_results[["result"]][, c("term_name", "p_value", "term_size", "intersection_size")]


# Basic bar plotting of the significant terms based on their p-values

# Extract the top 10 results based on p-value
top_results <- head(potential_results[order(potential_results$p_value), ], 10)
print(top_results)

# Plot the top results using a bar plot without y-axis labels
bar_midpoints <- barplot(-log10(top_results$p_value), 
                         las=2, 
                         main="Top GO terms by -log10(p-value)", 
                         col="lightblue", 
                         horiz=TRUE, 
                         border="white")

# Place the term names on the bars using the text function
text(x=-log10(top_results$p_value)/2, y=bar_midpoints, labels=top_results$term_name, col="black", pos=4, cex=0.7)


