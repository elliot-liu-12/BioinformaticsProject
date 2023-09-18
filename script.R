#Necessary dependencies
#install.packages("BiocManager", repos = "https://cloud.r-project.org")
#BiocManager::install(c("org.Hs.eg.db"))
#BiocManager::install("AnnotationDbi")
#install.packages("tidyverse")


#change working directory to correct folder - machine specific
#setwd("C:/Users/ellio/Downloads/Bioinformatics dataset/SRP192714")

#attach library:
library(org.Hs.eg.db)
library(magrittr)
library(AnnotationDbi)
library(tidyr)
library(dplyr)


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
head(mappedList)

#turn list into a data frame
mapped_dataframe <- mappedList %>% 
  tibble::enframe(name="Ensembl", value="HUGO") %>%
  tidyr::unnest(cols = HUGO)
#add HUGO to dataframe
hugo_table <- full_join(mapped_dataframe, data,
                        join_by("Ensembl" == "Gene"),
                        copy = FALSE,
                        suffix = c(".hugo", ".gene"),
                        keep = NULL,
                        na_matches = c("na", "never"),
                        relationship = NULL
                        )
#make a backup of the HUGO table before log scaling
backup_table <- hugo_table
#log scale that HUGO table
log_scaled_table <- hugo_table %>% mutate(across(SRR8907879:SRR8908899, ~ log2(.x)))
