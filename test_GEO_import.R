install.packages("BiocManager")
BiocManager::install("GEOquery")
library(tidyverse)


library(GEOquery)
gpl570_data = getGEO("GPL570")

gpl570_dataframe = Table(gpl570_data)

gpl570_expanded_dataframe = gpl570_dataframe %>%
  # Apply separate_rows to each column sequentially
  separate_rows(`Gene Symbol`, sep = " /// ") %>%
  separate_rows(`ENTREZ_GENE_ID`, sep = " /// ") %>%
  separate_rows(`RefSeq Transcript ID`, sep = " /// ") %>%
  # Optionally, you can now remove leading/trailing whitespace
  mutate_all(trimws)

gpl570_wide_dataframe = gpl570_dataframe %>%
  # This creates GeneSymbol_1, GeneSymbol_2
  separate(`Gene Symbol`, into = paste0("GeneSymbol_", 1:2), sep = " /// ", fill = "right") %>%
  # This creates EntrezID_1, EntrezID_2
  separate(`ENTREZ_GENE_ID`, into = paste0("EntrezID_", 1:2), sep = " /// ", fill = "right") %>%
  # This creates many columns for RefSeq IDs (adjust '10' if you need more)
  separate(`RefSeq Transcript ID`, into = paste0("RefSeqID_", 1:10), sep = " /// ", fill = "right")
