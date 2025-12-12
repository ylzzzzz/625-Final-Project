mydata <- read.csv("~/Documents/DNA/2025Fall/Biostat625/Final/Alzheimers_Disease_cleandata_final.csv",row.names = 1)
gene_data <- read.csv("~/Documents/DNA/2025Fall/Biostat625/Final/Batch_corrected_expression.csv",row.names = 1)
clinical_data <- mydata[, 1:6]
batch_labels <- clinical_data$batch_id
disease_labels <- clinical_data$Alzheimers_Disease
length(disease_labels)

library(limma)
design <- model.matrix(~ Alzheimers_Disease, data = clinical_data)
fit <- lmFit(t(gene_data), design)
fit <- eBayes(fit)
res <- topTable(fit,
                coef = "Alzheimers_Disease",
                adjust = "BH",
                number = Inf)


################################################################################
library(EnhancedVolcano)
library(cowplot) # Install if needed: install.packages("cowplot")

# 1. Filter: Keep only cohorts with significant genes
# We define significance as adj.P.Val < 0.05 and abs(logFC) > 1

# Check if there are ANY genes meeting the criteria
sig_genes = subset(res, adj.P.Val < 0.05 & abs(logFC) > 1)



# Create the plot object but DO NOT print it yet. Save it to a list.
p = EnhancedVolcano(res,
                    lab = rownames(res),
                    x = "logFC",
                    y = "adj.P.Val",
                    title = "Merged data before batch effect correction",
                    subtitle = paste(nrow(sig_genes), "sig genes"),
                    pCutoff = 0.05,
                    FCcutoff = 1,
                    caption = "")



final_grid = plot_grid(plotlist = p, ncol = , labels = "AUTO")
print(final_grid)









