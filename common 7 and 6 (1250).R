setwd("E:/BIOTECH 19_23/Graduation  project/data seys/FINAL Insaha'allah/Discovery/GSE106817 (dataset7)")
ds7 = read.csv('Top Table(1000)_GSE7trend.csv')
setwd("E:/BIOTECH 19_23/Graduation  project/data seys/FINAL Insaha'allah/Discovery/GSE211692 (dataset6)/script 1")
ds6 = read.csv('Top Table(1000)_GSE6trend.csv')
setwd("E:/BIOTECH 19_23/Graduation  project/data seys/FINAL Insaha'allah/Discovery/Common")

common_microRNAs <- intersect(ds7$miRNA_ID_LIST, ds6$X)
results <- data.frame(MicroRNA = common_microRNAs)

# Retrieve adjusted p-value, logFC, and label for each dataset
results$Adj_P_Val_GSE106817 <- ds7$adj.P.Val[match(common_microRNAs, ds7$miRNA_ID_LIST)]
results$LogFC_GSE106817 <- ds7$logFC[match(common_microRNAs, ds7$miRNA_ID_LIST)]
results$label_GSE106817 <- ds7$label[match(common_microRNAs, ds7$miRNA_ID_LIST)]

results$Adj_P_Val_GSE211692 <- ds6$adj.P.Val[match(common_microRNAs, ds6$X)]
results$LogFC_GSE211692 <- ds6$logFC[match(common_microRNAs, ds6$X)]
results$label_GSE211692 <- ds6$label[match(common_microRNAs, ds6$X)]

ourmir1 <- results[grepl("hsa-miR-155", results$MicroRNA), ]
ourmir2 <- results[grepl("hsa-miR-335", results$MicroRNA), ]
ourmir3 <- results[grepl("hsa-miR-373", results$MicroRNA), ]
ourmir4 <- results[grepl("hsa-miR-27a", results$MicroRNA), ]
ourmir5 <- results[grepl("hsa-miR-181", results$MicroRNA), ]
ourmir6 <- results[grepl("hsa-miR-146", results$MicroRNA), ]
ourmir7 <- results[grepl("hsa-miR-21", results$MicroRNA), ]
ourmir8 <- results[grepl("hsa-miR-222", results$MicroRNA), ]
ourmir9 <- results[grepl("hsa-miR-221", results$MicroRNA), ]
combined_results <- rbind(ourmir1, ourmir2, ourmir3, ourmir4, ourmir5, ourmir6, ourmir7, ourmir8, ourmir9)
# Print the result table
write.csv(results, file = 'common1000trend.csv')
write.csv(combined_results, file = 'ourmir1250trend.csv')




























