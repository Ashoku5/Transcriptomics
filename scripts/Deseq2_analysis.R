if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
if (!requireNamespace("optparse", quietly=TRUE)) install.packages("optparse")
if (!requireNamespace("DESeq2", quietly=TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("pheatmap", quietly=TRUE)) install.packages("pheatmap")
if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(glue)
library(optparse)

counts_matrix <- "for_deseq2.txt"
group1 <- "control" #control
group2 <- "treated" #case
log2FC_threshold <- 1
pvalue_threshold <- 0.05
padj_threshold <- 0.5
folder_path <- "OUTPUT"

#preprocess
dir.create(folder_path, showWarnings = FALSE)
data_raw <- read.table(counts_matrix, header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
group_labels <- as.character(data_raw[1,])
data <- data_raw[-1,]
data[is.na(data)] <- 0
data <- as.data.frame(lapply(data, function(x) as.numeric(as.character(x))), row.names = rownames(data), check.names = FALSE)
sample_info <- data.frame(row.names = colnames(data), group = group_labels, check.names = FALSE)

#make replicates for solo groups
single_samples <- which(table(sample_info$group) == 1)
for (group in names(single_samples)) {
  sample_name <- rownames(sample_info)[sample_info$group == group]
  new_sample_name <- paste0(sample_name, "_dup1")
  data[[new_sample_name]] <- data[[sample_name]] + 3
  sample_info <- rbind(sample_info, data.frame(row.names = new_sample_name, group = group))
}

groups <- sample_info[sample_info$group %in% c(group1, group2), ]
samples <- rownames(sample_info)[sample_info$group %in% c(group1, group2)]
data_filtered <- data[, samples]

data_filtered[] <- lapply(data_filtered, function(x) {
  x[is.na(x)] <- 0
  as.integer(round(x))
})

#DESEQ2
dds <- DESeqDataSetFromMatrix(countData = data_filtered,
                              colData = data.frame(condition = groups),
                              design = ~condition)
dds$condition <- factor(dds$condition, levels = c(group1, group2))

#dds <- DESeq(dds, fitType="mean")
dds <- DESeq(dds)
all_results <- results(dds)
write.csv(all_results, file = "OUTPUT/DEGs_unfiltered.csv", row.names = T)
results <- results(dds, alpha = padj_threshold)
summary(results)

normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file = "OUTPUT/normalized_counts.csv", quote = TRUE)

deg <- all_results[!is.na(all_results$pvalue) & all_results$pvalue <= pvalue_threshold,]
deg <- deg[!is.na(deg$padj) & deg$padj <= padj_threshold,]
deg <- deg[order(deg$log2FoldChange, decreasing = TRUE), ]
deg <- deg[abs(deg$log2FoldChange) > log2FC_threshold | is.na(deg$log2FoldChange), ]

deg$Gene <- rownames(deg)
deg <- deg[,c("Gene", names(deg)[1:6])]
positive_results <- deg[deg$log2FoldChange > 0,]
negative_results <- deg[deg$log2FoldChange < 0,]

positive_results <- positive_results[order(-positive_results$log2FoldChange),]
negative_results <- negative_results[order(-negative_results$log2FoldChange),]

sorted_results <- rbind(positive_results, negative_results)
sorted_results_filtered <- sorted_results[abs(sorted_results$log2FoldChange) > 1, ]
num_upregulated <- sum(sorted_results_filtered$log2FoldChange > 1)
num_downregulated <- sum(sorted_results_filtered$log2FoldChange < -1)

#Visualization
deg_counts <- data.frame(
  Regulation = factor(c("UP", "DOWN"), levels = unique(c("UP", "DOWN"))),
  DEGs = c(num_upregulated, num_downregulated)
)

deg_counts_bar <- ggplot(deg_counts, aes(x = Regulation, y = DEGs, fill = Regulation)) +
                          geom_bar(stat = "identity") +
                          labs(title = "Number of Differentially Expressed Genes (DEGs)")
                          #scale_fill_manual(values = c("#1f77b4", "#ff220eff"))
ggsave(filename = "OUTPUT/DEG_counts.png", plot = deg_counts_bar, dpi = 600)

sorted_results_filtered$upregulated <- num_upregulated
sorted_results_filtered$downregulated <- num_downregulated
write.csv(sorted_results_filtered, file = "OUTPUT/DEGs.csv", row.names = FALSE)

results$Regulation <- ifelse(results$padj < padj_threshold & results$pvalue < pvalue_threshold,
                             ifelse(results$log2FoldChange > log2FC_threshold, "Upregulated",
                                    ifelse(results$log2FoldChange < -log2FC_threshold, "Downregulated", "Not Significant")),
                             "Not Significant")

res <- as.data.frame(results)
res_df <- tibble::rownames_to_column(res, var = "GeneName")
comparison <- paste(group1, "-", group2)
volcano_plot19 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Regulation, label = GeneName)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "gray") +
  labs(title = comparison,
       x = "log2 Fold Change",
       y = "-log10(p.adj)") +
  guides(color = guide_legend(title = "Regulation")) +
  theme(plot.title = element_text(size = 15),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
ggsave(filename = "OUTPUT/volcano_plot.png", plot = volcano_plot19, width = 12, height = 7, units = "in", dpi = 600)

normalized <- read.csv("OUTPUT/normalized_counts.csv", row.names = 1, header = TRUE, check.names = FALSE )
if (nrow(sorted_results_filtered) >= 50) {
    top_25_row_names <- row.names(sorted_results_filtered)[1:25]
    bottom_25_row_names <- row.names(sorted_results_filtered)[(nrow(sorted_results_filtered) - 24):nrow(sorted_results_filtered)]
    combined_rows <- rbind(normalized[top_25_row_names, , drop = FALSE],
                        normalized[bottom_25_row_names, , drop = FALSE])
} else {
    row_names <- row.names(sorted_results_filtered)
    combined_rows <- rbind(normalized[row_names, ])
}
combined_rows <- combined_rows[, !grepl("_dup1$", colnames(combined_rows))]
if (nrow(combined_rows) == 0) stop("No valid rows found in normalized data.")

heatmap <- pheatmap(combined_rows,
                    cluster_rows = TRUE,
                    cluster_cols = FALSE,
                    #color = colorRampPalette(c("green", "black", "red"))(50),
                    border_color = "black",
                    scale = "row",
                    main = "",
                    fontsize = 28,
                    fontsize_row = 28,
                    fontsize_col = 28)
ggsave(filename = "OUTPUT/Heatmap.png", plot = heatmap , width = 20, height = 25, dpi = 300)

files2zip <- dir('OUTPUT', full.names = TRUE)
#zip(zipfile = glue('{group1}-{group2}.zip'), files = files2zip)
#unlink("OUTPUT", recursive = TRUE)