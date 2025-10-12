# ============================================================================
# Please note that this code is intended for reproducibility rather than instruction. 
# It allows interested users to recreate the figures and dataset used in this work.
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)  
})


pbmc <- pbmc3k
pbmc <- UpdateSeuratObject(object = pbmc)
counts_matrix <- GetAssayData(pbmc, assay = "RNA", slot = "counts")

bulk_david <- Matrix::rowSums(counts_matrix[, pbmc$patient == "David", drop = FALSE])
bulk_emily <- Matrix::rowSums(counts_matrix[, pbmc$patient == "Emily", drop = FALSE])
bulk_matrix <- cbind(bulk_david, bulk_emily)
colnames(bulk_matrix) <- c("David", "Emily")


immune_markers_clean <- list("T cell markers\n(CD3D, CD3E, CD3G)" = c("CD3D", "CD3E", "CD3G"), "CD8+ T markers\n(CD8A, CD8B)" = c("CD8A", "CD8B"),"B cell markers\n(MS4A1, CD79A, CD79B)" = c("MS4A1", "CD79A", "CD79B"))
results_h1_clean <- data.frame(patient = c("David", "Emily"))
for (marker_set in names(immune_markers_clean)) {
  genes <- immune_markers_clean[[marker_set]]
  genes_present <- genes[genes %in% rownames(bulk_matrix)]
  if (length(genes_present) > 0) {
    results_h1_clean[[marker_set]] <- colMeans(bulk_matrix[genes_present, , drop = FALSE])
  }
}


ratios_h1_clean <- results_h1_clean %>% summarise(across(-patient, ~.[patient == "Emily"] / .[patient == "David"]))
print(round(ratios_h1_clean, 2))
ratio_vals <- unlist(ratios_h1_clean)
all_balanced <- all(ratio_vals > 0.8 & ratio_vals < 1.2, na.rm = TRUE)


p1_clean <- ggplot(plot_data_h1, aes(x = patient, y = expression, fill = patient)) +
  geom_bar(stat = "identity", color = "black", size = 0.8, width = 0.6) +
  geom_text(aes(label = round(expression, 0)), vjust = 2, size = 5, fontface = "bold", color = "white") +
  facet_wrap(~marker_type, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("David" = "#4DAF4A", "Emily" = "#E41A1C")) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "lightblue"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 13, face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    panel.grid.major.x = element_blank()
  ) +
  labs(title = "Hypothesis 1: Both Patients Show Similar Marker Levels",
       subtitle = "Bulk RNA-seq cannot distinguish David from Emily",
       x = "Patient (Bulk RNA-seq sample)", 
       y = "Mean Marker Gene Expression")

composition <- pbmc@meta.data %>%
  group_by(patient, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patient) %>%
  mutate(percentage = count / sum(count) * 100)

# Calculate differences
comp_wide <- composition %>%
  select(patient, cell_type, percentage) %>%
  pivot_wider(names_from = patient, values_from = percentage, values_fill = 0) %>%
  mutate(Difference = Emily - David,
         AbsDiff = abs(Difference)) %>%
  arrange(desc(AbsDiff))

cat("Cell type composition differences:\n")
print(comp_wide)
cat("\n")

# Create the revelation plot
p3_clean <- ggplot(composition, aes(x = patient, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  geom_text(aes(label = ifelse(percentage > 3, paste0(round(percentage, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 4, fontface = "bold") +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "darkgreen"),
    panel.grid.major.x = element_blank()
  ) +
  labs(title = "Hypothesis 2: Single-Cell Reveals Different Cell Compositions!",
       subtitle = "Emily is monocyte-dominated, David is T-cell-rich",
       x = "Patient", y = "Percentage of Cells (%)", fill = "Cell Type")



# Add hypothesis labels to plots
p1_labeled <- p1_clean + 
  labs(tag = "A") +
  theme(plot.tag = element_text(size = 18, face = "bold"))

p2_labeled <- p2_clean + 
  labs(tag = "B") +
  theme(plot.tag = element_text(size = 18, face = "bold"))

p3_labeled <- p3_clean + 
  labs(tag = "C") +
  theme(plot.tag = element_text(size = 18, face = "bold"))

# Combine all three
combined_summary <- (p1_labeled / p2_labeled / p3_labeled) +
  plot_annotation(
    title = "Complete Hypothesis Testing: Bulk Fails, Single-Cell Succeeds",
    theme = theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5))
  )