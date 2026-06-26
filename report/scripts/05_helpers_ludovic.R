library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyverse)
library(Biostrings)

# ============================================================================
# Sourced ONCE from run_msstats.Rmd before the comparison loop.
# Builds the annotation table, defines helper functions, and initializes
# the accumulator used by the per-comparison child template.
# ============================================================================


# ----------------------------------------------------------------------------
# 1. Build protein annotation table ONCE from FASTA
# ----------------------------------------------------------------------------
message("Building protein annotation table from FASTA...")

.fasta <- readAAStringSet(fasta_file)
.headers <- names(.fasta)

annotation_table <- tibble(header = .headers) %>%
  mutate(
    db           = str_extract(header, "^(sp|tr)"),
    Protein_ID   = str_extract(header, "(?<=\\|)[A-Z0-9]+(?=\\|)"),
    Protein_entry = str_extract(header, "(?<=\\|)[A-Z0-9_]+(?= )"),
    Protein_name = str_extract(header, "(?<= )[^|]+?(?= OS=)"),
    Organism     = str_extract(header, "(?<=OS=).+?(?= OX=)"),
    TaxID        = str_extract(header, "(?<=OX=)\\d+"),
    Gene_name    = str_extract(header, "(?<=GN=)\\S+"),
    PE_level     = str_extract(header, "(?<=PE=)\\d"),
    Seq_version  = str_extract(header, "(?<=SV=)\\d+"),
    Length       = width(.fasta)
  ) %>%
  filter(!str_detect(header, "^rev_|^DECOY_")) %>%
  distinct(Protein_ID, .keep_all = TRUE)

rm(.fasta, .headers)

message(sprintf("Annotation table built: %d entries", nrow(annotation_table)))


# ----------------------------------------------------------------------------
# 2. Annotate a SINGLE comparison result
#    (no list logic, no FASTA re-read - uses the annotation_table built above)
# ----------------------------------------------------------------------------
annotate_single_result <- function(df, protein_col = "Protein") {

  # Handle DDA (sp|A8K2U0|A2ML1_HUMAN) and DIA (A8K2U0) formats
  if (any(str_detect(df[[protein_col]], "\\|"), na.rm = TRUE)) {
    df[[protein_col]] <- str_extract(df[[protein_col]], "(?<=\\|)[A-Z0-9]+(?=\\|)")
  }

  df_annotated <- merge(df, annotation_table,
                        by.x  = protein_col,
                        by.y  = "Protein_ID",
                        all.x = TRUE,
                        sort  = FALSE)

  # Put annotation columns right after Protein
  annotation_cols <- c("Protein_entry", "Protein_name", "Gene_name", "Organism")
  other_cols <- setdiff(names(df_annotated), c(protein_col, annotation_cols))
  df_annotated <- df_annotated[, c(protein_col, annotation_cols, other_cols)]

  # Sort by adjusted p-value if available
  if ("adj.pvalue" %in% names(df_annotated)) {
    df_annotated <- df_annotated[order(df_annotated$adj.pvalue), ]
  }

  return(df_annotated)
}


# ----------------------------------------------------------------------------
# 3. Export ALL accumulated comparison results to a single multi-sheet xlsx
#    Called once, after the comparison loop has finished.
# ----------------------------------------------------------------------------
export_results_to_excel <- function(results_list,
                                    output_file = "comparison_results.xlsx") {

  if (length(results_list) == 0) {
    message("No results to export; skipping xlsx generation.")
    return(invisible(NULL))
  }

  wb <- createWorkbook()

  headerStyle <- createStyle(fontSize = 11, fontColour = "#FFFFFF",
                             halign = "center", fgFill = "#4F81BD",
                             border = "TopBottom", borderColour = "#4F81BD",
                             textDecoration = "bold")

  for (label in names(results_list)) {
    sheet_name <- gsub("[:\\\\/*?\\[\\]]", "_", label)
    sheet_name <- substr(sheet_name, 1, 31)  # Excel sheet-name limit

    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, results_list[[label]])

    addStyle(wb, sheet_name, headerStyle, rows = 1,
             cols = 1:ncol(results_list[[label]]), gridExpand = TRUE)
    setColWidths(wb, sheet_name, cols = 1:ncol(results_list[[label]]), widths = "auto")
    freezePane(wb, sheet_name, firstRow = TRUE)
  }

  saveWorkbook(wb, output_file, overwrite = TRUE)
  message("Results exported to: ", output_file,
          " (", length(results_list), " sheets)")

  return(invisible(output_file))
}


# ----------------------------------------------------------------------------
# 4. Export each comparison to its own CSV
#    Useful for downstream tools (GSEA loaders, etc.) that prefer flat files.
# ----------------------------------------------------------------------------
export_results_to_csv <- function(results_list, output_dir = "results") {

  if (length(results_list) == 0) {
    message("No results to export; skipping CSV generation.")
    return(invisible(NULL))
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  for (label in names(results_list)) {
    safe_name <- gsub("[:\\\\/*?\\[\\]]", "_", label)
    csv_path <- file.path(output_dir, paste0(safe_name, ".csv"))
    write.csv(results_list[[label]], csv_path, row.names = FALSE)
  }

  message("Per-comparison CSVs written to: ", output_dir,
          " (", length(results_list), " files)")

  return(invisible(output_dir))
}


# ----------------------------------------------------------------------------
# 5. Per-comparison plotting + summary statistics
#    (Body unchanged from the original 05_differential_ludovic.R)
# ----------------------------------------------------------------------------
prepare_comparison_results <- function(comparison_data, comparison_name,
                                       pvalue_cutoff = 0.05, fc_cutoff = 2,
                                       highlight_genes = c()) {

  # Add significance categories
  comparison_data$Significance <- "NS"
  comparison_data$Significance[comparison_data$adj.pvalue < pvalue_cutoff &
                                 comparison_data$log2FC > log2(fc_cutoff)] <- "Up"
  comparison_data$Significance[comparison_data$adj.pvalue < pvalue_cutoff &
                                 comparison_data$log2FC < -log2(fc_cutoff)] <- "Down"

  # Mark genes to label
  comparison_data$Label <- ""
  comparison_data$Highlight <- FALSE

  if (length(highlight_genes) > 0) {
    comparison_data$Highlight <- comparison_data$Gene_name %in% highlight_genes
    comparison_data$Label[comparison_data$Highlight] <-
      comparison_data$Gene_name[comparison_data$Highlight]
  }

  # Label top-30 most-significant up- and down-regulated proteins
  top_up <- comparison_data %>%
    dplyr::filter(Significance == "Up") %>%
    arrange(adj.pvalue)
  if (nrow(top_up) > 30) top_up <- head(top_up, 30)

  top_down <- comparison_data %>%
    dplyr::filter(Significance == "Down") %>%
    arrange(adj.pvalue)
  if (nrow(top_down) > 30) top_down <- head(top_down, 30)

  comparison_data$Label[comparison_data$Protein %in%
                          c(top_up$Protein, top_down$Protein)] <-
    comparison_data$Gene_name[comparison_data$Protein %in%
                                c(top_up$Protein, top_down$Protein)]

  # MA plot baseMean
  protein_means <- tibble(
    Protein           = rownames(protein_matrix),
    meanLog2Abundance = rowMeans(protein_matrix, na.rm = TRUE)
  )
  protein_means <- protein_means %>%
    mutate(Protein_acc = ifelse(grepl("\\|", Protein),
                                sub(".*\\|([A-Z0-9]+)\\|.*", "\\1", Protein),
                                Protein))

  comparison_data <- comparison_data %>%
    left_join(protein_means %>% select(Protein_acc, meanLog2Abundance),
              by = c("Protein" = "Protein_acc"))
  data_highlight <- comparison_data[comparison_data$Label != "" &
                                      comparison_data$Highlight == TRUE, ]
  data_regular   <- comparison_data[comparison_data$Label != "" &
                                      comparison_data$Highlight == FALSE, ]

  # ---- MA plot ----
  ma_plot <- ggplot(comparison_data, aes(x = meanLog2Abundance, y = log2FC)) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "#d62728",
                                  "Down" = "#1f77b4",
                                  "NS" = "gray70")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = c(-log2(fc_cutoff), log2(fc_cutoff)),
               linetype = "dashed", color = "gray40") +
    labs(
      title = paste("MA Plot:", comparison_name),
      subtitle = paste0("Mean expression vs. fold change (FC cutoff: ", fc_cutoff, ")"),
      x = "log2(Mean Expression)",
      y = "log2(Fold Change)",
      color = "Regulation"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  if (nrow(data_regular) > 0) {
    ma_plot <- ma_plot +
      geom_text_repel(data = data_regular, aes(label = Label),
                      size = 3.5, fontface = "plain", max.overlaps = 20,
                      box.padding = 0.5, point.padding = 0.3,
                      segment.color = "gray50", segment.size = 0.3,
                      min.segment.length = 0)
  }

  if (nrow(data_highlight) > 0) {
    ma_plot <- ma_plot +
      geom_text_repel(data = data_highlight, aes(label = Label),
                      size = 3.5, fontface = "bold", color = "black",
                      max.overlaps = 20, box.padding = 0.6, point.padding = 0.4,
                      segment.color = "gray30", segment.size = 0.5,
                      min.segment.length = 0)
  }

  # ---- Volcano plot ----
  volcano_plot <- ggplot(comparison_data, aes(x = log2FC, y = -log10(adj.pvalue))) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "#d62728",
                                  "Down" = "#1f77b4",
                                  "NS" = "gray70")) +
    scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 1.5)) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.6)) +
    geom_hline(yintercept = -log10(pvalue_cutoff),
               linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-log2(fc_cutoff), log2(fc_cutoff)),
               linetype = "dashed", color = "gray40") +
    labs(
      title = paste("Volcano Plot:", comparison_name),
      subtitle = paste0("Cutoffs: adj.p-value < ", pvalue_cutoff,
                        ", |FC| > ", fc_cutoff),
      x = "log2(Fold Change)",
      y = "-log10(Adjusted P-value)",
      color = "Regulation"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    guides(size = "none", alpha = "none")

  if (nrow(data_regular) > 0) {
    volcano_plot <- volcano_plot +
      geom_text_repel(data = data_regular, aes(label = Label),
                      size = 3.5, fontface = "plain", max.overlaps = 20,
                      box.padding = 0.5, point.padding = 0.3,
                      segment.color = "gray50", segment.size = 0.3,
                      min.segment.length = 0)
  }

  if (nrow(data_highlight) > 0) {
    volcano_plot <- volcano_plot +
      geom_text_repel(data = data_highlight, aes(label = Label),
                      size = 3.5, fontface = "bold", color = "black",
                      max.overlaps = 20, box.padding = 0.6, point.padding = 0.4,
                      segment.color = "gray30", segment.size = 0.5,
                      min.segment.length = 0)
  }

  # ---- Stats ----
  n_up    <- sum(comparison_data$Significance == "Up",   na.rm = TRUE)
  n_down  <- sum(comparison_data$Significance == "Down", na.rm = TRUE)
  n_total <- nrow(comparison_data)
  n_sig   <- n_up + n_down

  up_proteins <- comparison_data %>%
    dplyr::filter(Significance == "Up") %>%
    dplyr::select(Gene_name, Protein, Protein_entry, log2FC, adj.pvalue, Protein_name) %>%
    arrange(adj.pvalue)

  down_proteins <- comparison_data %>%
    dplyr::filter(Significance == "Down") %>%
    dplyr::select(Gene_name, Protein, Protein_entry, log2FC, adj.pvalue, Protein_name) %>%
    arrange(adj.pvalue)

  return(list(
    volcano_plot   = volcano_plot,
    ma_plot        = ma_plot,
    n_up           = n_up,
    n_down         = n_down,
    n_total        = n_total,
    n_sig          = n_sig,
    up_proteins    = up_proteins,
    down_proteins  = down_proteins,
    pvalue_cutoff  = pvalue_cutoff,
    fc_cutoff      = fc_cutoff
  ))
}


# ----------------------------------------------------------------------------
# 6. Initialize the accumulator the child template will append to
# ----------------------------------------------------------------------------
all_annotated_results <- list()
