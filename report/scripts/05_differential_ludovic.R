library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)

# contrast matrix for the single comparison
contrast_mat <- matrix(c(-1, 1), nrow = 1)
colnames(contrast_mat) <- conditions
rownames(contrast_mat) <- comparison_name

# Fix direction: condition_A (first) gets +1, condition_B (second) gets -1
cond_A <- conditions[1]
cond_B <- conditions[2]
contrast_mat[1, cond_A] <- 1
contrast_mat[1, cond_B] <- -1

test_results <- groupComparison(contrast.matrix = contrast_mat,
                                data = summarized)



# Function to annotate proteins with gene names and Ensembl IDs using biomaRt
annotate_protein_results <- function(results_list, 
                                     protein_col = "Protein") {
  
  # Load required packages
  if (!require("tidyverse")) {
    stop("Package 'tidyverse' is required. Install with: install.packages('tidyverse')")
  }
  
  if (!require("Biostrings")) {
    stop("Package 'Biostrings' is required. Install with: BiocManager::install('Biostrings')")
  }
  message("Extracting annotations from FASTA database...")
  
  # Read FASTA
  fasta <- readAAStringSet("/home/jszymanski/Documents/nextflow_pipelines_proteomics/fasta_databases/2026-03-18-decoys-contam-uniprotkb_proteome_UP000005640_AND_revi_2026_03_18.fasta.fas")
  
  # Extract headers
  headers <- names(fasta)
  
  # Parse into annotation table
  annotation <- tibble(header = headers) %>%
    mutate(
      db          = str_extract(header, "^(sp|tr)"),
      Protein_ID  = str_extract(header, "(?<=\\|)[A-Z0-9]+(?=\\|)"),
      Protein_entry  = str_extract(header, "(?<=\\|)[A-Z0-9_]+(?= )"),
      Protein_name = str_extract(header, "(?<= )[^|]+?(?= OS=)"),
      Organism    = str_extract(header, "(?<=OS=).+?(?= OX=)"),
      TaxID       = str_extract(header, "(?<=OX=)\\d+"),
      Gene_name   = str_extract(header, "(?<=GN=)\\S+"),
      PE_level    = str_extract(header, "(?<=PE=)\\d"),
      Seq_version = str_extract(header, "(?<=SV=)\\d+"),
      Length      = width(fasta)
    ) %>%
    filter(!str_detect(header, "^rev_|^DECOY_")) %>%
    distinct(Protein_ID, .keep_all = TRUE)
  
  
  # Annotate each dataframe in the list
  annotated_list <- lapply(results_list, function(df) {
    
    # Handle both DDA (sp|A8K2U0|A2ML1_HUMAN) and DIA (A8K2U0) formats
    if (any(str_detect(df$Protein, "\\|"), na.rm = TRUE)) {
      df <- df %>%
        mutate(Protein = str_extract(df$Protein, "(?<=\\|)[A-Z0-9]+(?=\\|)"))
    } 
    
    df_annotated <- merge(df, annotation, 
                          by.x = protein_col, 
                          by.y = "Protein_ID", 
                          all.x = TRUE, 
                          sort = FALSE)
    
    # Reorder columns to put annotations right after Protein
    annotation_cols <- c("Protein_entry", "Protein_name", "Gene_name", "Organism")
    other_cols <- setdiff(names(df_annotated), c(protein_col, annotation_cols))
    df_annotated <- df_annotated[, c(protein_col, annotation_cols, other_cols)]
    
    # Sort by adjusted p-value if available
    if ("adj.pvalue" %in% names(df_annotated)) {
      df_annotated <- df_annotated[order(df_annotated$adj.pvalue), ]
    }
    
    return(df_annotated)
  })
  
  message("Annotation complete from FASTA")
  
  return(annotated_list)
}


# Function to split comparison results into list
split_comparison_results <- function(test_results) {
  
  # Extract comparison results
  results_df <- test_results$ComparisonResult
  
  # Get unique comparison labels
  comparison_labels <- unique(results_df$Label)
  
  # Split results into list of dataframes
  results_list <- list()
  
  for (label in comparison_labels) {
    results_list[[label]] <- results_df[results_df$Label == label, ]
  }
  
  message("Split results into ", length(comparison_labels), " comparisons")
  message("Comparison names: ", paste(names(results_list), collapse = ", "))
  
  return(results_list)
}


# Function to export results list to Excel
export_results_to_excel <- function(results_list, output_file = "comparison_results.xlsx") {
  
  # Load required package
  if (!require("openxlsx")) {
    stop("Package 'openxlsx' is required. Install with: install.packages('openxlsx')")
  }
  
  # Create workbook
  wb <- createWorkbook()
  
  # Add separate sheet for each comparison
  for (label in names(results_list)) {
    # Clean sheet name (Excel has restrictions)
    sheet_name <- gsub("[:\\\\/*?\\[\\]]", "_", label)
    sheet_name <- substr(sheet_name, 1, 31)  # Excel sheet name limit
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, results_list[[label]])
    
    # Format header
    headerStyle <- createStyle(fontSize = 11, fontColour = "#FFFFFF", 
                               halign = "center", fgFill = "#4F81BD",
                               border = "TopBottom", borderColour = "#4F81BD",
                               textDecoration = "bold")
    addStyle(wb, sheet_name, headerStyle, rows = 1, 
             cols = 1:ncol(results_list[[label]]), gridExpand = TRUE)
    
    # Auto-size columns
    setColWidths(wb, sheet_name, cols = 1:ncol(results_list[[label]]), widths = "auto")
    
    # Freeze top row
    freezePane(wb, sheet_name, firstRow = TRUE)
  }
  
  # Save workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  message("Results exported to: ", output_file)
  message("Number of sheets: ", length(results_list))
  
  return(invisible(output_file))
}

# Split results into list
results_list <- split_comparison_results(test_results)

# Annotate with biomaRt
annotated_results <- annotate_protein_results(
  results_list,
  protein_col = "Protein"
)

# Export to Excel
export_results_to_excel(
  annotated_results,
  output_file = "comparison_results.xlsx"
)





# Function to prepare comparison data and plots
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
    comparison_data$Label[comparison_data$Highlight] <- comparison_data$Gene_name[comparison_data$Highlight]
  }
  
  # Label top significant genes
  top_up <- comparison_data %>%
    dplyr::filter(Significance == "Up") %>%
    arrange(adj.pvalue)
  
  if (nrow(top_up) > 30) {
    top_up <- head(top_up, 30)
  }
  
  top_down <- comparison_data %>%
    dplyr::filter(Significance == "Down") %>%
    arrange(adj.pvalue)
  
  if (nrow(top_down) > 30) {
    top_down <- head(top_down, 30)
  }
  
  comparison_data$Label[comparison_data$Protein %in% c(top_up$Protein, top_down$Protein)] <- 
    comparison_data$Gene_name[comparison_data$Protein %in% c(top_up$Protein, top_down$Protein)]
  
  
  
  # Create MA plot
  comparison_data$baseMean <- 2^rowMeans(cbind(
    comparison_data$log2FC + log2(comparison_data$SE),
    log2(comparison_data$SE)
  ), na.rm = TRUE)
  
  # Prepare data subsets for labeling
  data_highlight <- comparison_data[comparison_data$Label != "" & comparison_data$Highlight == TRUE, ]
  data_regular <- comparison_data[comparison_data$Label != "" & comparison_data$Highlight == FALSE, ]
  
  ma_plot <- ggplot(comparison_data, aes(x = log10(baseMean), y = log2FC)) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "#d62728", "Down" = "#1f77b4", "NS" = "gray70")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = c(-log2(fc_cutoff), log2(fc_cutoff)), 
               linetype = "dashed", color = "gray40") +
    labs(
      title = paste("MA Plot:", comparison_name),
      subtitle = paste0("Mean expression vs. fold change (FC cutoff: ", fc_cutoff, ")"),
      x = "log10(Mean Expression)",
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
  
  # Add regular labels layer
  if (nrow(data_regular) > 0) {
    ma_plot <- ma_plot +
      geom_text_repel(
        data = data_regular,
        aes(label = Label),
        size = 3.5,
        fontface = "plain",
        max.overlaps = 20,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "gray50",
        segment.size = 0.3,
        min.segment.length = 0
      )
  }
  
  # Add highlighted labels layer (bold and larger)
  if (nrow(data_highlight) > 0) {
    ma_plot <- ma_plot +
      geom_text_repel(
        data = data_highlight,
        aes(label = Label),
        size = 3.5,
        fontface = "bold",
        color = "black",
        max.overlaps = 20,
        box.padding = 0.6,
        point.padding = 0.4,
        segment.color = "gray30",
        segment.size = 0.5,
        min.segment.length = 0
      )
  }
  
  
  # Create base volcano plot
  volcano_plot <- ggplot(comparison_data, aes(x = log2FC, y = -log10(adj.pvalue))) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "#d62728", "Down" = "#1f77b4", "NS" = "gray70")) +
    scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 1.5)) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.6)) +
    geom_hline(yintercept = -log10(pvalue_cutoff), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-log2(fc_cutoff), log2(fc_cutoff)), linetype = "dashed", color = "gray40") +
    labs(
      title = paste("Volcano Plot:", comparison_name),
      subtitle = paste0("Cutoffs: adj.p-value < ", pvalue_cutoff, ", |FC| > ", fc_cutoff),
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
  
  # Add regular labels layer
  if (nrow(data_regular) > 0) {
    volcano_plot <- volcano_plot +
      geom_text_repel(
        data = data_regular,
        aes(label = Label),
        size = 3.5,
        fontface = "plain",
        max.overlaps = 20,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "gray50",
        segment.size = 0.3,
        min.segment.length = 0
      )
  }
  
  # Add highlighted labels layer (bold and larger)
  if (nrow(data_highlight) > 0) {
    volcano_plot <- volcano_plot +
      geom_text_repel(
        data = data_highlight,
        aes(label = Label),
        size = 3.5,
        fontface = "bold",
        color = "black",
        max.overlaps = 20,
        box.padding = 0.6,
        point.padding = 0.4,
        segment.color = "gray30",
        segment.size = 0.5,
        min.segment.length = 0
      )
  }
  
  
  # Calculate statistics
  n_up <- sum(comparison_data$Significance == "Up", na.rm = TRUE)
  n_down <- sum(comparison_data$Significance == "Down", na.rm = TRUE)
  n_total <- nrow(comparison_data)
  n_sig <- n_up + n_down
  
  # Get up and down regulated proteins
  up_proteins <- comparison_data %>%
    dplyr::filter(Significance == "Up") %>%
    dplyr::select(Gene_name, Protein, Protein_entry, log2FC, adj.pvalue, Protein_name) %>%
    arrange(adj.pvalue)
  
  down_proteins <- comparison_data %>%
    dplyr::filter(Significance == "Down") %>%
    dplyr::select(Gene_name, Protein, Protein_entry, log2FC, adj.pvalue, Protein_name) %>%
    arrange(adj.pvalue)
  
  # Return list of results
  return(list(
    volcano_plot = volcano_plot,
    ma_plot = ma_plot,
    n_up = n_up,
    n_down = n_down,
    n_total = n_total,
    n_sig = n_sig,
    up_proteins = up_proteins,
    down_proteins = down_proteins,
    pvalue_cutoff = pvalue_cutoff,
    fc_cutoff = fc_cutoff
  ))
}