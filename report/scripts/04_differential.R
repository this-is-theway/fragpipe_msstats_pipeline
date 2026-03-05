library(openxlsx)
library(ggplot2)
library(ggrepel)
library(dplyr)


# Read the comparison file
# Uses `comparisons_file` variable set by run_msstats.Rmd from params
comparisons_df <- read.csv(comparisons_file, stringsAsFactors = FALSE, sep = ",")

# Function to create contrast matrix from CSV file
create_contrast_from_csv <- function(comparisons_df, conditions = NULL, 
                                     annotation = NULL) {
      
      # Check for required columns
      if (!all(c("condition_A", "condition_B") %in% colnames(comparisons_df))) {
        stop("CSV file must contain 'condition_A' and 'condition_B' columns")
      }
      
      # Clean up condition names (remove extra whitespace)
      comparisons_df$condition_A <- trimws(comparisons_df$condition_A)
      comparisons_df$condition_B <- trimws(comparisons_df$condition_B)
      
      # Get conditions from annotation if not provided
      if (is.null(conditions)) {
        if (is.null(annotation)) {
          # Extract all unique conditions from comparison file
          conditions <- unique(c(comparisons_df$condition_A, 
                                comparisons_df$condition_B))
        } else {
          # Get from annotation
          conditions <- unique(annotation$Condition)
        }
      }
      
      # Sort conditions for consistency
      conditions <- sort(conditions)
      
      # Initialize contrast matrix
      n_comparisons <- nrow(comparisons_df)
      n_conditions <- length(conditions)
      contrast_matrix <- matrix(0, nrow = n_comparisons, ncol = n_conditions)
      
      # Set column names
      colnames(contrast_matrix) <- conditions
      
      # Create row names and fill matrix
      row_names <- character(n_comparisons)
      
      for (i in 1:n_comparisons) {
        cond_A <- comparisons_df$condition_A[i]
        cond_B <- comparisons_df$condition_B[i]
        
        # Check if conditions exist
        if (!(cond_A %in% conditions)) {
          stop(paste("Condition", cond_A, "not found in available conditions"))
        }
        if (!(cond_B %in% conditions)) {
          stop(paste("Condition", cond_B, "not found in available conditions"))
        }
        
        # Create comparison name: ConditionA_vs_ConditionB
        row_names[i] <- paste0(cond_A, "_vs_", cond_B)
        
        # Set contrast values (A is numerator, B is denominator)
        contrast_matrix[i, cond_B] <- -1
        contrast_matrix[i, cond_A] <- 1
      }
      
      # Set row names
      rownames(contrast_matrix) <- row_names
      
      return(contrast_matrix)
}


contrast_mat <- create_contrast_from_csv(comparisons_df, 
                                         annotation = metadata)

test_results <- groupComparison(contrast.matrix = contrast_mat,
                                data = summarized)




# Function to annotate proteins with gene names and Ensembl IDs using biomaRt
annotate_protein_results <- function(results_list, 
                                     protein_col = "Protein",
                                     species = "human") {
        
        # Load required package
        if (!require("biomaRt")) {
          stop("Package 'biomaRt' is required. Install with: BiocManager::install('biomaRt')")
        }
        
        message("Fetching annotations from Ensembl via biomaRt...")
        
        # Set up Ensembl connection based on species
        if (species == "human") {
          ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
        } else if (species == "mouse") {
          ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
        } else if (species == "rat") {
          ensembl <- useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl")
        } else {
          stop("Supported species: 'human', 'mouse', or 'rat'")
        }
        
        # Get all unique proteins from all comparisons
        all_proteins <- unique(unlist(lapply(results_list, function(df) df[[protein_col]])))
        
        message("Found ", length(all_proteins), " unique proteins to annotate")
        
        # Extract UniProt IDs (remove isoform info if present)
        uniprot_ids <- gsub("-\\d+$", "", all_proteins)  # Remove -1, -2, etc.
        uniprot_ids <- gsub(";.*", "", uniprot_ids)      # Take first ID if multiple
        uniprot_ids <- unique(uniprot_ids)
        
        # Fetch annotations from Ensembl
        protein_annotation <- NULL
        tryCatch({
          protein_annotation <- getBM(
            attributes = c('uniprotswissprot', 'external_gene_name', 'ensembl_gene_id', 'description'),
            filters = 'uniprotswissprot',
            values = uniprot_ids,
            mart = ensembl
          )
          
          # Rename columns
          colnames(protein_annotation) <- c("Protein_ID", "Gene_name", "Ensembl_ID", "Description")
          
          # Remove gene description prefix and clean up
          protein_annotation$Description <- gsub("\\[Source:.*\\]$", "", protein_annotation$Description)
          protein_annotation$Description <- trimws(protein_annotation$Description)
          
        }, error = function(e) {
          message("Warning: biomaRt query failed. Trying without description...")
          protein_annotation <<- getBM(
            attributes = c('uniprotswissprot', 'external_gene_name', 'ensembl_gene_id'),
            filters = 'uniprotswissprot',
            values = uniprot_ids,
            mart = ensembl
          )
          colnames(protein_annotation) <<- c("Protein_ID", "Gene_name", "Ensembl_ID")
        })
        
        # Create mapping including original protein IDs (with isoforms)
        protein_mapping <- data.frame(
          Protein = all_proteins,
          Protein_ID_clean = gsub("-\\d+$", "", all_proteins),
          stringsAsFactors = FALSE
        )
        
        protein_mapping <- merge(protein_mapping, protein_annotation, 
                                 by.x = "Protein_ID_clean", 
                                 by.y = "Protein_ID", 
                                 all.x = TRUE)
        protein_mapping$Protein_ID_clean <- NULL  # Remove intermediate column
        
        # Annotate each dataframe in the list
        annotated_list <- lapply(results_list, function(df) {
          df_annotated <- merge(df, protein_mapping, 
                                by.x = protein_col, 
                                by.y = "Protein", 
                                all.x = TRUE, 
                                sort = FALSE)
          
          # Reorder columns to put annotations right after Protein
          if ("Description" %in% names(df_annotated)) {
            annotation_cols <- c("Gene_name", "Ensembl_ID", "Description")
          } else {
            annotation_cols <- c("Gene_name", "Ensembl_ID")
          }
          
          other_cols <- setdiff(names(df_annotated), c(protein_col, annotation_cols))
          df_annotated <- df_annotated[, c(protein_col, annotation_cols, other_cols)]
          
          # Sort by adjusted p-value if available
          if ("adj.pvalue" %in% names(df_annotated)) {
            df_annotated <- df_annotated[order(df_annotated$adj.pvalue), ]
          }
          
          return(df_annotated)
        })
        
        message("Annotation complete via biomaRt")
        message("Annotated ", nrow(protein_annotation), " out of ", length(all_proteins), " proteins")
        
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
# Uses `species` variable set by run_msstats.Rmd from params
annotated_results <- annotate_protein_results(
  results_list,
  protein_col = "Protein",
  species = species
)

# Export to Excel
export_results_to_excel(
  annotated_results,
  output_file = "comparison_results.xlsx"
)





# Function to prepare comparison data and plots
prepare_comparison_results <- function(comparison_data, comparison_name, 
                                       pvalue_cutoff = 0.05, fc_cutoff = 1.5,
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
        
        if (length(top_up) > 30) {
          top_up <- head(top_up, 30)
        }
        
        top_down <- comparison_data %>%
          dplyr::filter(Significance == "Down") %>%
          arrange(adj.pvalue)
        
        if (length(top_down) > 30) {
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
          geom_point(aes(color = Significance, size = Highlight, alpha = Highlight)) +
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
          dplyr::select(Gene_name, Protein, Ensembl_ID, log2FC, adj.pvalue, Description) %>%
          arrange(adj.pvalue)
        
        down_proteins <- comparison_data %>%
          dplyr::filter(Significance == "Down") %>%
          dplyr::select(Gene_name, Protein, Ensembl_ID, log2FC, adj.pvalue, Description) %>%
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
