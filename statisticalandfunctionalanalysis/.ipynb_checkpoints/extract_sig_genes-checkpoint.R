# Fixed timestamp and user info
timestamp <- "2025-05-29 11:21:38"
user <- "hingelman"

# Create output directory for significant genes
base_dir <- "condition_comparisons"
sig_genes_dir <- file.path(base_dir, paste0("significant_DE_genes", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(sig_genes_dir, recursive = TRUE, showWarnings = FALSE)

# Function to create safe filenames
make_safe_filename <- function(name) {
  gsub("[^a-zA-Z0-9]", "_", name)
}

# Get conditions and reference condition
conditions <- levels(dds$physiological_state)
reference_condition <- conditions[1]

# Create combinations: reference vs all others + all other pairs
reference_combinations <- lapply(conditions[-1], function(cond) c(cond, reference_condition))
other_combinations <- combn(conditions[-1], 2, simplify = FALSE)
all_combinations <- c(reference_combinations, other_combinations)

# Function to extract significant genes
analyze_deg <- function(condition_pair) {
  tryCatch({
    condition1 <- condition_pair[1]
    condition2 <- condition_pair[2]
    
    cat(sprintf("\nAnalyzing %s vs %s:\n", condition1, condition2))
    
    # Create safe filenames
    safe_name1 <- make_safe_filename(condition1)
    safe_name2 <- make_safe_filename(condition2)
    
    # Get results with padj threshold
    res <- results(dds, 
                  contrast = c("physiological_state", condition1, condition2),
                  alpha = 0.05)
    
    # Convert to data frame and filter
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    
    # Print diagnostic info
    cat(sprintf("Total genes: %d\n", nrow(res_df)))
    cat(sprintf("Genes with non-NA padj: %d\n", sum(!is.na(res_df$padj))))
    cat(sprintf("Genes with padj < 0.05: %d\n", sum(res_df$padj < 0.05, na.rm = TRUE)))
    cat(sprintf("Genes with |log2FC| >= 1: %d\n", sum(abs(res_df$log2FoldChange) >= 1, na.rm = TRUE)))
    
    # Filter for significant genes
    sig_genes <- res_df %>%
      filter(!is.na(padj), !is.na(log2FoldChange)) %>%
      filter(padj < 0.05, abs(log2FoldChange) >= 1) %>%
      arrange(padj)
    
    # Create output file path
    output_file <- file.path(sig_genes_dir, 
                            paste0("significant_genes_", safe_name1, "vs", safe_name2, ".csv"))
    
    # Process and save results
    if(nrow(sig_genes) > 0) {
      # Add comparison information
      sig_genes$comparison <- paste(condition1, "vs", condition2)
      
      # Save to CSV
      write.csv(sig_genes, output_file, row.names = FALSE)
      
      # Print summary
      cat("\nSignificant genes summary:\n")
      cat(sprintf("- Total significant: %d\n", nrow(sig_genes)))
      cat(sprintf("- Upregulated: %d\n", sum(sig_genes$log2FoldChange > 0)))
      cat(sprintf("- Downregulated: %d\n", sum(sig_genes$log2FoldChange < 0)))
      
      # Print top genes
      cat("\nTop 5 most significant genes:\n")
      print(head(sig_genes[, c("gene_id", "log2FoldChange", "padj")], 5))
    } else {
      # Create empty file with headers
      empty_df <- data.frame(
        gene_id = character(),
        baseMean = numeric(),
        log2FoldChange = numeric(),
        lfcSE = numeric(),
        stat = numeric(),
        pvalue = numeric(),
        padj = numeric(),
        comparison = character(),
        stringsAsFactors = FALSE
      )
      write.csv(empty_df, output_file, row.names = FALSE)
      cat("\nNo genes passed significance thresholds (padj < 0.05 & |log2FC| >= 1)\n")
    }
    
    cat(sprintf("\nResults saved to: %s\n", output_file))
    
    return(list(
      comparison = paste(condition1, "vs", condition2),
      total_genes = nrow(res),
      significant_genes = nrow(sig_genes),
      upregulated = if(nrow(sig_genes) > 0) sum(sig_genes$log2FoldChange > 0) else 0,
      downregulated = if(nrow(sig_genes) > 0) sum(sig_genes$log2FoldChange < 0) else 0
    ))
    
  }, error = function(e) {
    cat(sprintf("\nError in analyzing %s vs %s: %s\n", 
                condition_pair[1], condition_pair[2], as.character(e)))
    return(list(
      comparison = paste(condition_pair[1], "vs", condition_pair[2]),
      total_genes = 0,
      significant_genes = 0,
      upregulated = 0,
      downregulated = 0,
      error = as.character(e)
    ))
  })
}

# Create log file
log_file <- file.path(base_dir, paste0("significant_DE_genes_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_log.txt"))
sink(log_file, split = TRUE)

# Print header
cat("=== Significant Genes Analysis Log ===\n")
cat(sprintf("Date/Time (UTC): %s\n", timestamp))
cat(sprintf("User: %s\n", user))
cat(sprintf("Output Directory: %s\n\n", normalizePath(sig_genes_dir)))

# Print DESeq2 info
cat("DESeq2 Object Information:\n")
cat(sprintf("- Number of genes: %d\n", nrow(dds)))
cat(sprintf("- Number of samples: %d\n", ncol(dds)))
cat("- Design formula:", deparse(design(dds)), "\n")
cat("- Reference level:", reference_condition, "\n\n")

# Print conditions
cat("All conditions:\n")
for(cond in conditions) {
  cat(sprintf("- %s%s\n", cond, if(cond == reference_condition) " (reference)" else ""))
}
cat("\n")

# Print analysis parameters
cat("Analysis Parameters:\n")
cat("- Significance threshold (padj): < 0.05\n")
cat("- Log2 fold change threshold: >= 1 (absolute value)\n\n")

# Run analysis
cat("Processing comparisons...\n")
summaries <- lapply(all_combinations, analyze_deg)

# Print overall summary
cat("\nOverall Summary:\n")
cat(sprintf("Total number of comparisons: %d\n", length(all_combinations)))
cat("Results by comparison:\n")
for(summary in summaries) {
  cat(sprintf("- %s: %d significant genes (Up: %d, Down: %d)\n",
              summary$comparison, summary$significant_genes,
              summary$upregulated, summary$downregulated))
}

# Print completion
cat(sprintf("\nAnalysis completed at %s UTC\n", timestamp))

# Close log file
sink()

# Print final message to console
cat(sprintf("\nAnalysis complete!\n"))
cat(sprintf("- Results saved in: %s\n", normalizePath(sig_genes_dir)))
cat(sprintf("- Log file: %s\n", normalizePath(log_file)))

# Return the directory paths invisibly
invisible(list(
  results_dir = sig_genes_dir,
  log_file = log_file
))