# Fixed timestamp and user as provided
timestamp <- "2025-05-29 01:00:05"
user <- "hingelman"

# Create output directory for significant genes
base_dir <- "statisticalandfunctionalanalysis"
sig_genes_dir <- file.path(base_dir, paste0("significant_genes_", format(Sys.time(), "%Y%m%d_%H%M%S")))
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

# Function to extract significant genes with debugging
analyze_deg <- function(condition_pair) {
  condition1 <- condition_pair[1]
  condition2 <- condition_pair[2]
  
  # Print debug info
  cat(sprintf("\nAnalyzing %s vs %s:\n", condition1, condition2))
  
  # Create safe filenames
  safe_name1 <- make_safe_filename(condition1)
  safe_name2 <- make_safe_filename(condition2)
  
  # Get results with debugging info
  cat("- Getting DESeq2 results... ")
  res <- results(dds, contrast = c("physiological_state", condition1, condition2))
  cat("done\n")
  
  # Print initial counts
  cat(sprintf("- Total genes in results: %d\n", nrow(res)))
  cat(sprintf("- NA values in padj: %d\n", sum(is.na(res$padj))))
  cat(sprintf("- NA values in log2FoldChange: %d\n", sum(is.na(res$log2FoldChange))))
  
  # Convert to data frame and add gene_id
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id")
  
  # Print filtering steps
  cat("Filtering steps:\n")
  
  # First filter: padj < 0.05
  genes_sig_p <- sum(res_df$padj < 0.05, na.rm = TRUE)
  cat(sprintf("- Genes with padj < 0.05: %d\n", genes_sig_p))
  
  # Second filter: |log2FC| >= 1
  genes_sig_fc <- sum(abs(res_df$log2FoldChange) >= 1, na.rm = TRUE)
  cat(sprintf("- Genes with |log2FC| >= 1: %d\n", genes_sig_fc))
  
  # Both filters
  sig_genes <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%  # Remove NA values explicitly
    filter(padj < 0.05, abs(log2FoldChange) >= 1) %>%
    arrange(padj)
  
  cat(sprintf("- Genes passing both filters: %d\n", nrow(sig_genes)))
  
  # Add comparison column and save results
  if(nrow(sig_genes) > 0) {
    sig_genes$comparison <- paste(condition1, "vs", condition2)
    output_file <- file.path(sig_genes_dir, 
                            paste0("significant_genes_", safe_name1, "vs", safe_name2, ".csv"))
    write.csv(sig_genes, output_file, row.names = FALSE)
    cat(sprintf("- Results saved to: %s\n", output_file))
    
    # Print top 5 most significant genes
    cat("Top 5 most significant genes:\n")
    print(head(sig_genes[, c("gene_id", "log2FoldChange", "padj")], 5))
  }
  
  # Create summary stats
  summary_stats <- list(
    comparison = paste(condition1, "vs", condition2),
    total_genes = nrow(res),
    significant_genes = nrow(sig_genes),
    upregulated = sum(sig_genes$log2FoldChange > 0),
    downregulated = sum(sig_genes$log2FoldChange < 0)
  )
  
  return(summary_stats)
}

# Create log file with debugging info
log_file <- file.path(base_dir, paste0("significant_genes_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_log.txt"))
sink(log_file, split = TRUE)

# Print header with DESeq2 info
cat("=== Significant Genes Analysis Log ===\n")
cat(sprintf("Date/Time (UTC): %s\n", timestamp))
cat(sprintf("User: %s\n", user))
cat(sprintf("Output Directory: %s\n\n", normalizePath(sig_genes_dir)))

# Print DESeq2 object info
cat("DESeq2 Object Information:\n")
cat(sprintf("- Number of genes: %d\n", nrow(dds)))
cat(sprintf("- Number of samples: %d\n", ncol(dds)))
cat("- Design formula:", deparse(design(dds)), "\n")
cat("- Reference level:", reference_condition, "\n\n")

# Print all conditions
cat("All conditions:\n")
for(cond in conditions) {
  cat(sprintf("- %s%s\n", cond, if(cond == reference_condition) " (reference)" else ""))
}
cat("\n")

# Print analysis parameters
cat("Analysis Parameters:\n")
cat("- Significance cutoff (padj): < 0.05\n")
cat("- Log2 fold change cutoff: >= 1 (absolute value)\n\n")

# Run analysis with detailed output
cat("Processing comparisons...\n")
summaries <- lapply(all_combinations, analyze_deg)

# Create overall summary
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