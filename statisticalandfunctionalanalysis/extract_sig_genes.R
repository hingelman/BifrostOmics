# Fixed timestamp and user
timestamp <- "2025-05-29 01:08:11"
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

# Function to extract significant genes
analyze_deg <- function(condition_pair) {
  condition1 <- condition_pair[1]
  condition2 <- condition_pair[2]
  
  cat(sprintf("\nAnalyzing %s vs %s:\n", condition1, condition2))
  
  # Create safe filenames
  safe_name1 <- make_safe_filename(condition1)
  safe_name2 <- make_safe_filename(condition2)
  
  # Get results with alpha=1 to get all p-values
  res <- results(dds, 
                contrast = c("physiological_state", condition1, condition2),
                alpha = 1)  # Get all results without filtering
  
  # Convert to data frame and add gene_id
  sig_genes <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    filter(padj < 0.05, abs(log2FoldChange) >= 1) %>%
    arrange(padj)
  
  # Add comparison information
  sig_genes$comparison <- paste(condition1, "vs", condition2)
  
  # Save results
  output_file <- file.path(sig_genes_dir, 
                          paste0("significant_genes_", safe_name1, "vs", safe_name2, ".csv"))
  write.csv(sig_genes, output_file, row.names = FALSE)
  
  # Print summary
  cat(sprintf("- Total genes analyzed: %d\n", nrow(res)))
  cat(sprintf("- Significant genes (padj < 0.05 & |log2FC| >= 1): %d\n", nrow(sig_genes)))
  cat(sprintf("- Upregulated: %d\n", sum(sig_genes$log2FoldChange > 0)))
  cat(sprintf("- Downregulated: %d\n", sum(sig_genes$log2FoldChange < 0)))
  
  # Print top 5 most significant genes
  if(nrow(sig_genes) > 0) {
    cat("\nTop 5 most significant genes:\n")
    print(head(sig_genes[, c("gene_id", "log2FoldChange", "padj")], 5))
  }
  
  return(list(
    comparison = paste(condition1, "vs", condition2),
    total_genes = nrow(res),
    significant_genes = nrow(sig_genes),
    upregulated = sum(sig_genes$log2FoldChange > 0),
    downregulated = sum(sig_genes$log2FoldChange < 0)
  ))
}

# Create log file
log_file <- file.path(base_dir, paste0("significant_genes_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_log.txt"))
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
cat("- Significance cutoff (padj): < 0.05\n")
cat("- Log2 fold change cutoff: >= 1 (absolute value)\n")
cat("- Using alpha=1 in results() to get all genes\n\n")

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