# Get timestamp and user info
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
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
reference_condition <- conditions[1]  # The first condition is the reference

# Create combinations: reference vs all others + all other pairs
reference_combinations <- lapply(conditions[-1], function(cond) c(cond, reference_condition))
other_combinations <- combn(conditions[-1], 2, simplify = FALSE)
all_combinations <- c(reference_combinations, other_combinations)

# Function to extract significant genes
analyze_deg <- function(condition_pair) {
  condition1 <- condition_pair[1]
  condition2 <- condition_pair[2]
  
  # Create safe filenames
  safe_name1 <- make_safe_filename(condition1)
  safe_name2 <- make_safe_filename(condition2)
  
  # Get results
  res <- results(dds, contrast = c("physiological_state", condition1, condition2))
  
  # Extract significant genes using dplyr
  sig_genes <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    filter(padj < 0.05, abs(log2FoldChange) >= 1) %>%
    arrange(padj)  # Sort by adjusted p-value
  
  # Add comparison information
  sig_genes$comparison <- paste(condition1, "vs", condition2)
  
  # Save results
  output_file <- file.path(sig_genes_dir, 
                          paste0("significant_genes_", safe_name1, "vs", safe_name2, ".csv"))
  write.csv(sig_genes, output_file, row.names = FALSE)
  
  # Create summary statistics
  summary_stats <- list(
    comparison = paste(condition1, "vs", condition2),
    total_genes = nrow(res),
    significant_genes = nrow(sig_genes),
    upregulated = sum(sig_genes$log2FoldChange > 0),
    downregulated = sum(sig_genes$log2FoldChange < 0)
  )
  
  # Print progress
  cat(sprintf("Processed %s vs %s:\n", condition1, condition2))
  cat(sprintf("- Total genes analyzed: %d\n", summary_stats$total_genes))
  cat(sprintf("- Significant genes: %d\n", summary_stats$significant_genes))
  cat(sprintf("- Upregulated: %d\n", summary_stats$upregulated))
  cat(sprintf("- Downregulated: %d\n", summary_stats$downregulated))
  cat(sprintf("- Results saved to: %s\n\n", output_file))
  
  return(summary_stats)
}

# Create log file
log_file <- file.path(base_dir, paste0("significant_genes_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_log.txt"))
sink(log_file, split = TRUE)

# Print header
cat("=== Significant Genes Analysis Log ===\n")
cat(sprintf("Date/Time (UTC): %s\n", timestamp))
cat(sprintf("User: %s\n", user))
cat(sprintf("Output Directory: %s\n\n", normalizePath(sig_genes_dir)))

# Print reference condition
cat(sprintf("Reference condition: %s\n\n", reference_condition))

# Print conditions
cat("All conditions:\n")
for(cond in conditions) {
  cat(sprintf("- %s%s\n", cond, if(cond == reference_condition) " (reference)" else ""))
}
cat("\n")

# Print analysis parameters
cat("Analysis Parameters:\n")
cat("- Significance cutoff (padj): < 0.05\n")
cat("- Log2 fold change cutoff: >= 1 (absolute value)\n\n")

# Run analysis and collect summaries
cat("Processing comparisons...\n\n")
summaries <- lapply(all_combinations, analyze_deg)

# Create overall summary
cat("\nOverall Summary:\n")
cat(sprintf("Total number of comparisons: %d\n", length(all_combinations)))
total_sig_genes <- sum(sapply(summaries, function(x) x$significant_genes))
cat(sprintf("Total significant genes across all comparisons: %d\n", total_sig_genes))

# Print completion message
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