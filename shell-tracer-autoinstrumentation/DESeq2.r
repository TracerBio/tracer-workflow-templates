# Load necessary libraries
if (!requireNamespace("edgeR", quietly = TRUE)) {
  install.packages("edgeR", repos='http://cran.rstudio.com/')
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos='http://cran.rstudio.com/')
}

if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2", repos='http://cran.rstudio.com/')
}

library(edgeR)
library(ggplot2)
library(reshape2)

# Set working directory to where the count files are
setwd("/workspace/tracer-workflow-templates/data")

# List all count files
files <- list.files(pattern = "*_counts.txt$", full.names = TRUE)

# Check if files are loaded correctly
print("Count files found:")
print(files)

# Define sample names based on the number of count files detected
sampleNames <- c("control", "test")

# Read count data into edgeR format
countData <- lapply(files, function(x) {
  data <- read.table(x, header = TRUE, row.names = 1, sep = "\t")  # Ensure correct delimiter
  counts <- data[, ncol(data)]  # Take the count column
  names(counts) <- rownames(data)
  return(counts)
})

# Combine into a single matrix
countMatrix <- do.call(cbind, countData)
print("Count matrix dimensions:")
print(dim(countMatrix))

# Verify consistency of row names across all files
consistent_rows <- Reduce(intersect, lapply(countData, names))
print("Number of consistent rows (genes) across samples:")
print(length(consistent_rows))

# Ensure the count matrix is constructed using consistent row names
countMatrix <- countMatrix[consistent_rows, , drop=FALSE]

# Visualize distribution of counts to understand data before filtering
boxplot(countMatrix, main = "Boxplot of raw counts", las = 2)

# Define group for each sample
group <- factor(sampleNames)

# Create a DGEList
dge <- DGEList(counts = countMatrix, group = group)

# Check the DGEList object
print("DGEList object:")
print(dge)

# Explore filtering threshold for low counts
# Print number of genes before filtering
print(paste("Number of genes before filtering:", nrow(dge$counts)))

# Filter lowly expressed genes with a less stringent threshold
keep <- rowSums(cpm(dge) > 0.5) >= 1  # At least 1 sample should have CPM > 0.5
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Check dimensions after filtering
print("Dimensions after filtering:")
print(dim(dge$counts))

# Normalize the data
dge <- calcNormFactors(dge)

# Manually set the common dispersion for data without replicates
common_dispersion <- 0.1  # You can adjust this value based on prior knowledge or literature
dge$common.dispersion <- common_dispersion

# Perform exact test using the manually set common dispersion
et <- exactTest(dge, dispersion = common_dispersion)

# Obtain results
res <- topTags(et, n = Inf)$table

# Print summary of results
print(summary(decideTests(et)))

# Save the results to a CSV file
write.csv(res, file = "edger_results.csv")

# Read the diff.txt file for heatmap generation
diff_data <- read.table("diff.txt", header = TRUE, sep = ",", row.names = 1)

# Print the data to verify
print("Data from diff.txt:")
print(diff_data)

# Ensure diff_data has the correct column names
colnames(diff_data) <- c("logFC")

# Plot the heatmap using ggplot2
ggplot(diff_data, aes(x = rownames(diff_data), y = logFC, fill = logFC)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Heatmap of Differential Expression", x = "Gene", y = "logFC") +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("diff_heatmap.png", width = 6, height = 4)