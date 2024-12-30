# Load required libraries
library(ggplot2)      # For data visualization
library(reshape2)     # For reshaping data (e.g., melt function)
library(grid)         # For adding custom grid elements (e.g., text on plot)
library(tidyverse)    # For data manipulation and visualization
library(psych)        # For statistical analysis
library(pheatmap)     # For generating heatmaps

# Import data files
env <- read.csv("CGenv.csv")  # Load environmental data into 'env'
sp <- read.csv("CGsp.csv")    # Load species data into 'sp'

# Initialize two matrices: one for correlation values and the other for significance markers
correlation_matrix <- matrix(NA, nrow = ncol(env), ncol = ncol(sp))  # Empty matrix for correlation coefficients
significance_matrix <- matrix("", nrow = ncol(env), ncol = ncol(sp)) # Empty matrix for storing significance levels
rownames(correlation_matrix) <- colnames(env)  # Assign row names from environmental variable names
colnames(correlation_matrix) <- colnames(sp)   # Assign column names from species variable names

# Calculate correlations and significance, and add custom significance markers
for (i in 1:ncol(env)) {  # Loop through columns of 'env'
  for (j in 1:ncol(sp)) {  # Loop through columns of 'sp'
    result <- cor.test(env[, i], sp[, j], method = "spearman", use = "complete.obs") # Perform Pearson correlation test
    correlation_matrix[i, j] <- result$estimate  # Store correlation coefficient
    p_value <- result$p.value  # Extract p-value
    # Assign significance markers based on p-value thresholds
    if (p_value < 0.001) {
      significance_matrix[i, j] <- "***"
    } else if (p_value < 0.01) {
      significance_matrix[i, j] <- "**"
    } else if (p_value < 0.05) {
      significance_matrix[i, j] <- "*"
    }
  }
}

# Reshape the correlation matrix into a long format
melted_correlation_matrix <- melt(correlation_matrix)  # Reshape correlation matrix for plotting
colnames(melted_correlation_matrix) <- c("Var1", "Var2", "Correlation")  # Rename columns
melted_correlation_matrix$Correlation <- as.numeric(as.character(melted_correlation_matrix$Correlation))  # Convert correlation values to numeric

# Reshape the significance matrix into a long format
melted_significance_matrix <- melt(significance_matrix)  # Reshape significance matrix for plotting
colnames(melted_significance_matrix) <- c("Var1", "Var2", "Significance")  # Rename columns
melted_significance_matrix$Significance <- as.character(melted_significance_matrix$Significance)  # Ensure significance values are characters

# Ensure the column names of both data frames are consistent
colnames(melted_correlation_matrix) <- c("Var1", "Var2", "Correlation")  # Confirm column names for correlation data
colnames(melted_significance_matrix) <- c("Var1", "Var2", "Significance")  # Confirm column names for significance data

# Convert variable names (Var1, Var2) to characters for compatibility during merging
melted_correlation_matrix[, 1:2] <- lapply(melted_correlation_matrix[, 1:2], as.character)  # Convert Var1, Var2 to character
melted_significance_matrix[, 1:2] <- lapply(melted_significance_matrix[, 1:2], as.character)  # Convert Var1, Var2 to character

# Merge correlation and significance data frames
merged_data <- cbind(melted_correlation_matrix, melted_significance_matrix)  # Combine both data frames
merged_data <- merged_data[, -c(4, 5)] %>% na.omit()  # Remove redundant columns and handle missing values

# Plot heatmap with space for custom text
heatmap_plot <- ggplot(merged_data, aes(Var1, Var2, fill = Correlation)) +  # Create heatmap with correlation values
  geom_tile() +  # Add tiles to represent correlations
  geom_text(aes(label = Significance), color = "white", size = 12) +  # Add significance markers as text
  scale_fill_gradient2(low = "#367DB0", high = "#3D9F3C", mid = "white", midpoint = 0) +  # Customize color gradient
  theme_minimal() +  # Use minimal theme for simplicity
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  # Rotate x-axis text
        axis.text.y = element_text(hjust = 1, size = 20)) +  # Customize y-axis text
  labs(title = "Heatmap of Correlation",  # Add title
       x = "",  # Remove x-axis label
       y = "")  # Remove y-axis label

# Print the heatmap plot
print(heatmap_plot)

# Add custom significance legend text to the plot
grid.text(" *: P < 0.05\n **: P < 0.01\n ***: P < 0.001",  # Define text for significance levels
          x = 0.88, y = 0.4, just = "left", gp = gpar(col = "black", fontsize = 10))  # Place text in a specific position
