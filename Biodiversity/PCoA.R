#PCOA

library(vegan)      # Load the 'vegan' package for ecological analysis, including PCoA and distance matrix calculations.
library(ggplot2)    # Load the 'ggplot2' package for data visualization.

# Import data, first row.name is sample name, and the first column name is OTU.
otu_raw <- read.table(file="CGEnv.txt", sep="\t", header=TRUE, check.names=FALSE, row.names=1)  
# Read OTU abundance data from a tab-delimited text file. 
# `header=TRUE` indicates the first row contains column names.
# `check.names=FALSE` prevents automatic conversion of column names.
# `row.names=1` sets the first column as row names.

# Matrix transpose
otu <- t(otu_raw)  
# Transpose the OTU matrix so that rows represent samples and columns represent OTUs.

# Calculate Bray-Curtis distance
otu.distance <- vegdist(otu)  
# Compute the Bray-Curtis dissimilarity matrix for the OTU abundance data using the `vegdist` function.

# Perform Principal Coordinates Analysis (PCoA)
pcoa <- cmdscale(otu.distance, eig=TRUE)  
# Perform PCoA on the Bray-Curtis distance matrix.
# `eig=TRUE` extracts eigenvalues for variance explanation.

pc12 <- pcoa$points[,1:2]  
# Extract the first two principal coordinates (PC1 and PC2).

pc <- round(pcoa$eig / sum(pcoa$eig) * 100, digits=2)  
# Calculate the percentage of variance explained by PC1 and PC2, rounded to two decimal places.

# Transform pc12 (matrix) into a data frame
pc12 <- as.data.frame(pc12)  
# Convert the matrix of principal coordinates into a data frame for easier manipulation.

# Add sample names as a new column to the pc12 data frame
pc12$samples <- row.names(pc12)  
# Add a column `samples` containing the row names (sample names).

head(pc12)  
# Display the first six rows of the `pc12` data frame.

# Plot the PCoA result
p <- ggplot(pc12, aes(x=V1, y=V2)) +   
  geom_point(size=5) +            # Add scatter points to the plot with size 5.
  theme_bw()                      # Apply a clean, white-background theme.
p                                 # Display the plot.

# Import group data
group <- read.table("CGEnvgroup.txt", sep='\t', header=TRUE)  
# Read group metadata from a tab-delimited text file.
# `header=TRUE` indicates the first row contains column names.

# Modify column names
colnames(group) <- c("samples", "group")  
# Rename the columns of the group data frame to "samples" and "group".

# Merge sample and group data
df <- merge(pc12, group, by="samples")  
# Merge the PCoA data frame (`pc12`) with the group metadata based on the `samples` column.

head(df)  
# Display the first six rows of the merged data frame.

df  
# Display the entire merged data frame.

# Define custom color palette
color <- c("#3D9F3C", "#367DB0", "#0099FF", "#FF6666")  
# Define a vector of custom colors for the groups.

# Create an enhanced PCoA plot with group information
p1 <- ggplot(data=df, aes(x=V1, y=V2, color=group, shape=group)) +  
  theme_bw() +                                                        # Apply a clean, white-background theme.
  geom_point(size=6) +                                                # Add scatter points to the plot with size 6.
  theme(panel.grid = element_blank()) +                               # Remove grid lines from the plot.
  geom_vline(xintercept = 0, lty="dashed") +                          # Add a dashed vertical line at x=0.
  geom_hline(yintercept = 0, lty="dashed") +                          # Add a dashed horizontal line at y=0.
  geom_text(aes(label=samples, y=V2+0.03, x=V1+0.03, vjust=0), size=6) +  
  # Add sample labels slightly offset from the scatter points for better visibility.
  
  guides(color=guide_legend(title=NULL)) +                            # Remove the title from the legend.
  labs(x=paste0("PC1 ", pc[1], "%"), y=paste0("PC2 ", pc[2], "%")) +  # Add axis labels showing the percentage of variance explained.
  stat_ellipse(data=df, geom="polygon", level=0.9,                    # Add confidence ellipses around groups.
               linetype=2, linewidth=1, aes(fill=group),              # Set ellipse style (dashed, fill, etc.).
               alpha=0.2, show.legend=TRUE) +                         # Adjust transparency of ellipses and include legend.
  
  scale_color_manual(values=color) +                                  # Apply the custom color palette for point colors.
  scale_fill_manual(values=color) +                                   # Apply the same custom colors for ellipse fills.
  
  theme(axis.title.x=element_text(size=20, face="bold"),              # Customize x-axis title (size and bold text).
        axis.title.y=element_text(size=20, angle=90, face="bold"),    # Customize y-axis title (size, angle, and bold text).
        axis.text.y=element_text(size=19, face="bold"),               # Customize y-axis text labels (size and bold text).
        axis.text.x=element_text(size=19, face="bold"),               # Customize x-axis text labels (size and bold text).
        panel.grid=element_blank())                                   # Remove panel grid lines from the plot.

p1  
# Display the enhanced PCoA plot with group information.
