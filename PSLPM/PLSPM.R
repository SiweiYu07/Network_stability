
# Load the plspm package for Partial Least Squares Path Modeling
library(plspm) 

# Read the dataset "CG11.csv" into a data frame
test_otu <- read.csv("CG11.csv") 

# Select rows 1 to 26 and specific columns for the subset dataset
df <- test_otu[1:26, c(1:6, 7:12)]

# Rename the columns for easier interpretation
colnames(df) <- c("Precip", "Temp", "TOC", "TN", "TP", "IP", 
                  "Prichness", "Zrichness", "Pturnover", 
                  "Zturnover", "Degree", "Robustness")

# (Optional) Standardize the data by centering and scaling
# df = scale(df[, 1:12], center = TRUE, scale = TRUE)

# Define the relationships (inner model) between latent variables
Precip <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
Temp <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
TOC <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
TN <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
TP <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
Prichness <- c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
Zrichness <- c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0)
Pdissimilarity <- c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0)
Zdissimilarity <- c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0)
Complexity <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0)
Stability <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0)

# Combine the relationships into a path matrix for the inner model
df_path <- rbind(Precip, Temp, TOC, TN, TP, Prichness, Zrichness,
                 Pdissimilarity, Zdissimilarity, Complexity, Stability)

# Visualize the inner model structure
innerplot(df_path)

# Define the outer model by associating observed variables with latent variables
df_blocks = list(1:1, 2:2, 3:3, 4:4, 5:6, 7:7, 8:8, 9:9, 10:10, 11:11, 12:12)

# Specify measurement mode for each block; "A" indicates a reflective model
df_modes = rep("A", 11)

# Perform the PLSPM analysis with bootstrap validation enabled
df_pls = plspm(df, df_path, df_blocks, modes = df_modes, boot.val = TRUE)

# Plot loadings to evaluate the relationship between latent and observed variables
plot(df_pls, what = "loadings", arr.width = 0.1, show.values = TRUE, lcol = 'gray')

# Visualize the inner model with path coefficients and significance
innerplot(df_pls, colpos = '#CA5023', colneg = '#457CC3', show.values = TRUE, 
          lcol = 'gray20', box.lwd = 0)

# Display path coefficients (inner model) and their significance
df_pls$inner_model

# Compute the goodness-of-fit (GOF) index; GOF > 0.7 indicates a reliable model
df_pls$gof

# Output a detailed summary of the PLSPM results
summary(df_pls)
