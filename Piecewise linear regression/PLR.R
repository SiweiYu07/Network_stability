# Load the required libraries
library(SiZer) # For piecewise linear regression analysis
library(ggplot2) # For data visualization

# Read the dataset containing CGEnvYear and CGEnvPCoA values
seg <- read.csv("CGEnvPCoA.csv")

# Plot the relationship between CGEnvYear and CGEnvPCoA with a smoothed trend line
ggplot(seg, aes(CGEnvYear, CGEnvPCoA)) +
  geom_point(size = 10, color = "#33CCCC") + # Add points to the scatterplot
  geom_smooth(size = 3, color = "blue", fill = "#FF9999") + # Add a smooth fit curve
  xlab("CGEnvYear") + # Label x-axis
  ylab("CGEnvPCoA") + # Label y-axis
  theme(axis.text.x = element_text(size = 27), # Customize x-axis text size
        axis.text.y = element_text(size = 27), # Customize y-axis text size
        axis.title.x = element_text(size = 27), # Customize x-axis title size
        axis.title.y = element_text(size = 27), # Customize y-axis title size
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) # Add a border around the plot

# Perform piecewise linear regression to find a threshold time (breakpoint)
# The model will automatically identify the breakpoint and calculate a 95% confidence interval
model <- piecewise.linear(x = seg$CGEnvYear, y = seg$CGEnvPCoA, CI = TRUE, bootstrap.samples = 1000, sig.level = 0.05)

# Display the results of the piecewise linear regression model
model

# Plot the piecewise linear regression results
plot(model, xlab = 'CGEnvYear', ylab = 'CGEnvPCoA', pch = 19, col = "#33CCCC", cex = 2, # Customize point size, color, etc.
     cex.axis = 2, cex.lab = 2, cex.sub = 2, cex.main = 2, lwd = 4) # Customize axis and line widths
lines(model, col = "blue", lwd = 4) # Add the fitted piecewise linear regression line

# Calculate R-squared to evaluate the goodness of fit
fit1 <- predict(model, seg$CGEnvYear) # Predicted values from the piecewise model
SSre <- sum((seg$CGEnvPCoA - fit1)^2) # Residual sum of squares
SStot <- sum((seg$CGEnvPCoA - mean(seg$CGEnvPCoA))^2) # Total sum of squares
R2 <- round(1 - SSre / SStot, 3) # Compute R-squared
R2

# Load the segmented package for alternative breakpoint analysis
library(segmented)

# Fit a simple linear regression model for the relationship
fit_lm <- lm(CGEnvYear ~ CGEnvPCoA, data = seg) # Linear regression
summary(fit_lm) # Summarize the linear model

# Identify possible breakpoints using the segmented() function
# First method: Specify the number of breakpoints (e.g., 1 breakpoint)
lm_seg1 <- segmented(fit_lm, seg.Z = ~CGEnvPCoA, npsi = 1) # Automatically search for 1 breakpoint
summary(lm_seg1) # Summarize the segmented model
plot(lm_seg1, xlab = 'CGEnvYear', ylab = 'CGEnvPCoA', col = "blue", lwd = 8) # Plot the segmented regression results
points(CGEnvPCoA ~ CGEnvYear, data = seg, pch = 19, col = "#33CCCC", cex = 2, # Add data points to the plot
       cex.axis = 2, cex.lab = 2, cex.sub = 2, cex.main = 2, lwd = 20)

# Second method: Specify an approximate initial breakpoint position
# For example, if 2012 is a possible breakpoint, set it as an initial value
lm_seg2 <- segmented(fit_lm, seg.Z = ~CGZPCoA, psi = 2012) # Refine search near 2012
summary(lm_seg2) # Summarize the segmented model

# Plot the segmented regression results using the second method
plot(lm_seg2, xlab = 'CGZYear', ylab = 'CGZPCoA') # Plot the results
points(CGZYear ~ CGZPCoA, data = seg) # Add points to the plot

# Note: Results may vary between different methods and packages due to differences in algorithms.
