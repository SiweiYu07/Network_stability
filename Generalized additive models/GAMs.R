# Load required libraries
library("mgcv") # For generalized additive models (GAMs)
#> Loading required package: nlme
#> This is mgcv 1.8-24. For overview type 'help("mgcv-package")'.
library("scam") # For shape-constrained additive models
library("ggplot2") # For data visualization
library("cowplot") # For combining ggplot objects into a single plot
library("tidyr") # For data manipulation and reshaping

# Install and load the 'gratia' package (for model diagnostics and visualization)
install.packages("devtools") # Installs the 'devtools' package for GitHub package installations
devtools::install_github("gavinsimpson/gratia") # Install the 'gratia' package from GitHub
library("gratia") # Load 'gratia' for GAM diagnostics and visualization
## Note: Ensure the package installation works; replace with CRAN alternatives if needed.

# Set default ggplot theme
theme_set(theme_bw()) # Set a black-and-white theme for consistent plotting

# Load source data for analysis
Dc <- read.csv("CGGAM.csv", sep = ',', header = TRUE) # Read CSV file containing GAM data
head(Dc) # Display the first few rows of the data

# Load additional data set for braya analysis
TPTurnover <- read.csv("TPTurnover.csv", sep = ',', header = TRUE) # Read CSV file containing turnover data

# Clean up variable names for the turnover dataset
names(TNTurnover) <- c("Depth", "DepthUpper", "DepthLower", "Year", "YearYoung", # Rename columns for clarity
                       "YearOld", "Turnover")

# Add a variable for the time span of each sediment sample
braya <- transform(braya, sampleInterval = YearYoung - YearOld) # Calculate the interval and add it as a new column
head(braya) # Display the first few rows of the modified dataset

# GAM linear fit analysis
small <- read.csv("CGGAM.csv") # Read the GAM data file
Complexity <- expression(Complexity) # Define the variable 'Stability' for plotting
small_plt <- ggplot(small, aes(x = Temperature, y = Complexity)) + # Create a scatterplot for 'Link' vs 'Stability'
  geom_point(size = 10) + # Add points to the plot with size 10
  labs(y = Complexity, x = "Link") # Add axis labels
plot(small_plt, ncol = 1, labels = "auto", align = "hv", axis = "lr") # Arrange the plot layout

# Fit a GAM model
m <- gam(Complexity ~ s(Temperature, k = 10), data = small, method = "REML") # Fit a GAM with cubic splines and REML estimation
mod <- gamm(Complexity ~ s(Temperature, k = 15), data = small, # Fit a GAMM model with correlation structure
            correlation = corCAR1(form = ~ Depth), method = "REML")
smallPhi <- intervals(mod$lme, which = "var-cov")$corStruct # Extract CAR1 correlation structure
summary(mod$gam) # Summarize the GAM model

# Generate a sequence of values for the CAR1 correlation structure
S <- seq(0, 50, length = 150) # Generate a sequence for time lag (S)
car1 <- setNames(as.data.frame(t(outer(smallPhi, S, FUN = `^`)[1, , ])), # Compute correlation decay
                 c("Lower", "Correlation", "Upper"))
car1 <- transform(car1, S = S) # Add 'S' values to the dataframe
car1Plt <- ggplot(car1, aes(x = S, y = Correlation)) + # Plot the CAR1 correlation structure
  geom_ribbon(aes(ymax = Upper, ymin = Lower), fill = "black", alpha = 0.2) + # Add uncertainty ribbon
  geom_line() # Add the line for the correlation decay

# Predict values for new data points and compute confidence intervals
newYear <- with(small, data.frame(Temperature = seq(min(Temperature), max(Temperature), length.out = 200))) # Generate prediction data
newYear <- cbind(newYear, data.frame(predict(mod$gam, newYear, se.fit = TRUE))) # Predict values and standard errors
crit.t <- qt(0.975, df = df.residual(mod$gam)) # Compute the critical t-value for 95% confidence intervals
newYear <- transform(newYear, upper = fit + (crit.t * se.fit), lower = fit - (crit.t * se.fit)) # Compute confidence bounds

# Plot the fitted GAM model with confidence intervals
sfittt <- ggplot(newYear, aes(x = Temperature, y = fit)) + # Base plot with predicted values
  geom_ribbon(aes(ymin = lower, ymax = upper, x = Temperature), alpha = 0.4, inherit.aes = FALSE, fill = "#66CCCC") + # Add confidence intervals
  geom_point(data = small, size = 7, color = "#3366CC", mapping = aes(x = Temperature, y = Complexity), inherit.aes = FALSE) + # Add data points
  geom_line(size = 2, color = "#FF3366") + # Add the fitted GAM line
  labs(y = Complexity, x = Temperature) + # Add axis labels
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, size = 15), # Customize x-axis text
        axis.text.y = element_text(size = 15)) # Customize y-axis text
sfittt # Display the final plot
