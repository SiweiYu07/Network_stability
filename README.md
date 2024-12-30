This repository contains code, data, and scripts for measuring temporal biodiversity and constructing network models in ecological studies. The project uses multiple statistical and modeling approaches to investigate temporal changes in biodiversity and the stability of ecological networks.


Project Overview
The primary aim of this project is to evaluate temporal biodiversity patterns and construct ecological network models. The methods applied include:

Biodiversity indices (e.g., alpha and beta diversity).
Correlation analysis between environmental and community variables.
Generalized additive models (GAM) to explore non-linear relationships.
Piecewise linear regressions to identify thresholds or breakpoints.
Partial Least Squares Path Modeling (PLSPM) to identify causal relationships.
Ecological network construction to assess community stability and interactions.
Folder Structure
The repository is organized into the following directories:

1. Biodiversity/
Scripts for calculating alpha diversity, beta diversity, and functional diversity indices.
Temporal biodiversity trends analysis.
2. Correlation analysis/
Code for running correlation analyses between environmental variables (e.g., temperature, nutrients) and community metrics (e.g., species richness, turnover).
Includes visualization scripts for heatmaps and correlation matrices.
3. Data/
Contains input datasets used for the analysis, such as:
Community composition data (species/genus-level abundance tables).
Environmental variables (e.g., temperature, pH, nutrients).
Network interaction datasets.
4. Generalized additive models/
Scripts for building and analyzing Generalized Additive Models (GAMs).
Used to explore non-linear relationships between biodiversity metrics and environmental variables over time.
5. Network construction/
Contains scripts for constructing ecological networks based on species co-occurrence and interaction patterns.
Includes methods to evaluate network complexity, degree centrality, and stability metrics.
6. Piecewise linear regression/
Code for identifying ecological thresholds using piecewise linear regression.
Helpful for pinpointing significant regime shifts in biodiversity and environmental conditions.
7. PLSPM/
Scripts for Partial Least Squares Path Modeling (PLSPM).
Evaluates causal relationships between environmental variables, biodiversity metrics, and community stability.

Key Features
Temporal Analysis: Tracks biodiversity changes over time using alpha, beta, and functional diversity indices.
Non-linear Modeling: Applies GAMs to uncover non-linear relationships and temporal trends.
Threshold Identification: Detects critical ecological thresholds using piecewise regressions.
Network Models: Constructs and analyzes ecological networks to understand interactions and stability.
Path Modeling: Uses PLSPM to study causal relationships and identify key drivers of biodiversity patterns.

Dependencies
The following R packages are required to run the analyses:
vegan: For diversity metrics.
plspm: For Partial Least Squares Path Modeling.
mgcv: For Generalized Additive Models.
igraph: For network construction and analysis.
segmented: For piecewise regression models.
corrplot: For visualizing correlation matrices.
