### Network models construction

# Load required libraries for network analysis
library(WGCNA)              # For correlation and network analysis.
library(dynamicTreeCut)     # For hierarchical clustering and dynamic tree cutting.
library(fastcluster)        # Faster hierarchical clustering algorithms.
library(psych)              # For statistical and psychometric functions.
library(reshape2)           # For reshaping data frames.
library(igraph)             # For network creation and analysis.

# Load the OTU (Operational Taxonomic Unit) abundance table
otu_table = read.csv('otu1938-1995.csv', header = TRUE, row.names = 1)  
# Reads a CSV file with OTU data where the first column is used as row names.

# Calculate relative abundances for each OTU
rel_abundance = apply(otu_table, 2, function(x) x / sum(x))  
# Converts raw counts into relative abundances by dividing each value by its column sum.

rel_abundance  # Display the relative abundance matrix.

# Calculate the mean relative abundance for each OTU
mean_rel_abundance = rowMeans(rel_abundance)  
# Computes the average relative abundance across all samples for each OTU.

# Identify OTUs with very low mean relative abundance (< 0.0001)
low_rel_abundance_otu = rownames(otu_table)[mean_rel_abundance < 0.0001]

# Filter out low-abundance OTUs from the OTU table
otu_table_filtered = otu_table[!(rownames(otu_table) %in% low_rel_abundance_otu), ]

# Calculate the frequency of occurrence for each OTU across samples
freq = apply(otu_table_filtered, 1, function(x) sum(x > 0) / length(x))

# Retain only OTUs present in at least 1/6 of the samples
keep = freq >= 1 / 6  
otu = otu_table_filtered[keep, ]  # Filter OTUs based on frequency.

otu = otu_table  # (Optional) Reset to original table, potentially for debugging.

# Calculate Spearman correlation coefficients and p-values between OTUs
cor = corAndPvalue(t(otu_table), y = NULL, use = "pairwise.complete.obs", 
                   alternative = 'two.sided', method = 'spearman')

r = cor$cor  # Extract the correlation coefficients.
p = cor$p    # Extract the p-values.

# Adjust the p-values using the Benjamini-Hochberg method for multiple comparisons
p = p.adjust(p, method = 'BH')

# Filter correlations: Set correlations to 0 if p > 0.001 or |r| < 0.60
r[p > 0.001 | abs(r) < 0.60] = 0  

# Save the filtered correlation matrix as a CSV file
write.csv(data.frame(r, check.names = FALSE), 'corr.matrix1938-1995.csv')  

# Create a weighted undirected graph from the correlation matrix
g = graph_from_adjacency_matrix(r, mode = "undirected", weighted = TRUE, diag = FALSE)

# Remove isolated nodes (degree = 0) from the graph
g = delete_vertices(g, names(degree(g)[degree(g) == 0]))

# Assign correlation values to the 'corr' edge attribute
E(g)$corr = E(g)$weight  

# Set the edge weights to the absolute value of the correlations
E(g)$weight = abs(E(g)$weight)

# Load taxonomic classification data for nodes
tax = read.csv('taxa1938-1995.csv', row.names = 1, header = TRUE)  

# Match the taxonomy data to the nodes in the graph
tax = tax[as.character(V(g)$name), ]  

# Assign taxonomic attributes to the graph nodes
V(g)$Kingdom = tax$Kingdom  
V(g)$Phylum = tax$Phylum  
V(g)$Class = tax$Class  
V(g)$Order = tax$Order  
V(g)$Family = tax$Family  
V(g)$Genus = tax$Genus  
V(g)$Species = tax$Species  

# Create a node list with taxonomic information
node_list = data.frame(
  label = names(V(g)), 
  kingdom = V(g)$Kingdom, 
  phylum = V(g)$Phylum, 
  class = V(g)$Class, 
  order = V(g)$Order, 
  family = V(g)$Family, 
  genus = V(g)$Genus, 
  species = V(g)$Species)

head(node_list)  # Display the first few rows of the node list.

# Save the node list as a CSV file
write.csv(node_list, 'network.node_list1938-1995.csv')

# Create an edge list from the graph
edge = data.frame(as_edgelist(g))  
edge_list = data.frame(
  source = edge[[1]], 
  target = edge[[2]], 
  weight = E(g)$weight, 
  correlation = E(g)$corr)

head(edge_list)  # Display the first few rows of the edge list.

# Save the edge list as a CSV file
write.csv(edge_list, 'network.edge_list1938-1995.csv')

# Export the graph in GraphML format for visualization in Gephi
write_graph(g, 'network1938-1995.graphml', format = 'graphml')

# Calculate basic network properties
nodes_num = length(V(g))  # Number of nodes in the graph.
edges_num = length(E(g))  # Number of edges in the graph.

positive.cor_num = sum(E(g)$corr > 0)  # Number of positive correlations.
negative.cor_num = sum(E(g)$corr < 0)  # Number of negative correlations.

average_degree = mean(degree(g))  # Average degree of nodes.
average_path_length = mean_distance(g, directed = FALSE)  # Average path length.
network_diameter = diameter(g, directed = FALSE)  # Network diameter.
network_density = edge_density(g)  # Network density.
clustering_coefficient = transitivity(g)  # Clustering coefficient.

# Create a summary of the network parameters
network_parameter = data.frame(
  nodes_num, 
  edges_num, 
  positive.cor_num, 
  negative.cor_num, 
  average_degree, 
  average_path_length, 
  network_diameter, 
  network_density, 
  clustering_coefficient)

network_parameter  # Display the network parameters.

# Save the network parameters as a CSV file
write.csv(network_parameter, 'network_parameter1938-1995.csv')

# Create a binary adjacency matrix (presence/absence) from the OTU table
otu1 = otu
otu1[otu1 > 0] = 1

# Save the binary adjacency matrix as a CSV file
write.csv(otu1, 'adjacent_matrix1938-1995.csv')
