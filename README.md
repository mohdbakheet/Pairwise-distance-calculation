### **Pipeline for Genetic Distance and Phylogenetic Tree Analysis**

### **1. Set Up Environment and Load Libraries**

```r
# Load required libraries
library(ape)
library(geiger)
library(phytools)
library(phangorn)
library(adegenet)
library(Biostrings)

# Set working directory to where your data is stored
setwd("your working directory")
```

---

### **2. Load Sequence Data**

#### 2.1 Nucleotide Sequences
```r
# Load nucleotide sequence in phyDat format
alignmentPhyDat <- read.phyDat(file = "sequences.fasta", format = "fasta", type = "DNA")

# Load nucleotide sequences in DNAbin format
nuceseq <- fasta2DNAbin("aligned-seqs.fasta")
```

#### 2.2 Amino Acid Sequences
```r
# Load amino acid sequence in Phylip format
aaseq <- read.aa("result.phylipi", format = "interleaved")

# Convert amino acid sequence to phyDat format for further analysis
aaseq_phyDat <- as.phyDat(aaseq)
```

#### 2.3 Phylogenetic Tree
```r
# Load the phylogenetic tree in Newick format
smithtree <- read.tree("tree.nwk")
```

---

### **3. Pairwise Nucleotide Distance Calculation**

#### 3.1 Using DNAbin Format
```r
# Pairwise nucleotide distances (pure pairwise)
pure_pairwise <- dist.gene(x = nuceseq)
write.csv(as.matrix(pure_pairwise), "pairwise_nuc.csv")

# Logdet distance using the logdet model
logdet_dist <- dist.dna(x = nuceseq, model = "logdet", as.matrix = TRUE)
write.csv(logdet_dist, "nuclogdet.csv")
```

#### 3.2 Using phyDat Format
```r
# p-distance calculation
p_dist <- dist.p(x = alignmentPhyDat)
write.csv(as.matrix(p_dist), "p_dist_nuc.csv")

# Hamming distance calculation
hamming_dist <- dist.hamming(x = alignmentPhyDat)
write.csv(as.matrix(hamming_dist), "hammingdist_nuc.csv")
```

---

### **4. Pairwise Amino Acid Distance Calculation**

#### 4.1 Different Amino Acid Models
```r
models <- c("WAG", "JTT", "LG", "Dayhoff", "Blosum62", "Dayhoff_DCMut", "JTT_DCMut")

# Jukes-Cantor (JC69) model for amino acids
aadist_JC69 <- dist.ml(x = aaseq_phyDat, model = "JC69", exclude = "none")
write.csv(as.matrix(aadist_JC69), "aadis_JC69.csv")

# F81 model for amino acids
aadist_F81 <- dist.ml(x = aaseq_phyDat, model = "F81", exclude = "none")
write.csv(as.matrix(aadist_F81), "aadis_F81.csv")

# FLU Model
aadist_FLU <- dist.ml(x = aaseq_phyDat, model = "FLU", exclude = "none")
write.csv(as.matrix(aadist_FLU), "aadis_FLU.csv")

# WAG Model
aadist_WAG <- dist.ml(x = aaseq_phyDat, model = "WAG", exclude = "none")
write.csv(as.matrix(aadist_WAG), "aadis_WAG.csv")

# JTT Model
aadist_JTT <- dist.ml(x = aaseq_phyDat, model = "JTT", exclude = "none")
write.csv(as.matrix(aadist_JTT), "aadis_JTT.csv")

# LG Model
aadist_LG <- dist.ml(x = aaseq_phyDat, model = "LG", exclude = "none")
write.csv(as.matrix(aadist_LG), "aadis_LG.csv")

# Dayhoff Model
aadist_Dayhoff <- dist.ml(x = aaseq_phyDat, model = "Dayhoff", exclude = "none")
write.csv(as.matrix(aadist_Dayhoff), "aadis_Dayhoff.csv")
```

#### 4.2 LogDet Distance (Amino Acids)
```r
# LogDet distance for amino acids (Hadamard conjugation)
aalogdet <- dist.logDet(x = aaseq_phyDat)
write.csv(as.matrix(aalogdet), "aalogdet.csv")
```

---

### **5. Phylogenetic Tree-Based Distance Calculation**

```r
# Cophenetic distance from the phylogenetic tree (tree-based distances)
cophenetictree <- cophenetic.phylo(smithtree)
write.csv(as.matrix(cophenetictree), "cophenetic_tree.csv")
```

---

### **6. IQTree Analysis for Codon and Protein Models**

#### 6.1 Prepare the Input
Once your sequence files and models are ready, you can move to IQTree for more complex analysis.

1. Use the `IQTree` command-line tool to run codon and protein-mixing models:
   - Example command:
   ```bash
   iqtree -s sequences.fasta -m GTR+I+G -bb 1000 -alrt 1000 -st AA
   ```
   This example runs a GTR+I+G model on your sequences with amino acids (`-st AA`) and bootstraps (`-bb 1000`).

2. Explore other codon and protein substitution models suitable for your analysis. IQTree supports various complex models beyond the basic ones used above.

---

### **7. Visualization and Analysis**

- Use `ape` and `phytools` for tree visualization.
- You can also use `ggtree` or `ggplot2` for more customized and publication-quality plots.

---

###  **K-Means Clustering and Saving Cluster Members to CSV**

---

### **1. Load Required Libraries**

```r
# Load the required libraries
library(tidyverse)  # For data manipulation and visualization
library(cluster)    # For clustering algorithms
library(factoextra)  # For clustering visualization and finding the optimal number of clusters
library(readr)      # For reading CSV files
```

---

### **2. Read Data and Handle Missing Values**

```r
# Set working directory
setwd("your working directory")

# Read in the data
data <- read_csv("aalogdet_coord.csv")

# Check the structure and head of the data to understand its content
str(data)
head(data)

# Handle missing values by removing rows with NA (alternative: consider imputation)
data <- na.omit(data)  # This removes rows with any NA values
```

---

### **3. Feature Selection and Scaling**

```r
# Select the relevant features for clustering (dim.1 and dim.2)
features <- select(data, dim.1, dim.2)

# Scale the selected features to standardize the values for clustering
scaled_features <- scale(features)
```

---

### **4. Determining the Optimal Number of Clusters Using Elbow Method**

```r
# Visualize the Elbow Method to find the optimal number of clusters
fviz_nbclust(scaled_features, kmeans, method = "wss") +
  geom_vline(xintercept = 12, linetype = 2) +  # Adjust based on observation
  labs(subtitle = "Elbow Method for Optimal Clusters")
```
- **Note:** Inspect the plot and adjust the `geom_vline(xintercept)` based on where the "elbow" is located, i.e., where the reduction in WSS starts to slow down.
  
---

### **5. K-Means Clustering**

```r
# Set a seed for reproducibility
set.seed(42)

# Determine the optimal number of clusters from the Elbow Method visualization
optimal_clusters <- 11  # Adjust this based on the Elbow Method result

# Perform K-Means clustering
kmeans_result <- kmeans(scaled_features, centers = optimal_clusters, nstart = 25)

# Append cluster labels to the original data
data$Cluster <- kmeans_result$cluster
```

---

### **6. Cluster Visualization**

```r
# Visualize the clustering results
fviz_cluster(kmeans_result, data = scaled_features, geom = "point", stand = FALSE) +
  ggtitle("K-Means Cluster Visualization") +
  theme_minimal() +  # Apply a cleaner theme
  scale_color_brewer(palette = "Set2")  # Optional: improve visualization with color palettes
```

- **Improvement:** I added `theme_minimal()` and a color palette `scale_color_brewer()` to make the plot visually clearer.

---

### **7. Centroid Distance Calculation**

```r
# Extract the centroids from the K-Means result
centroids <- kmeans_result$centers

# Calculate the pairwise distances between centroids
centroid_distances <- dist(centroids)

# Convert the distance object to a matrix for easier interpretation
centroid_distance_matrix <- as.matrix(centroid_distances)

# Print the centroid distance matrix
print(centroid_distance_matrix)
```

---

### **8. Save Cluster Members to Separate CSV Files**

```r
# Check if the directory for saving cluster files exists, create if not
if (!dir.exists("new_clusters_files2")) {
  dir.create("new_clusters_files2")
}

# Loop through each cluster and save members to separate CSV files
for (i in 1:optimal_clusters) {
  cluster_members <- data[data$Cluster == i, ]
  write_csv(cluster_members, paste0("new_clusters_files2/cluster_", i, "_members.csv"))
}
```

---
