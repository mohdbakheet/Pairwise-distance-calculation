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
