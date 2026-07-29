# 🧬 DCEmiR

**DCEmiR: A Method for Identifying Cell-Specific microRNA Causal Regulatory Networks**

# 📰 Background

In the gene regulation field, miRNA regulation has attracted broad attention due to its potential for clinical translation. However, cancer cells exhibit significant heterogeneity, and miRNA regulatory functions can differ markedly across cell types or even within the same cell under varying conditions. To decipher the molecular mechanisms underlying this heterogeneity, researchers have developed numerous methods for constructing gene regulatory networks (GRNs) using bulk transcriptomic data. While these approaches are valuable, they only capture average regulatory levels across cell populations, failing to reveal cell-specific regulatory heterogeneity. Emerging single-cell methods have begun to address this gap, yet they remain limited by their reliance on correlation-based analyses and lack of causal inference capabilities, with the inherent noise and sparsity of single-cell data further compromising network accuracy.

To overcome these limitations, we propose DCEmiR, a novel method for identifying cell-specific miRNA-mRNA causal regulatory networks at single-cell resolution. By integrating single-cell transcriptomic data with prior knowledge of miRNA-target interactions and accounting for inter-cellular heterogeneity, DCEmiR constructs personalized causal regulatory networks at the individual cell level. We applied this method to single-cell datasets from hepatocellular carcinoma and leukemia, systematically characterizing miRNA regulatory networks and revealing key causal pathways across cancer cells. Our findings provide theoretical and methodological support for personalized cancer diagnosis and therapy, with potential applications in early detection, prognosis, drug target discovery, and regenerative medicine.

A schematic illustration of **DCEmiR** is shown in the folowing. ![A schematic illustration of DCEmiR](/DCEmiR_schematic_illustration.png) For single-cell transcriptomic data, whether or not they contain prior information on miRNA–mRNA interactions, DCEmiR first screens genes based on existing prior knowledge (where available). Subsequently, using a perturbation strategy, DCEmiR removes cells one by one to construct a control group (background data) and n experimental groups (perturbed data). These n+1 datasets are then processed using a differential causal effects model to generate n+1 causal matrices. By subtracting the corresponding coefficients from the background data from the causal effect coefficients of each perturbation dataset, DCEmiR obtains a differential causal matrix for each cell and visualises it as n directed single-cell gene regulatory network topologies for individual cells.

# 📁 Description of each file in R and Data folders

-   **K562_19_single-cell_matched_miR_mR.RData:** miRNA and mRNA expression data from 19 single leukaemia (K562) cells

-   **HCC_32_single-cell_matched_miR_mR.RData:** miRNA and mRNA expression data from 32 single cells of hepatocellular carcinoma (HCC)

-   **TargetScan_9621627_v8.0.csv:** Predicted miRNA–target gene interactions (a priori data)

-   **miRTarBase_v10.0+TarBase_v9.0.csv:** Gold standard validation data

-   **DCEmiR_HCC.R:** Code for using DECmiR in HCC

-   **DCEmiR_K562.R:** Code for using DECmiR in K562

# 📌 The usage of DCEmiR

Paste the miRNA and mRNA expression data, along with the predicted miRNA–target gene interactions (optional), into a folder (set this folder as the directory for the R environment).

# 🚀 Quick example to use DCEmiR

To identify cell-specific microRNA causal regulatory networks, users need to prepare matched miRNA and mRNA expression data, along with putative miRNA-target interactions (optional). Place the datasets and our source script (DCEmiR_HCC.R and DCEmiR_K562.R) into a single folder, and set this folder as the R working directory. Then, run the following scripts to perform the analysis. For convenience, we have prepared single-cell transcriptomics datasets and putative miRNA-target interactions for users, which can be downloaded [here](Data/).

``` r
# DCEmiR application in HCC dataset:
# Load the HCC datasets.
load("Data/HCC_32_single-cell_matched_miR_mR.RData")

# Load required packages
library(igraph)
library(doParallel)
library(dce)

## Function definitions ##
Averg_Duplicate <- function(Exp_scRNA){
  uniqueNameList <- unique(colnames(Exp_scRNA))
  noOfgenes <- length(uniqueNameList)
  temp <- matrix(0, nrow = nrow(Exp_scRNA), ncol = noOfgenes)
  colnames(temp) <- uniqueNameList
  rownames(temp) <- rownames(Exp_scRNA)
  for(c in 1:noOfgenes){
    GeneList <- which(colnames(Exp_scRNA) == colnames(temp)[c])
    for(r in 1:nrow(temp)) {
      temp[r, c] <- mean(as.numeric(Exp_scRNA[r, GeneList]))  
    }
  }
  return(temp)
}

Redice <- function(priori_graph, 
                   ExpData, 
                   adjustment_type = "parents", 
                   effect_type = "total", 
                   p_method = "hmp", 
                   test = "wald", 
                   deconfounding = TRUE, 
                   conservative = FALSE,
                   p.value = 0.05,
                   num.cores = 1){
  
  cl <- makeCluster(num.cores)
  registerDoParallel(cl) 
  res.single.cell.pvalue <- foreach(i = seq(nrow(ExpData)), .packages = c("dce")) %dopar% {
    dce(priori_graph, ExpData, ExpData[-i, ], adjustment_type = adjustment_type, effect_type = effect_type, p_method = p_method, test = test, deconfounding = deconfounding, conservative = conservative)$dce_pvalue
  }
  stopCluster(cl)
  stopImplicitCluster() 
  res.single.cell <- lapply(seq(res.single.cell.pvalue), function (i) abs(res.single.cell.pvalue[[i]] < p.value))
  for (i in seq(res.single.cell)) {
    res.single.cell[[i]][is.na(res.single.cell[[i]])] <- 0
  }
  res.single.cell.graph <- lapply(seq(res.single.cell), function (i) graph_from_adjacency_matrix(res.single.cell[[i]]))
  return(res.single.cell.graph)
}

## Data preprocessing ##
# Log transformation
miRNA_raw_norm <- log2(miRNA_raw_count + 1)
mRNA_raw_norm <- log2(mRNA_raw_count + 1)

# Remove genes with zero expression in more than 80% of cells
zero_ratio_miRNA <- rowSums(miRNA_raw_norm == 0) / ncol(miRNA_raw_norm)
miRNA_raw_norm <- miRNA_raw_norm[zero_ratio_miRNA <= 0.8, ]
zero_ratio_mRNA <- rowSums(mRNA_raw_norm == 0) / ncol(mRNA_raw_norm)
mRNA_raw_norm <- mRNA_raw_norm[zero_ratio_mRNA <= 0.8, ]

# Average expression values for duplicate columns (samples)
miRNA_scRNA_norm_average <- Averg_Duplicate(miRNA_raw_norm)
mRNA_scRNA_norm_average <- Averg_Duplicate(mRNA_raw_norm)

# Merge expression matrices (rows = genes, columns = samples)
fullExpr <- t(cbind(t(miRNA_scRNA_norm_average), t(mRNA_scRNA_norm_average)))

## Reduce computational cost using prior network ##
# With priori: construct directed prior network
TargetScan <- read.csv("Data/TargetScan_9621627_v8.0.csv")

# Build directed prior network
priori_graph <- igraph::graph_from_data_frame(
  d = TargetScan,
  directed = TRUE,
  vertices = unique(c(TargetScan[, 1], TargetScan[, 2]))
)

# Filter expression matrix to keep only genes present in the prior network
net_nodes <- igraph::V(priori_graph)$name
expr_genes <- colnames(fullExpr)
common_genes <- intersect(net_nodes, expr_genes)
priori_graph <- igraph::induced_subgraph(priori_graph, common_genes)
fullExpr <- fullExpr[, common_genes, drop = FALSE]   # Preserve matrix structure

# ## Without priori: create pseudo-prior network
# # Extract miRNA and mRNA gene names from raw expression matrices
# miRNA_genes <- rownames(miRNA_scRNA_norm_average)
# mRNA_genes <- rownames(mRNA_scRNA_norm_average)
# # Create all possible miRNA-mRNA combinations as pseudo-prior network
# pseudo_prior <- expand.grid(miRNA = miRNA_genes, mRNA = mRNA_genes)
# # Build directed prior network
# priori_graph <- igraph::graph_from_data_frame(
#   d = pseudo_prior,
#   directed = TRUE,
#   vertices = unique(c(pseudo_prior$miRNA, pseudo_prior$mRNA))
# )

## Construct DCEmiR gene regulatory networks ##
HCC_timestart <- Sys.time()
set.seed(123)
DCEmiR_TargetScan_HCC <- Redice(
  priori_graph = priori_graph,
  ExpData = fullExpr,
  adjustment_type = "parents",
  effect_type = "total",
  p_method = "hmp",
  test = "wald",
  deconfounding = TRUE,
  conservative = FALSE,
  p.value = 0.05,
  num.cores = 70
)
HCC_timeend <- Sys.time()

# Count edges in each network
edge_counts <- sapply(DCEmiR_TargetScan_HCC, function(net) {
  ifelse(is.null(net), 0, ecount(net))
})

# Identify non-empty networks
non_empty_nets <- which(edge_counts >= 0)

# Print non-empty network IDs and their edge counts
cat("Non-empty networks (ID: edge count):\n")
print(data.frame(ID = non_empty_nets, Edges = edge_counts[non_empty_nets]))
```
