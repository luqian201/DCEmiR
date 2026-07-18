
library(pracma)
library(WGCNA)
library(igraph)
library(energy)
library(Hmisc)
library(parmigene)
library(minerva)
library(glmnet)
library(pcalg)
library(doParallel)
library(philentropy)
library(StatMatch)
library(propr)
library(gtools)
library(pbapply)
library(pcaPP)

 


## 函数：计算矩阵的z-score标准化值
matrixzscore <- function(mat){
  mat.mean <- mean(mat[!is.na(mat)])  # 计算均值(忽略NA值)
  mat.sd <- sd(mat[!is.na(mat)])      # 计算标准差
  mat.zscore <- (mat - mat.mean)/mat.sd  # z-score标准化
  return(mat.zscore)
}

## 函数：通过结合基因表达数据和miRNA-靶标相互作用数据查询基因表达数据
querydata <- function(ExpData, miRTarget) {
  
  ExpDataNames <- colnames(ExpData)
  miRTarget <- as.matrix(miRTarget)
  
  miRTargetQuery <- miRTarget[intersect(which(miRTarget[, 1] %in% ExpDataNames),
                                        which(miRTarget[, 2] %in% ExpDataNames)), ] 
  
  ExpDataQuery <- ExpData[, union(which(ExpDataNames %in% unique(miRTargetQuery[, 1])),
                                  which(ExpDataNames %in% unique(miRTargetQuery[, 2])))]
  
  return(list(ExpDataQuery, miRTargetQuery))
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
## Load R packages
library(dce)
library(igraph)
library(doParallel)
library(WGCNA)
library(pracma)
library(glmnet)
library(bnlearn)
library(pcalg)
library(parallel)
library(ParallelPC)
library(InvariantCausalPrediction)
library(infotheo)

miRNA_filtered <- t(miRNA_log)
mRNA_filtered <- t(mRNA_log)

# ================== 1. 数据预处理（添加前缀） ==================
# 为表达矩阵列名添加类型前缀
colnames(miRNA_filtered) <- paste0("ceRNA_", colnames(miRNA_filtered))
colnames(mRNA_filtered)  <- paste0("mRNA_",  colnames(mRNA_filtered))

dim(miRNA_filtered)
dim(mRNA_filtered)

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
miRNA_scRNA_norm_average <- Averg_Duplicate(miRNA_filtered)
mRNA_scRNA_norm_average <- Averg_Duplicate(mRNA_filtered)

dim(miRNA_scRNA_norm_average)
dim(mRNA_scRNA_norm_average)
# 合并并转置（行=基因，列=样本）
fullExpr <- cbind(miRNA_scRNA_norm_average, mRNA_scRNA_norm_average)
#fullExpr <- t(fullExpr)   # 若需要转置，请根据实际数据结构调整

# ================== 2. 构建先验网络（含前缀） ==================
# 读取 TargetScan 数据
TargetScan <- read.csv("~/LQ/TargetScan_all_9622805.csv", stringsAsFactors = FALSE)
cat("原始数据行数:", nrow(TargetScan), "\n")

# 去除第二列（mRNA）的版本号后缀
TargetScan[, 2] <- gsub("\\.[0-9]+$", "", TargetScan[, 2])
# 若第一列也有版本号，可取消注释：
# TargetScan[, 1] <- gsub("\\.[0-9]+$", "", TargetScan[, 1])

# 基于两列去重（保留唯一调控关系）
TargetScan_unique <- unique(TargetScan)
cat("去重后行数:", nrow(TargetScan_unique), "\n")
cat("去除重复行数:", nrow(TargetScan) - nrow(TargetScan_unique), "\n")

# 为两列添加前缀
TargetScan_unique[, 1] <- paste0("ceRNA_", TargetScan_unique[, 1])
TargetScan_unique[, 2] <- paste0("mRNA_", TargetScan_unique[, 2])

# 构建有向网络
priori_graph <- igraph::graph_from_data_frame(
  d = TargetScan_unique,
  directed = TRUE,
  vertices = unique(c(TargetScan_unique[, 1], TargetScan_unique[, 2]))
)
cat("先验网络节点数:", igraph::vcount(priori_graph), "\n")
cat("先验网络边数:", igraph::ecount(priori_graph), "\n")

# ================== 3. 匹配共有基因（带前缀） ==================
net_nodes <- igraph::V(priori_graph)$name
expr_genes <- colnames(fullExpr)
common_genes <- intersect(net_nodes, expr_genes)

if (length(common_genes) == 0) {
  stop("错误：网络节点与表达矩阵基因无交集，请检查命名一致性。")
}
cat(sprintf("\n匹配统计：\n"))
cat(sprintf("  网络节点总数     : %d\n", length(net_nodes)))
cat(sprintf("  表达矩阵基因数   : %d\n", length(expr_genes)))
cat(sprintf("  共有基因数       : %d\n", length(common_genes)))

# ================== 4. 过滤（仅保留共有基因） ==================
priori_graph <- igraph::induced_subgraph(priori_graph, common_genes)
fullExpr     <- fullExpr[, common_genes, drop = FALSE]   # 保持矩阵结构

# ================== 5. 统计（带前缀） ==================
num_miRNA <- sum(grepl("^ceRNA_", colnames(fullExpr)))
num_mRNA  <- sum(grepl("^mRNA_", colnames(fullExpr)))
cat("\n=== 数据过滤完成（含前缀） ===\n")
cat(sprintf("  最终网络节点数   : %d\n", igraph::vcount(priori_graph)))
cat(sprintf("  最终表达矩阵基因 : %d\n", ncol(fullExpr)))
cat(sprintf("  其中 miRNA       : %d\n", num_miRNA))
cat(sprintf("  其中 mRNA        : %d\n", num_mRNA))
cat(sprintf("  网络边数          : %d\n", igraph::ecount(priori_graph)))

# ================== 6. 清除所有前缀 ==================
# 去除表达矩阵列名的前缀（ceRNA_ 和 mRNA_）
colnames(fullExpr) <- gsub("^(ceRNA_|mRNA_)", "", colnames(fullExpr))

# 去除网络节点名的前缀
igraph::V(priori_graph)$name <- gsub("^(ceRNA_|mRNA_)", "", igraph::V(priori_graph)$name)

cat("\n=== 前缀已清除 ===\n")
cat(sprintf("  清除后网络节点数 : %d\n", igraph::vcount(priori_graph)))
cat(sprintf("  清除后表达矩阵基因数 : %d\n", ncol(fullExpr)))
cat(sprintf("  网络边数          : %d\n", igraph::ecount(priori_graph)))
cat("\n✓ 所有基因名称已恢复为原始符号（无前缀），可继续后续分析。\n")
# --------------- 运行Redice算法 ---------------
K562_timestart <- Sys.time()
set.seed(123)
cellSpecificNets <- Redice(
  priori_graph = priori_graph,  # 使用转换后的igraph对象
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
DCEmiR_TargetScan_K562 <- cellSpecificNets
K562_timeend <- Sys.time()

save(DCEmiR_TargetScan_K562, file = " K562_96_11190_DCEmiR_prior_9621627.RData")


if(length(cellSpecificNets) == nrow(miRNA_filtered)) {
  cat("成功生成", length(cellSpecificNets), "个细胞特异性网络\n")
  cat("首个网络信息：\n")
  print(cellSpecificNets[[2]])
} else {
  warning("网络生成数量与细胞数不一致")
}
# 统计每个网络的边数量
edge_counts <- sapply(cellSpecificNets, function(net) {
  ifelse(is.null(net), 0, ecount(net))
})


# 找出所有非空网络
non_empty_nets <- which(edge_counts >= 0)

# 打印非空网络ID及其边数量
cat("Non-empty networks (ID: edge count):\n")
print(data.frame(ID = non_empty_nets, Edges = edge_counts[non_empty_nets]))


edge_counts <- sapply(cellSpecificNets,
                      ecount)
print(table(edge_counts))#查看边数分布
