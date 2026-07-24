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
library(Seurat)       # 单细胞分析
library(dplyr)        # 数据操作
library(ggplot2)      # 数据可视化
library(patchwork)    # 图形排版
library(cowplot)      # 图形增
library(dce)

miRNA_filtered <- t(miRNA_log)
mRNA_filtered <- t(mRNA_log)

################################################################################


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

Averg_Duplicate <- function(Exp_scRNA){
  
  # 获取唯一的基因名列表
  uniqueNameList <- unique(colnames(Exp_scRNA))
  noOfgenes <- length(uniqueNameList)
  
  # 初始化结果矩阵
  temp <- matrix(0, nrow = nrow(Exp_scRNA), ncol = noOfgenes)
  colnames(temp) <- uniqueNameList
  rownames(temp) <- rownames(Exp_scRNA)
  
  # 对每个唯一基因计算平均表达值
  for(c in 1:noOfgenes){
    GeneList <- which(colnames(Exp_scRNA) == colnames(temp)[c])
    for(r in 1:nrow(temp)) {
      temp[r, c] <- mean(as.numeric(Exp_scRNA[r, GeneList]))  
    }
  }
  return(temp)
}
# ================== 1. 数据预处理（添加前缀） ==================
# 为表达矩阵列名添加类型前缀
colnames(miRNA_filtered) <- paste0("ceRNA_", colnames(miRNA_filtered))
colnames(mRNA_filtered)  <- paste0("mRNA_",  colnames(mRNA_filtered))
dim(miRNA_filtered)
dim(mRNA_filtered)
miRNA_scRNA_norm_average <- Averg_Duplicate(miRNA_filtered)
mRNA_scRNA_norm_average <- Averg_Duplicate(mRNA_filtered)
dim(miRNA_scRNA_norm_average)
dim(mRNA_scRNA_norm_average)
# 合并并转置（行=基因，列=样本）
fullExpr <- cbind(miRNA_filtered, mRNA_filtered)
fullExpr <- t(fullExpr)   # 若需要转置，请根据实际数据结构调整

# ================== 2. 构建先验网络（含前缀） ==================
# 从fullExpr中获取所有基因名
all_genes <- colnames(fullExpr)

# 创建伪先验网络数据 - 所有可能的miRNA-mRNA组合
# 分离miRNA和mRNA
miRNA_genes <- grep("^ceRNA_", all_genes, value = TRUE)
mRNA_genes <- grep("^mRNA_", all_genes, value = TRUE)

# 创建所有可能的组合
pseudo_prior <- expand.grid(ceRNA = miRNA_genes, mRNA = mRNA_genes)

# 转换为与原始TargetScan相同的格式
pseudo_prior <- data.frame(ceRNA = pseudo_prior$ceRNA, mRNA = pseudo_prior$mRNA)

# 现在使用这个伪先验数据替换原来的TargetScan
TargetScan_clean <- pseudo_prior

# 转换为igraph对象 
priori_graph <- igraph::graph_from_data_frame(
  d = TargetScan_clean,                                                                  
  directed = TRUE,
  vertices = unique(c(TargetScan_clean[,1], TargetScan_clean[,2]))
)

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
set.seed(123)
HCC_timestart <- Sys.time()
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
  num.cores = 40
)

DCEmiR_Null_res <- cellSpecificNets
HCC_timeend <- Sys.time()

save(cellSpecificNets, file = " HCC_70_6898_onlyDCEmiR_Null.RData")


if(length(cellSpecificNets) == nrow(fullExpr)) {
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
