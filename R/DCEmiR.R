# Load the K562 and HCC datasets along with putative miRNA-target interactions.
load("Data/K562_19_single-cell_matched_miR_mR.RData")
source("R/Scan.interp.R")


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
library(Seurat)       
library(dplyr)        
library(ggplot2)      
library(patchwork)    
library(cowplot)     
library(dce)


