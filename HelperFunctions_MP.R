# Helper funtions

## Libraries
# ----------------------------------------------------------------------- #

library(data.table)
library("WGCNA")
library("RColorBrewer")
library("tidyverse")
library("gtools")
library('plyr')
library("ggplot2")
library("reshape2")


# Calculate Counts from log10(x+1)
# ----------------------------------------------------------------------- #

CalcCount <- function(data){
  counts <- apply(data, 2, function(x){
    (10^x)-1
  })
  rownames(counts) <- rownames(data)
}


# Determine beta and plot sft
# ----------------------------------------------------------------------- #

filterData <- function(data){
  multExpr <- data
  # WGCNA use dataframes with samples in rows and genes/probes in columns
  
  # Filtering of NA-Entries in Genes and Samples
  k <- WGCNA::goodSamplesGenes(multExpr, minFraction = 0.1)
  k$allOK # if TRUE the next command will nothing change, else genes and/or samples will be deleted
  multExpr <- multExpr[which(k$goodSamples),which(k$goodGenes)]
  
  #print(dim(multExpr)) 
}


# Determine beta and plot sft
# ----------------------------------------------------------------------- #
FindBetaAndPlot <- function(data, powers, data_set_name){
  
  #powers <- c(c(1:10), seq(from = 12, to=20, by=2))
  sft <- WGCNA::pickSoftThreshold(data, powerVector = powers, verbose = 5, networkType = "signed")
  
  svg(file=paste("R_01_", project, "_softThresholdBeta_", data_set_name ,".svg", sep = ""), height = 5, width = 6)
  par(mfrow=c(1,2))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
       ylim = c(0,1), xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
       type="n",main = paste("Scale independence")) 
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=1.0,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,6], xlab="Soft Threshold (power)",ylab="Median Connectivity", type="n", main = paste("Median connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,6], labels=powers, cex=1.0,col="red")
  dev.off()
  
  # Also plot it in md
  par(mfrow=c(1,2))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
       ylim = c(0,1), xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
       type="n",main = paste("Scale independence")) 
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=0.9,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,6], xlab="Soft Threshold (power)",ylab="Median Connectivity", type="n", main = paste("Median connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,6], labels=powers, cex=0.9,col="red")
}


# Find data files and names
# ----------------------------------------------------------------------- #
FindDataFiles <- function(pattern, what2replace){
  
  files <<- list.files(pattern = pattern)
  names_data <<- gsub(files, pattern = what2replace, replacement = "")
  
}

LoadData <- function(data_path){
  x <- load(data_path)
  y <<- get(x)
}


CalcModulePres <- function(ref_data_path, test_data_path, moduleVector_path, data_name, nPerms, transpose_test = FALSE){
  if (data_name == "Gold"){
    # load gold data and color vector
    Ref_data <- LoadData(ref_data_path)
    Test_data <- LoadData(test_data_path)
    Module_data <- LoadData(moduleVector_path)
    #assign 2 MP framework
    setLabels = c("Ref", "Test");
    multiExpr = list(Ref = list(data = Ref_data), 
                     Test = list(data = Test_data));
    multiColor = list(Ref = Module_data)
    #call the function
    svg(file=paste("R_",Analysis_Part,"_",project,"_ModPres_Gold2",data_name,".svg", sep=""))
    check_preservation(multiExpr=multiExpr, 
                       multiColor=multiColor, 
                       nPerm=nPerms)
    dev.off()
    
    CompData <<- cbind(moduleSizes,mp$preservation$Z[[1]][[2]][, 2])
    
    save(CompData, file = paste0(rdata, "03_",project,"_CompData_Gold_MP.RData"))
  } else if (transpose_test == TRUE){
    Ref_data <- LoadData(ref_data_path)
    Test_data <- LoadData(test_data_path)
    Module_data <- LoadData(moduleVector_path)
    #assign 2 MP framework
    setLabels = c("Ref", "Test");
    multiExpr = list(Ref = list(data = Ref_data), 
                     Test = list(data = t(Test_data)));
    multiColor = list(Ref = Module_data)
    #call the function
    svg(file=paste("R_",Analysis_Part,"_",project,"_ModPres_Gold2",data_name,".svg", sep=""))
    check_preservation(multiExpr=multiExpr, 
                       multiColor=multiColor, 
                       nPerm=nPerms)
    dev.off()
  } else {
    Ref_data <- LoadData(ref_data_path)
    Test_data <- LoadData(test_data_path)
    Module_data <- LoadData(moduleVector_path)
    #assign 2 MP framework
    setLabels = c("Ref", "Test");
    multiExpr = list(Ref = list(data = Ref_data), 
                     Test = list(data = Test_data));
    multiColor = list(Ref = Module_data)
    #call the function
    svg(file=paste("R_",Analysis_Part,"_",project,"_ModPres_Gold2",data_name,".svg", sep=""))
    check_preservation(multiExpr=multiExpr, 
                       multiColor=multiColor, 
                       nPerm=nPerms)
    dev.off()
  } 
}

# Post processing of module preservation comparison plots

PostProcessMP <- function(comp_data, dataset_name, names_vec, module_names){
  rownames(comp_data) <- module_names
  colnames(comp_data) <- c("Module Size", "Gold", names_vec)
  save(comp_data, file = paste0(Analysis_Part, "_", project, "_Comp_Gold2", dataset_name, ".RData"))
}

# Infer module preservation
#--------------------------------------------------------------------------------------------------------

check_preservation= function(multiExpr, multiColor, nPerm){
  mp <<- modulePreservation(multiExpr, multiColor, referenceNetworks = 1, nPermutations = nPerm, 
                            randomSeed = 1, dataIsExpr =T,
                            maxModuleSize = 5000,
                            quickCor = 0,
                            verbose = 3)
  
  ref = 1
  test = 2
  statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
  statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
  # Module labels and module sizes are also contained in the results
  modColors <<- rownames(mp$preservation$observed[[ref]][[test]])
  moduleSizes <<- mp$preservation$Z[[1]][[test]][, 1];
  # leave grey and gold modules out
  plotMods = !(modColors %in%  "gold");
  # Text labels for points
  text = modColors[plotMods];
  # Auxiliary convenience variable
  plotData <<- cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
  # Main titles for the plot
  mains = c("Preservation Median rank", "Preservation Zsummary");
  #sizeGrWindow(10, 5);
  par(mfrow = c(1,2))
  par(mar = c(4.5,4.5,2.5,1))
  for (p in 1:2)
  {
    min = min(plotData[, p], na.rm = TRUE);
    max = max(plotData[, p], na.rm = TRUE);
    # Adjust ploting ranges appropriately
    if (p==2)
    {
      if (min > -max/10) min = -max/10
      ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
    } else
      ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
    plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
         main = mains[p],
         cex = 2.4,
         ylab = mains[p], xlab = "Module size", log = "x",
         ylim = c(-1,100),
         xlim = c(10, max(moduleSizes[plotMods])), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
    text(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 0.5, offs = 0.08);
    # For Zsummary, add threshold lines
    if (p==2)
    {
      abline(h=0)
      abline(h=2, col = "blue", lty = 2)
      abline(h=10, col = "darkgreen", lty = 2)
    }
  }
}


# Build WGCNA tree
#--------------------------------------------------------------------------------------------------------

build_tree <- function( multiExpr, beta, data_type, minClusterSize){
  
  adjacency <- WGCNA::adjacency(multiExpr, type = "signed" , power = beta, corFnc = "bicor")
  TOM <<- WGCNA::TOMsimilarity(adjacency);
  dissTOM = 1-TOM
  geneTree <- hclust(as.dist(dissTOM),method="average");
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                              method="hybrid", deepSplit = 2, 
                              pamRespectsDendro = FALSE, minClusterSize = minClusterSize);
  dynamicColors <<- labels2colors(dynamicMods)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                      dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                      guideHang = 0.05)
  
  svg(file=paste0("R_02_", project,"_", data_type ,"_Clusterdendro.svg"), width = 8 , height = 5)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                      dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                      guideHang = 0.05, main =paste("Cluster Dendrogram of", data_type))
  dev.off()
}

# Calculate fold change
#--------------------------------------------------------------------------------------------------------


calc_Foldchange <- function(comp){
  remove <- "gold"
  dat <- comp[!rownames(comp) %in% remove, ]
  
  dat <- dat[order(dat[,1], decreasing = T),]
  dat_final <- foldchange(abs(dat[,3:dim(comp)[2]]), dat[,2])
  fc_comparison <<- cbind(dat[,1], dat_final)
}


calc_LogFoldchange <- function(comp){
  remove <- "gold"
  dat <- comp[!rownames(comp) %in% remove, ]
  
  dat <- dat[order(dat[,1], decreasing = T),]
  dat_final <- foldchange(abs(dat[,3:dim(comp)[2]]), dat[,2])
  log_fold <-foldchange2logratio(dat_final, base =2)
  fc_comparison <<- cbind(dat[,1], log_fold)
}

# Visualize Difference
#--------------------------------------------------------------------------------------------------------
# Vis_Foldchange <- function(data, method, filename){
#   svg(filename)
#   par(mfrow=c(2,3))
#   for (i in 2:dim(data)[2]){
#     if(min(data[,i]) < -10){
#       plot(y = data[,i], 
#            x = data[,1],
#            ylim = c(-180,2), 
#            col = c(rownames(data)),
#            pch = rep(19,dim(data)[1]),
#            cex = 2,
#            xlab = "Module Size", ylab = paste("Foldchange Zsummary (", method, "- Gold)"),
#            main = colnames(data)[i])
#       abline(h = 0, col = "grey", lty =2)}
#     else{
#       plot(y = data[,i], 
#            x = data[,1],
#            ylim = c(-10,2), 
#            col = c(rownames(data)),
#            pch = rep(19,dim(data)[1]),
#            cex = 2,
#            xlab = "Module Size", ylab = paste("Foldchange Zsummary (", method, "- Gold)"),
#            main = colnames(data)[i])
#       abline(h = 0, col = "grey", lty =2)}
#   }
#   dev.off()
# }

