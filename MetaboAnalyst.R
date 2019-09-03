
metanr_packages <- function(){
  
  metr_pkgs <- c("Rserve", "ellipse", "scatterplot3d", "Cairo", "randomForest", "caTools", "e1071", "som", "impute", "pcaMethods", "RJSONIO", "ROCR", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "pheatmap", "SSPA", "sva", "Rcpp", "pROC", "data.table", "limma", "car", "fitdistrplus", "lars", "Hmisc", "magrittr", "methods", "xtable", "pls", "caret", "lattice", "igraph", "gplots", "KEGGgraph", "reshape", "RColorBrewer", "tibble", "siggenes", "plotly", "xcms", "CAMERA", "fgsea", "MSnbase", "BiocParallel", "metap", "reshape2", "scales")
  
  list_installed <- installed.packages()
  
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs)!=0){
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(new_pkgs, version = "3.8")
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

metanr_packages()

# Step 1: install devtools
#install.packages("devtools")
#library(devtools)

# Step 2: Install MetaboAnalystR with documentation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("impute", "pcaMethods", "globaltes", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes", "xcms", "CAMERA", "fgsea", "MSnbase", "BiocParallel"))

BiocManager::install("plotly")

devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))

#load package but install no updates
library(MetaboAnalystR)

# Step 3: Analyze

mSet<-InitDataObjects("conc", "stat", FALSE)
mSet<-Read.TextData(mSet, "../../Downloads/MetaboAnalyst_data.csv - Sheet1.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
mSet<-ANOVA.Anal(mSet, F, 0.05, "fisher", FALSE)
mSet<-PlotANOVA(mSet, "aov_0_", "png", 72, width=NA)
mSet<-ANOVA.Anal(mSet, F, 0.05, "fisher", FALSE)
mSet<-PlotANOVA(mSet, "aov_1_", "png", 72, width=NA)
mSet<-PlotCorrHeatMap(mSet, "corr_0_", "png", 72, width=NA, "col", "pearson", "bwm", "overview", F, F, F, 100)
mSet<-PlotCorrHeatMap(mSet, "corr_1_", "png", 72, width=NA, "row", "spearman", "bwm", "detail", F, F, F, 999)
mSet<-PCA.Anal(mSet)
mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
mSet<-PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
mSet<-PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1,2);
mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", "png", 72, width=NA, 1,2)
mSet<-PlotPCA3DScoreImg(mSet, "pca_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
mSet<-PlotPCA3DLoading(mSet, "pca_loading3d_0_", "json", 1,2,3)
mSet<-PLSR.Anal(mSet, reg=TRUE)
mSet<-PlotPLSPairSummary(mSet, "pls_pair_0_", "png", 72, width=NA, 5)
mSet<-PlotPLS2DScore(mSet, "pls_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
mSet<-PlotPLS3DScoreImg(mSet, "pls_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
mSet<-PlotPLSLoading(mSet, "pls_loading_0_", "png", 72, width=NA, 1, 2);
mSet<-PlotPLS3DLoading(mSet, "pls_loading3d_0_", "json", 1,2,3)
mSet<-PLSDA.CV(mSet, "L",5, "Q2")
mSet<-PlotPLS.Classification(mSet, "pls_cv_0_", "png", 72, width=NA)
mSet<-PlotPLS.Imp(mSet, "pls_imp_0_", "png", 72, width=NA, "vip", "Comp. 1", 15,FALSE)
mSet<-SAM.Anal(mSet, "d.stat", FALSE, TRUE)
mSet<-PlotSAM.FDR(mSet, 0.6, "sam_view_0_", "png", 72, width=NA)
mSet<-SetSAMSigMat(mSet, 0.6)
mSet<-PlotSAM.Cmpd(mSet, "sam_imp_0_", "png", 72, width=NA)
mSet<-PlotHCTree(mSet, "tree_0_", "png", 72, width=NA, "euclidean", "ward.D")
mSet<-PlotHCTree(mSet, "tree_1_", "png", 72, width=NA, "euclidean", "ward.D")
mSet<-PlotHCTree(mSet, "tree_2_", "png", 72, width=NA, "spearman", "ward.D")
mSet<-PlotHCTree(mSet, "tree_3_", "png", 72, width=NA, "pearson", "ward.D")
mSet<-PlotHCTree(mSet, "tree_4_", "png", 72, width=NA, "spearman", "average")
mSet<-PlotHeatMap(mSet, "heatmap_0_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "overview", T, T, NA, T, F)
mSet<-PlotSubHeatMap(mSet, "heatmap_1_", "png", 72, width=NA, "norm", "row", "euclidean", "average","bwm", "tanova", 25, "overview", T, T, T, F)
mSet<-SOM.Anal(mSet, 1,3,"linear","gaussian")
mSet<-PlotSOM(mSet, "som_0_", "png", 72, width=NA)
mSet<-Kmeans.Anal(mSet, 3)
mSet<-PlotKmeans(mSet, "km_0_", "png", 72, width=NA)
mSet<-RF.Anal(mSet, 500,7,1)
mSet<-PlotRF.Classify(mSet, "rf_cls_0_", "png", 72, width=NA)
mSet<-PlotRF.VIP(mSet, "rf_imp_0_", "png", 72, width=NA)
mSet<-PlotRF.Outlier(mSet, "rf_outlier_0_", "png", 72, width=NA)
mSet<-RF.Anal(mSet, 500,11,1)
mSet<-PlotRF.Classify(mSet, "rf_cls_1_", "png", 72, width=NA)
mSet<-PlotRF.VIP(mSet, "rf_imp_1_", "png", 72, width=NA)
mSet<-PlotRF.Outlier(mSet, "rf_outlier_1_", "png", 72, width=NA)
mSet<-SaveTransformedData(mSet)


PreparePDFReport(mSetObj = mSet, usrName = "Derreck Carter-House")
