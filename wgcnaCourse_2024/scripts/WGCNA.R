source("InitializeProject.R")
library(tidyverse)
# Load the package
library(WGCNA)
library(Biobase)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Read in the female liver data set


eset <- rio::import(here::here(datadir, "GSE21034.ExpressionSet.RDS"))
datExpr <- t(exprs(eset))

# outlier detection
sampleTrees = hclust(dist(datExpr), method = "average")

# define the params for WGCNA
params <- list(networkType = "unsigned",
               corType = "pearson",
               maxBlockSize = 24000,
               # TOMType = "signed",
               minModuleSize = 30,
               reassignThreshold = 1e-6,
               detectCutHeight = .998,
               mergeCutHeight = 0.15,
               deepSplit = 4,
               numericLabels = TRUE,
               pamStage = FALSE,
               pamRespectsDendro = TRUE,
               verbose = 6,
               datExpr = datExpr)


# enable parallel computing
enableWGCNAThreads()

# Choose a set of soft-thresholding powers
powers <- c(seq(1,10,by=1), seq(12,20, by=2))

powerTables <- pickSoftThreshold(params$datExpr,
                                 RsquaredCut = 0.75,
                                 powerVector=powers,
                                 networkType = params$networkType,
                                 verbose = 6)
collectGarbage()

params$beta <- powerTables$powerEstimate

# original call
#net <- blockwiseModules(datExpr, 
#                        power = beta,
#                        maxBlockSize = 25000,
#                        minModuleSize = 30, 
#                        deepSplit = 3,
#                        pamRespectsDendro = FALSE,
#                        mergeCutHeight = 0.15, 
#                        numericLabels = TRUE,
#                        minKMEtoStay = 0,
#                        saveTOMs = FALSE,
#                        loadTOMS = TRUE,
#                        pamStage = FALSE,
#                        verbose = 5)
# modofied call

net <- do.call(blockwiseModules, c(params))

# attach the underlying parameters and data
net$params <- as.list(args(blockwiseModules))
net$params[names(params)] <- params

# attach the underlying eset
net$eset <- eset

# save the net object and we are finished
net %>%
  rio::export(here::here(datadir, "WGCNA_net.RDS"))

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
dendro = net$dendrograms[[1]]
#pdf(file = "Plots/ConsensusDendrogram-auto.pdf", wi = 8, he = 6)
plotDendroAndColors(dendro, 
                    moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
pdata <- pData(eset)

Traits <- pdata %>%
  select(SurvTime, BCR_FreeTime, Type) %>%
  mutate(Type = case_when(Type == "PRIMARY" ~ 0,
                          TRUE ~ 1))

moduleTraitCor = cor(net$MEs, Traits, use = "p")
moduleTraitPvalue = corPvalueFisher(moduleTraitCor, nrow(pdata))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor),
               yLabels = rownames(moduleTraitCor),
               #               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships across\n"))


