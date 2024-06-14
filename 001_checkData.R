library(fs)
# this are the vst transformed counts
c_vst <- readRDS("data/c_vst.rds")
str(assay(c_vst))
m_vst <- assay(c_vst)
outDir <- fs::dir_create("014_WGCNA_networkCourseThomas_out")
saveRDS(m_vst, paste(outDir, "m_vst.rds", sep = "/"))
# this is the colData
dt_colDataNew <- readRDS("/scratch/rgraus29/projects/rna_seq_neha/analysis/002_prepareAnnoDataToCheckOnPcaHm_out/dt_colDataNew.rds")
table(dt_colDataNew$RNACluster)
