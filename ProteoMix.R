pcname <- system('uname -n',intern=T)
if (pcname=='crosswood.stanford.edu') {
  Root = "~/Proteomics"
}else{
  Root = "/Users/mmagzoub/Desktop/Proteomics"
}
study.dir = paste0(Root, "/Study")
correlation.dir = paste0(Root, "/PostProcessing/DataComparison")

strip_clusters <- function(genes){
  clusters <- genes[grepl("---", genes)]
  genes[grepl("---", genes)] <- unlist(lapply(clusters, function(x) unlist(strsplit(x, "---"))[[1]]))
  return(genes)
}

#cancer = "BRCA"
#cancer = "COADREAD"
cancer = "OV"

####################
#load all data
####################

#MET
load(paste0(study.dir, "/data/", cancer, "/methylmix_data.Rdata"))
MET <- METcancer
MET.normal <- METnormal
all.cpg_clusters <- intersect(rownames(MET), rownames(MET.normal))

#MRNA
MRNA <- MAcancer
MRNA.normal <- MAnormal

#PROT
load(paste0(study.dir, "/data/", cancer, "/protein_abundance_data.Rdata"))

rm(MAcancer)
rm(MAnormal)
rm(METcancer)
rm(METnormal)
rm(ProbeMapping)

# study data
samples <- Reduce(intersect, list(colnames(MET), colnames(MRNA), colnames(PROT)))
genes <- Reduce(intersect, list(rownames(MRNA), rownames(PROT),
                                strip_clusters(all.cpg_clusters)))
cpg.clusters <- all.cpg_clusters[strip_clusters(all.cpg_clusters) %in% genes]
# match data
PROT <- PROT[genes, samples]
MRNA <- MRNA[genes, samples]
MET <- MET[cpg.clusters, ]
MET.normal <- MET.normal[cpg.clusters,]

study.data <- list('MRNA' = MRNA, 'MET' = MET, 'MET.normal' = MET.normal, 'PROT' = PROT)
save(study.data, file = paste0(study.dir, '/data/', cancer, '/study_data.Rdata'))

####################
#run models
####################

library(MethylMixBeta)
source(paste0(study.dir, "/Update_MethylMix.R"))
MethylMix_Update()

proteomix <- MethylMix(MET, PROT, MET.normal)
save(proteomix, file = paste0(study.dir, "/results/", cancer, "/proteomix_output.Rdata" ))
methylmix <- MethylMix(MET, MRNA, MET.normal)
save(methylmix, file = paste0(study.dir, "/results/", cancer, "/methylmix_output.Rdata" ))

proteomix.cpg_clusters <- proteomix$MethylationDrivers
methylmix.cpg_clusters <- methylmix$MethylationDrivers

proteomix.missing <- MethylMix(MET, PROT, MET.normal, setdiff(methylmix.cpg_clusters, proteomix.cpg_clusters), F)
save(proteomix.missing, file = paste0(study.dir, "/results/", cancer, "/proteomix_missing.Rdata"))
methylmix.missing <- MethylMix(MET, MRNA, MET.normal, setdiff(proteomix.cpg_clusters, methylmix.cpg_clusters), F)
save(methylmix.missing, file = paste0(study.dir, "/results/", cancer, "/methylmix_missing.Rdata"))