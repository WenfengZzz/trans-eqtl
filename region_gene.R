rm(list=ls())
load("sample_gene.RData")
gene_snp_inter[[6]]
gene_snp_intra[[5]]
gene_intra <- c()
for (i in 1:87) {
  gene_snp_intra[[i]] <- gene_snp_intra[[i]][gene_snp_intra[[i]]$hgnc_symbol != "",]
  gene_intra[i] <- gene_snp_intra[[i]][1,2]
}
intra <- rbind(chr1_intra_region,  gene_intra)
rownames(intra) <- NULL

load("sample_region.RData")
chr1_inter_region[1,] <- chr1_inter_region[1,] + 1

chr1_intra_region
chr1_inter_region

# chr1
# ENSG00000120948_Liver_elnt.Rdata
# ENSG00000116688_Liver_elnt.Rdata
# ENSG00000162494_liver_elnt.Rdata
# 
# chr3 ENSG00000244674_Liver_elnt.Rdata
# chr9 ENSG00000196873_Liver_elnt.Rdata
# chr16 ENSG00000260518_Liver_elnt.Rdata
# chr20 NA
# chr21 NA


gene_intra

file = "snp_info_paste.sh"
sink(file)
i <- 1
for (i in 1:length(gene_intra)){
  # annotation
  sent <- paste0("cp -r /home/fas/zhao/wz284/scratch60/data/GTEx/cis_snp_by_gene/chr1/",gene_intra[i],"/",
                 " /home/fas/zhao/wz284/scratch60/move_to_mac/;")
  cat(sent, "\n")
}
