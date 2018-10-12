rm(list = ls())
load("sample_gene.RData")

for (i in 1:length(gene_snp_inter)){
  gene_snp_inter[[i]] <- gene_snp_inter[[i]][gene_snp_inter[[i]]$hgnc_symbol != "",]
  #gene_snp_inter[[i]] <- gene_snp_inter[[i]][1,]
}

for (i in 1:length(gene_snp_intra)){
  gene_snp_intra[[i]] <- gene_snp_intra[[i]][gene_snp_intra[[i]]$hgnc_symbol != "",]
  #gene_snp_intra[[i]] <- gene_snp_intra[[i]][1,]
}

gene_snp_inter[[1]]


paste0("chr", gene_snp_inter[[1]]$chromosome_name)

gene2 <- append(gene_snp_inter, gene_snp_intra)

file = paste0("paste_files.sh")
sink(file)

for (i in 1:length(gene2)){
  snp_gene <- as.data.frame(gene2[i])[,2]
  chr_snp <- paste0("chr", gene2[[i]]$chromosome_name)
  # annotation
  sent <- paste0("mkdir -p /home/wz284/scratch60/data_to_grace/GTEx/annotation/", chr_snp, "; ",
    "cp -r /ysm-gpfs/pi/zhao/wl382/snpPred_epi/Annotation/data/annotation/",chr_snp,"/",snp_gene,"*/",
                 " /home/wz284/scratch60/data_to_grace/GTEx/annotation/", chr_snp,"/;")
  cat(sent, "\n")
  # snp info
  sent <- paste0("mkdir -p /home/wz284/scratch60/data_to_grace/GTEx/cis_snp_by_gene/", chr_snp, "; ",
    "cp -r /ysm-gpfs/pi/zhao/from_louise/jw2372/GTEx/cis_snp_by_gene/",chr_snp,"/",snp_gene,"*/", 
                 " /home/wz284/scratch60/data_to_grace/GTEx/cis_snp_by_gene/", chr_snp,"/;")
  cat(sent, "\n")
}

# annotation
paste0("cp /ysm-gpfs/pi/zhao/wl382/snpPred_epi/Annotation/data/annotation/",chr_snp,"/",gene_snp,"*/"," /home/fas/zhao/wz284/scratch60/data/GTEx/annotation/")

# snp info
paste0("/ysm-gpfs/pi/zhao/from_louise/jw2372/GTEx/cis_snp_by_gene/",chr_snp,"/",gene_snp,"*/", " /home/fas/zhao/wz284/scratch60/data/GTEx/cis_snp_by_gene/")

# expr
paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/adjusted_expr1/",chr_exp,"/",gene_exp,"*/"," /home/fas/zhao/wz284/scratch60/data/GTEx/adjusted_expr1/")
