# rm(list=ls())
load("sample_gene.RData")

for (i in 1:length(gene_snp_inter)){
  gene_snp_inter[[i]] <- gene_snp_inter[[i]][gene_snp_inter[[i]]$hgnc_symbol != "",]
  #gene_snp_inter[[i]] <- gene_snp_inter[[i]][1,]
}

for (i in 1:length(gene_snp_intra)){
  gene_snp_intra[[i]] <- gene_snp_intra[[i]][gene_snp_intra[[i]]$hgnc_symbol != "",]
  #gene_snp_intra[[i]] <- gene_snp_intra[[i]][1,]
}

paste0("chr", gene_snp_inter[[1]]$chromosome_name)

i <- 1
gene1 <- gene_exp
gene2 <- gene_snp_inter

sh.generate <- function(gene1, gene2){
  for (i in 1:length(gene2)){
    exp_gene <- as.data.frame(gene1[1])[,2]
    snp_gene <- as.data.frame(gene2[i])[,2]
    chr_exp <- paste0("chr", gene1[[1]]$chromosome_name[1])
    chr_snp <- paste0("chr", gene2[[i]]$chromosome_name)
    for (j in 1:length(exp_gene)){
      for (k in 1:length(snp_gene)){
        result.dir = paste0("/home/fas/zhao/wz284/scratch60/result/trans/",chr_exp,"/",exp_gene[j],"/")
        sent <- paste0("mkdir -p ",result.dir,
                       "; module load Apps/R/3.4.3-generic; Rscript /home/fas/zhao/wz284/project/hic/script/transeqtl_elnt.R ",
                       chr_exp," ",exp_gene[j], " ",chr_snp[k], " ", snp_gene[k])
        cat(sent,"\n")
      }
    }
  }
}

file = paste0("sample_inter.sh")
sink(file)
sh.generate(gene_exp, gene_snp_inter)

file = paste0("sample_intra.sh")
sink(file)
sh.generate(gene_exp, gene_snp_intra)



# load("chr1_genes.Rdata")
# gene2[190]
# file = paste0("chr1_391.sh")
# sink(file)
# for (i in 1:length(gene1)){
#   exp_gene <- as.data.frame(gene1[i])[,2]
#   snp_gene <- as.data.frame(gene2[i])[,2]
#   for (j in 1:length(exp_gene)){
#     for (k in 1:length(snp_gene)){
#       chr = "chr1"
#       result.dir = paste0("/home/fas/zhao/wz284/scratch60/result/trans/",chr,"/",exp_gene[j],"/")
#       sent <- paste0("mkdir -p ",result.dir,"; module load Apps/R/3.4.3-generic; Rscript /home/fas/zhao/wz284/project/hic/script/transeqtl_elnt.R ",chr," ",exp_gene[j], " ",chr, " ", snp_gene[k])
#       cat(sent,"\n")
#     }
#   }
# }
# 











# rm(list=ls())
load("sample_gene.RData")

for (i in 1:length(gene_snp_inter)){
  gene_snp_inter[[i]] <- gene_snp_inter[[i]][gene_snp_inter[[i]]$hgnc_symbol != "",]
  gene_snp_inter[[i]] <- gene_snp_inter[[i]][1,]
}

for (i in 1:length(gene_snp_intra)){
  gene_snp_intra[[i]] <- gene_snp_intra[[i]][gene_snp_intra[[i]]$hgnc_symbol != "",]
  gene_snp_intra[[i]] <- gene_snp_intra[[i]][1,]
}

paste0("chr", gene_snp_inter[[1]]$chromosome_name)

i <- 1
sh.generate <- function(gene1, gene2){
  for (i in 1:length(gene2)){
    exp_gene <- as.data.frame(gene1[1])[,2]
    snp_gene <- as.data.frame(gene2[i])[,2]
    chr_exp <- paste0("chr", gene1[[1]]$chromosome_name[1])
    chr_snp <- paste0("chr", gene2[[i]]$chromosome_name)
    for (j in 1:length(exp_gene)){
      for (k in 1:length(snp_gene)){
        result.dir = paste0("/home/fas/zhao/wz284/scratch60/result/trans/",chr_exp,"/",exp_gene[j],"/")
        sent <- paste0("mkdir -p ",result.dir,"; module load Apps/R/3.4.3-generic; Rscript /home/fas/zhao/wz284/project/hic/script/transeqtl_elnt.R ",chr_exp," ",exp_gene[j], " ",chr_snp, " ", snp_gene[k])
        cat(sent,"\n")
      }
    }
  }
}

file = paste0("sample_inter.sh")
sink(file)
sh.generate(gene_exp, gene_snp_inter)

file = paste0("sample_intra.sh")
sink(file)
sh.generate(gene_exp, gene_snp_intra)


