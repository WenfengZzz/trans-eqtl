# form snp bins to snp location
args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
snp_bin_file <- args[2]
load("/home/wz284/scratch60/trans_eqtl/data/trans_snp/SNP.RData")

bin_to_location <- function(snp_bin_file){
  load(snp_bin_file)
  snps <- list()
  for (i in 1:dim(snp_region)[1]) {
    snps[[i]] <- as.matrix(SNPinfo[SNPinfo[,1] == snp_region[i,1] & 
                                     SNPinfo[,2]>snp_region[i,2] & SNPinfo[,2]<snp_region[i,3],][,1:2])
    rownames(snps[[i]]) <- c()
    colnames(snps[[i]]) <- c()  
  }
  snps <- as.matrix(do.call(rbind, snps))
  write.table(snps, file=paste0("/home/wz284/scratch60/trans_eqtl/data/trans_snp/snp_position/", chr, "/",  exp_gene, "_snp_position"), row.names=FALSE, col.names=FALSE)  
}

bin.dir <- paste0("/home/wz284/scratch60/trans_eqtl/data/trans_snp/exp_gene_snp_region/",chr,"/")

# snp_bin_files <- list.files(bin.dir)
# file <- snp_bin_files[id]
file <- paste0(bin.dir, snp_bin_file)
bin_to_location(file)

