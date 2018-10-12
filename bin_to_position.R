# form snp bins to snp location

load("/Users/wenfeng/Box Sync/Graduate/Research/Hi-C/SNP.RData")

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
  write.table(snps, file=paste0("snp_locations/", exp_gene, "_snp_position"), row.names=FALSE, col.names=FALSE)  
}

snp_bin_files <- list.files("/Users/wenfeng/Box Sync/Graduate/Research/Hi-C/exp_gene_snp_region/chr1")

for (i in 1:length(snp_bin_files)) {
  file <- paste0("/Users/wenfeng/Box Sync/Graduate/Research/Hi-C/exp_gene_snp_region/chr1/", snp_bin_files[i])
  bin_to_location(file)
}
file <- paste0("/Users/wenfeng/Box Sync/Graduate/Research/Hi-C/exp_gene_snp_region/chr1/", snp_bin_files[i])
bin_to_location(file)