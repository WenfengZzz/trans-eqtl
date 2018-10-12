# used to find out gene_exp and gene_snp (intra and inter)
args <- commandArgs()
library(biomaRt)
load("/Users/wenfeng/Box Sync/Graduate/Research/Hi-C/ENCFF720ZUE.h5_interaction.Rdata")

# set NAs to 0
interactions[is.na(interactions)] <- 0
# remove lower triangular matrix
interactions[lower.tri(interactions)] <- 0
# set diagonal to 0
diag(interactions) <- 0


# interactions[1:10,1:10]
# chr_bin_range
# chrs
# bin_positions[1:3,1:10]

# find location of significant interactions top 0.1%
t <- quantile(interactions, probs = 0.99)
sig_interactions <- which(interactions>t, arr.ind = TRUE)
n <- dim(sig_interactions)[1]
bin_positions_sig <- matrix(nrow = 6, ncol = n)


for (i in 1:n){
rw <- sig_interactions[i,1]
cl <- sig_interactions[i,2]
bin_positions_sig[1,i] <- bin_positions[1,rw]
bin_positions_sig[2,i] <- bin_positions[2,rw]
bin_positions_sig[3,i] <- bin_positions[3,rw]

bin_positions_sig[4,i] <- bin_positions[1,cl]
bin_positions_sig[5,i] <- bin_positions[2,cl]
bin_positions_sig[6,i] <- bin_positions[3,cl]
}


# distance between trans-regions must exceed 5Mb
trans_bin_positions_sig <- bin_positions_sig[,abs(bin_positions_sig[2,]-bin_positions_sig[5,])>5000000|
                                               bin_positions_sig[1,]!=bin_positions_sig[4,]]
# retain only autosomes
trans_bin_positions_sig <- trans_bin_positions_sig[, !(trans_bin_positions_sig[4,] %in% c(22,23,24))]

# find out regions
search_region <- function(chr, begin, end, intra){
  target <- trans_bin_positions_sig[,trans_bin_positions_sig[1,] == chr-1
                                    & trans_bin_positions_sig[2,] == begin
                                    & trans_bin_positions_sig[3,] == end]
  if (intra == TRUE){
  target <- target[4:6,target[4,] == chr-1]
  }
  else {
  target <- target[4:6,target[4,] != chr-1]
  }
}

strt <- args[3]
end <- args[4]
chr <- args[2]
chr1_intra_region <- search_region(1, strt, end, TRUE)
chr1_inter_region <- search_region(1, strt, end, FALSE)

save(chr1_inter_region, chr1_intra_region, file = "sample_region.RData")

n_intra <- dim(chr1_intra_region)[2]
n_inter <- dim(chr1_inter_region)[2]

gene_exp <- list()
gene_snp_intra <- list()
gene_snp_inter <- list()


ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
# expression genes
gene_exp[[1]] <- getBM(c("chromosome_name","ensembl_gene_id","hgnc_symbol","start_position","end_position"),
                    filters = c("chromosome_name", "start","end"),
                    values = list(chr, strt, end), mart = ensembl)

# intra genes
for (i in 1:n_intra){
  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

  chr2 <- chr1_intra_region[1,i]
  start2 <- chr1_intra_region[2,i]
  end2 <- chr1_intra_region[3,i]
  
  gene_snp_intra[[i]] <- getBM(c("chromosome_name","ensembl_gene_id","hgnc_symbol","start_position","end_position"),
                      filters = c("chromosome_name", "start","end"),
                      values = list(chr2+1, start2, end2), mart = ensembl)
}
# inter genes

for (i in 1:n_inter){
  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  
  chr2 <- chr1_inter_region[1,i]
  start2 <- chr1_inter_region[2,i]
  end2 <- chr1_inter_region[3,i]
  
  gene_snp_inter[[i]] <- getBM(c("chromosome_name","ensembl_gene_id","hgnc_symbol","start_position","end_position"),
                               filters = c("chromosome_name", "start","end"),
                               values = list(chr2+1, start2, end2), mart = ensembl)
}


save(gene_exp, gene_snp_inter, gene_snp_intra, file = "sample_gene.RData")

# # merge regions
# merged_region <- matrix(nrow = 3, ncol = 20)
# 
# merged_region[1,] <- chr1_intra_region[1,1]
# j <- 1
# merged_region[2,j] <- chr1_intra_region[2,1]
# while (chr1_intra_region[3,j] - chr1_intra_region[2,j+1] == -1) {
#   end <- chr1_intra_region[3,j+1]
#   j <- j + 1
# }  
# end

# dim(chr1_inter_region)
# chr1_inter_region[1:3, 1:6]
# chr1_intra_region[1:3, 1:6]
# 
# # find corresponding genes
# gene1 <- list()
# gene2 <- list()
# n_sig <- dim(trans_bin_positions_sig)[2]
# chr1_sig <- trans_bin_positions_sig[,trans_bin_positions_sig[1,]==0 & 
#                                       !(trans_bin_positions_sig[4,] %in% c(22,23,24))]
# n_chr1 <- dim(chr1_sig)[2]
# 
# levels(as.factor(chr1_sig[4,]))
# 
# for (i in 1:300){
#   ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
#   chr1 <- chr1_sig[1,i]
#   start1 <- chr1_sig[2,i]
#   end1 <- chr1_sig[3,i]
#   
#   chr2 <- chr1_sig[4,i]
#   start2 <- chr1_sig[5,i]
#   end2 <- chr1_sig[6,i]
#   
#   gene1[[i]] <- getBM(c("chromosome_name","ensembl_gene_id","hgnc_symbol","start_position","end_position"),
#                       filters = c("chromosome_name", "start","end"),
#                       values = list(chr1+1, start1, end1), mart = ensembl)
#   gene2[[i]] <- getBM(c("chromosome_name","ensembl_gene_id","hgnc_symbol","start_position","end_position"),
#                       filters = c("chromosome_name", "start","end"),
#                       values = list(chr2+1, start2, end2), mart = ensembl)
# }




