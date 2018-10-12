#this script is used to find out gene_exp and gene_snp (intra and inter)

library(biomaRt)
load("/Users/wenfeng/Box Sync/Graduate/Research/Hi-C/ENCFF720ZUE.h5_interaction.Rdata")
#-----------------------------------------------------------------------------------------


# data preprocessing ---------------------------------------------------------------------

# set NAs to 0
interactions[is.na(interactions)] <- 0
# remove lower triangular matrix
interactions[lower.tri(interactions)] <- 0
# set diagonal to 0
diag(interactions) <- 0

# find location of significant interactions top 0.1%
t <- quantile(interactions, probs = 0.99)
sig_interactions <- which(interactions>t, arr.ind = TRUE)
n <- dim(sig_interactions)[1]
bin_positions_sig <- matrix(nrow = 6, ncol = n)

# create dataframe to store regions
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

trans_bin_positions_sig[1,] <- trans_bin_positions_sig[1,] + 1
trans_bin_positions_sig[4,] <- trans_bin_positions_sig[4,] + 1
#-----------------------------------------------------------------------------------------------------

#function used to find out snp_regions for a specific exp_gene region
search_region <- function(chr, begin, end){
  target <- trans_bin_positions_sig[4:6,trans_bin_positions_sig[1,] == chr
                                    & trans_bin_positions_sig[2,] == begin
                                    & trans_bin_positions_sig[3,] == end]
  return(target)
}

# bins of expression genes --------------------------------------------------------------------------
exp_gene_bin <- unique(t(trans_bin_positions_sig[1:3, ]))
exp_gene_bin[,1] <- exp_gene_bin[,1]
head(exp_gene_bin)
# expression genes in exp_gene bins
dim(exp_gene_bin)
for (j in 5173:5180){
  chr <- exp_gene_bin[j,1]
  strt <- exp_gene_bin[j,2]
  end <- exp_gene_bin[j,3]
  
  gene_exp <- list()
  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  gene_exp[[1]] <- getBM(c("chromosome_name","ensembl_gene_id","hgnc_symbol","start_position","end_position"),
                         filters = c("chromosome_name", "start","end"),
                         values = list(chr, strt, end), mart = ensembl)
  
  # bins of corresponding snps
  snp_region <- t(search_region(chr, strt, end))
  
  exp_genes <- gene_exp[[1]][,2]
  exp_genes[1]
  for (i in 1:length(exp_genes)){
    result <- paste0("/Users/wenfeng/Box Sync/Graduate/Research/Hi-C/exp_gene_snp_region/chr",chr,"/",exp_genes[i], "_snpregion.Rdata")
    exp_gene=exp_genes[i]
    save(exp_gene, snp_region, file = result)
  }
}

