source("https://bioconductor.org/biocLite.R")
install.packages('yaml')
library(IRanges)
library(biomaRt)
library(rhdf5)

h5_file <- paste0("/Users/wayne/Box Sync/Graduate/Research/Hi-C/ENCFF007MJG.h5") 

# import h5 file

h5_group <- h5ls(h5_file)[,2]

if(length(h5_group)==4){
  bin_positions <- h5read(h5_file, h5_group[1])
  chr_bin_range <- h5read(h5_file, h5_group[2])
  chrs <- h5read(h5_file, h5_group[3])
  interactions <- h5read(h5_file, h5_group[4])
} else 
{
  balance_factors <- h5read(h5_file, h5_group[1])
  bin_positions <- h5read(h5_file, h5_group[2])
  chr_bin_range <- h5read(h5_file, h5_group[3])
  chrs <- h5read(h5_file, h5_group[4])
  interactions <- h5read(h5_file, h5_group[5])
}

#interactions[1:10, 1:10]

#bin_positions[1:3, 495:500]

interactions[is.na(interactions)] <- 0

# balancing
#while(TRUE){
#  test <- interactions[100,100]
#  interactions <- interactions/(rowMeans(interactions))
#  interactions <- interactions/(colMeans(interactions))
#  if (abs(interaction[100,100]-t) < 0.00001) break
#}

# remove lower triangular matrix
interactions[lower.tri(interactions)] <- 0
interactions[1:10, 1:10]

# find location of significant interactions top 1%
t <- quantile(interactions, probs = 0.99, na.rm = TRUE)
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


trans_bin_positions_sig <- bin_positions_sig[,abs(bin_positions_sig[2,]-bin_positions_sig[5,])>5000000|
                                               bin_positions_sig[1,]!=bin_positions_sig[4,]]

dim(interactions)

sum(trans_bin_positions_sig[1,]==trans_bin_positions_sig[4,])

dim(trans_bin_positions_sig)
1280/1310

heatmap(interactions[1:30,1:30],Rowv=NA, Colv=NA)
?heatmap

# find corresponding genes
gene1 <- list()
gene2 <- list()


for (i in 1:10){ #dim(bin_positions)[2]){
  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  chr1 <- trans_bin_positions_sig[1,i]
  start1 <- trans_bin_positions_sig[2,i]
  end1 <- trans_bin_positions_sig[3,i]

  chr2 <- trans_bin_positions_sig[4,i]
  start2 <- trans_bin_positions_sig[5,i]
  end2 <- trans_bin_positions_sig[6,i]

  gene1[[i]] <- getBM(c("chromosome_name","ensembl_gene_id","hgnc_symbol","start_position","end_position"), 
      filters = c("chromosome_name", "start","end"), 
      values = list(chr1+1, start1, end1), mart = ensembl)
  gene2[[i]] <- getBM(c("chromosome_name","ensembl_gene_id","hgnc_symbol","start_position","end_position"), 
      filters = c("chromosome_name", "start","end"), 
      values = list(chr2+1, start2, end2), mart = ensembl)
}


chr = "chr1"

gene_exp <- gene1[[1]][,2]
gene_snp <- gene2[[2]][,2]


expr_dir = paste0("D:/Yale/",chr,"/",gene_exp,"/")
gene_dir = paste0("D:/Yale/snp/",chr,"/",gene_snp,"/")


# read in expression data

tissue = "Whole_Blood"

expr.file <- paste0(expr_dir,tissue,".adj_expr")
dosage_file <- paste0(gene_dir,gene_snp,".mach.dose.gz")
result_dir = paste0("D:/Yale/")


for (i in 1:length(gene_exp)){
  for (j in 1:length(gene_snp)){
  i = 1
  j = 1
  expr.file <- expr.file[i]
  dosage_file <-  dosage_file[j]
  # dosage_file <- "D:/Yale/snp/chr1/ENSG00000230021/ENSG00000230021.mach.dose.gz"
  expr <- read.table(expr.file, header= F, stringsAsFactors=F, fill=T)
  sample.names <- as.vector(expr[,1])
  M = 1
  N = length(sample.names)
  y.array <- array(data=expr[,2], dim = c(N,1,M),
                     dimnames = list(c(sample.names),c("expr"),
                                     c(tissue)))
    
#################load genotype data for trans-SNPs#####################
    
  gen_dat = read.table(gzfile(dosage_file),header=F,stringsAsFactors=F,fill=T)
  gen_dat.1 = gen_dat[complete.cases(gen_dat),]
  sample.names.gen <- gen_dat.1[,1]
  sample.names.gen.1 <- c()
  sample.strsp <- strsplit(sample.names.gen,"->")
    
    
  for(i in 1:length(sample.strsp)){
      name.tmp <- sample.strsp[[i]][1]
      sample.names.gen.1 <- c(sample.names.gen.1,name.tmp)
  }
    
  p <- dim(gen_dat.1)[2]
    
  ## find the sample intersection between genotype data and expr data
  sample.names.inter <- intersect(sample.names.gen.1,sample.names)
  
  #gen.sample = dim(gen_dat)[2]-6
  x.array<-array(0,dim=c(length(sample.names.inter),p,M))
  
  gen_sort = gen_dat.1[match(sample.names.inter, sample.names.gen.1),-c(1,2)]
  # only extract snps which also have annotation information
  
  gen_sort_1 <- gen_sort
  
  for(i in 1:dim(gen_sort_1)[2]){
    # x.array[,,i] = scale(gen_sort_1, scale = FALSE)
    # x.array[,i,1] = gen_sort_1[,i]
    x.array[,i,1] = scale(gen_sort_1[,i],center=T, scale = FALSE)
    if(is.na(x.array[1,i,1])){
      x.array[,i,1] = gen_sort_1[,i]
    }
  }
  
  ### extract y.array by sample intersections
  
  
  #expr.sample <- scale(expr[match(sample.names.inter,sample.names),2],center=T, scale=T)
  expr.sample <- as.array(expr[match(sample.names.inter,sample.names),2])
  exprs.sample <- c()
  for(i in 1:(dim(expr.sample)[1])){
    exprs.sample <- c(exprs.sample,expr.sample[i])
  }
  

  # sample size for gene expression
  n.expr = length(sample.names.inter)
  
  # let y be the expression array
  y.array <- array(data=exprs.sample, dim = c(length(sample.names.inter),1,M),
                   dimnames = list(c(sample.names.inter),c("expr"),
                                   c(tissue)
                   )
  )
  
  
  # elnt model
  
  
  y = y.array[,1,1]
  x = x.array[,,1]
  x[1:10, 1:10]
  dim(x)
  length(y)
  coeff <- numeric()
  p.val <- numeric()
  r2 <- numeric()
  for (k in 1:1300){
    m.lm <- summary(lm(y~x[,k]))
    coeff[k] <- m.lm$coefficients[2,1]
    p.val[k] <- m.lm$coefficients[2,4]
    r2[k] <- m.lm$r.squared
  }
  p.val[which(p.val < 0.01)]
  r2[which(p.val < 0.01)]
  min(p.val)
  max(r2)
  which.min(p.val)
  which.max(r2)
  
  #1_5413282_G_T_b37	G	T	-		-	-	Typed_Only	-	-	-	-	-
  #1_5413994_T_C_b37	T	C	-		-	-	Typed_Only	-	-	-	-	-
  #1_5415884_C_T_b37	C	T	-		-	-	Typed_Only	-	-	-	-	-
  
  
  cv.elnt = cv.glmnet(x,y, alpha=0.5)
  rsq <- 1-cv.elnt$cvm/var(y)
  betas.elnt = coef(cv.elnt,s="lambda.1se")
  
  sum(betas.elnt!=0)
  save(betas.elnt, file=paste0(result_dir,"trans",i,j,".Rdata"))
  } 
}


