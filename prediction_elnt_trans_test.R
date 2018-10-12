library(foreach, lib = "/gpfs/loomis/scratch60/fas/zhao/wz284/R")
library(glmnet, lib = "/gpfs/loomis/scratch60/fas/zhao/wz284/R")

args <- commandArgs()
# args[6]: chr; args[7]: exp_gene; args[8]: snp_gene

chr = args[6]
gene_snp = args[7] 
gene_exp = args[8]
tissue = "Liver"

gene_exp = "ENSG00000157873"
gene_snp = "ENSG00000226487"
chr = "chr1"

# get tissues fo the gene
annot.dir <- paste0("/home/fas/zhao/wz284/scratch60/data/GTEx/annotation/",chr,"/",gene_snp,"/")
# files = list.files(annot.dir)

gene_dir = paste0("/home/fas/zhao/wz284/scratch60/data/GTEx/cis_snp_by_gene/",chr,"/",gene_snp,"/")
expr_dir = paste0("/home/fas/zhao/wz284/scratch60/data/GTEx/adjusted_expr1/",chr,"/",gene_exp,"/")


if(!file.exists(gene_dir)){
  next()
}
snp_files = list.files(gene_dir)
snp_info_file = paste0(gene_dir, grep("mach.info", snp_files, value = TRUE))
if(!file.exists(snp_info_file)){
  next()
}
snp_info = read.table(snp_info_file, header = T, stringsAsFactors = F,fill=T)
result_dir = paste0("/home/fas/zhao/wz284/scratch60/result/trans/",chr,"/",gene_exp,"/")

  #################data loading#####################


  ################load annotation information########################

  if(!file.exists(annot_file)){
    annot_file = paste0(annot.dir,"annotation_",tissue,".txt")
    if(!file.exists(annot_file)){
    next()
   }
  }

  annot = data.frame(read.table(annot_file,header=T,stringsAsFactors=F))


  annot.list <- c()
  for(j in 4:dim(annot)[2]){
    for (i in 1:dim(annot)[1]){
      annot.list <- c(annot.list,annot[i,j])
    }
  }

  
  ################## for one tissue
  annotat.array <- array(annot.list,dim=c(dim(annot)[1],dim(annot)[2]-3))
  # delete na's
  na.rows <- which(is.na(annotat.array), arr.ind=TRUE)[,1]
  if(length(na.rows) != 0){
    annotat.array.noNa <- array(annotat.array[-as.numeric(na.rows),], dim=c(dim(annot)[1]-length(na.rows),dim(annot)[2]-3)) # for a single tissue
    annotat.array <- annotat.array.noNa
  }
  na.cols <- na.rows

  
  
  ################read expr data######################
  expr.file <- paste0(expr_dir,tissue,".adj_expr")
  if(!file.exists(expr.file)){
    next()
  }
  
  expr <- read.table(expr.file, header= F, stringsAsFactors=F, fill=T)
  sample.names <- as.vector(expr[,1])
  M = 1
  N = length(sample.names)
  y.array <- array(data=expr[,2], dim = c(N,1,M),
                   dimnames = list(c(sample.names),c("expr"),
                                   c(tissue)))
  

  
  #################load genotype data for cis-SNPs#####################
  # gen_file = paste0(gene_dir,"union_eqtl_allTissue_gen_AD.txt")

  dosage_file = paste0(gene_dir, grep("mach.dose", snp_files, value = TRUE))
#   if(!file.exists(dosage_file)){
#  #   next()
# 	dosage_file = paste0(gene_dir,gene_snp,".mach.dose")
#   }

  gen_dat = read.table(gzfile(dosage_file),header=F,stringsAsFactors=F,fill=T)
  gen_dat.1 = gen_dat[complete.cases(gen_dat),]
  sample.names.gen <- gen_dat.1[,1]
  sample.names.gen.1 <- c()
  sample.strsp <- strsplit(sample.names.gen,"->")
  
  
  for(i in 1:length(sample.strsp)){
    name.tmp <- sample.strsp[[i]][1]
    sample.names.gen.1 <- c(sample.names.gen.1,name.tmp)
  }
  
  p <- dim(annotat.array)[1]
  

  ## find the sample intersection between genotype data and expr data
  sample.names.inter <- intersect(sample.names.gen.1,sample.names)
  
  #gen.sample = dim(gen_dat)[2]-6
  x.array<-array(0,dim=c(length(sample.names.inter),p,M))
  
  gen_sort = gen_dat.1[match(sample.names.inter, sample.names.gen.1),-c(1,2)]
  # only extract snps which also have annotation information
  if(length(na.cols) != 0 ){
    gen_sort_1 = gen_sort[,-as.numeric(na.cols)]
  }else{
    gen_sort_1 = gen_sort
  }
  
  
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
  
  

  #################parameter set up#####################
  
  #number of tissue
  M<-dim(x.array)[3]
  
  #number of obs in each dataset
  n <- rep(dim(x.array)[1],M)
  
  #number of covariates
  p<-dim(x.array)[2]
  
  
  y = y.array[,1,1]
  x = x.array[,,1]
  annotat.array = cbind(rep(1,dim(x)[2]),annotat.array)
  ######## for annotation selected models #########
  # ind.H3K4me1 = which(annot$H3K4me1>0)
  # ind.H3K4me3 = which(annot$H3K4me3>0)
  # ind.H3K9ac = which(annot$H3K9ac>0)
  # ind.1 = union(ind.H3K4me1,ind.H3K4me3)
  # ind.2 = union(ind.1,ind.H3K9ac)
  #################################################
  # mde.varbvs = varbvs(X = x,y=y,Z=NULL,family="gaussian")
  mde.elnt = cv.glmnet(x,y, alpha=0.5)
  # mde.varbvs.annot = varbvs(X = x[,ind.2],y=y,Z=NULL,family="gaussian")
  # mde.elnt.annot = cv.glmnet(x[,ind.2],y, alpha=0.5)
  
  ######################################################
  ############### parameters retrival ##################
  ######################################################
  betas.elnt = coef(mde.elnt,s="lambda.min")
  # betas.elnt.annot = coef(mde.elnt.annot,s="lambda.min")
  # top.vars.vb = summary(mde.varbvs)$top.vars
  # top.vars.vb.annot = summary(mde.varbvs.annot)$top.vars
  # betas.vb = mde.varbvs$beta
  # betas.vb.annot = mde.varbvs.annot$beta

  save(x.array,y.array,
    betas.elnt,
    file=paste0(result_dir,gene_snp,"_",tissue,"_elnt.Rdata"))  
#}
#save(x.array,y.array,annotat.array,ind.2,
#    betas.elnt,betas.elnt.annot,
#    file=paste0(result_dir,gene,"_",tissue,"_elnt_model.Rdata"))



