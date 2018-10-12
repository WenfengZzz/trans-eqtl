args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
gene <- args[2]
tissue <- "Liver"

print(chr)
print(gene)

# chr <- "chr1"
# gene <- "ENSG00000007341.14"
# tissue <- "Liver"


#-------------------------- expression data --------------------------------
expr_dir <- paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/adjusted_expr1/",
	chr, "/", gene, "/")
expr.file <- paste0(expr_dir,tissue,".adj_expr")
if(!file.exists(expr.file)){
	print("Error: expr.file")
	next()
}
expr <- read.table(expr.file, header= F, stringsAsFactors=F, fill=T)
expr <- expr[complete.cases(expr),]
sample.names <- as.vector(expr[,1])
M = 1
N = length(sample.names)


#-------------------------- snp data --------------------------------
snp_gene <- substr(gene, start = 1, stop = 15)
gene_dir <- paste0("/ysm-gpfs/pi/zhao/wz284/trans_eqtl/data/trans_snp/snp_info/",
	chr, "/", snp_gene, "/")
if(!file.exists(gene_dir)){
  print("Error: gene_dir")
  next()
}
snp_files <- list.files(gene_dir)

dosage_file = paste0(gene_dir, grep("mach.dose", snp_files, value = TRUE))
gen_dat = read.table(gzfile(dosage_file),header=F,stringsAsFactors=F,fill=T)
gen_dat.1 = gen_dat[complete.cases(gen_dat),]
sample.names.gen <- gen_dat.1[,1]
sample.names.gen.1 <- c()
sample.strsp <- strsplit(sample.names.gen,"->")


for(i in 1:length(sample.strsp)){
	name.tmp <- sample.strsp[[i]][1]
   	sample.names.gen.1 <- c(sample.names.gen.1,name.tmp)
}

gen_dat.1[,1] <- sample.names.gen.1
  

#----- find the sample intersection between genotype data and expr data -------
sample.names.inter <- intersect(sample.names.gen.1,sample.names)
x <- gen_dat.1[match(sample.names.inter, sample.names.gen.1),-c(1,2)]
y <- expr[match(sample.names.inter, sample.names),-1]
x.raw <- x
# remove cols in x with all NAs
x <- x[, colSums(is.na(x)) != nrow(x)]




#--------------------- standardize x matrix -----------------------------
x <- scale(x)
x <- x[, colSums(is.na(x)) != nrow(x)]

#---------------------- linear regression --------------------------------

p.val <- numeric()
eff <- numeric()
r.square <- numeric()

for (i in 1:dim(x)[2]){
	smry.lm <- summary(lm(y~x[,i]))
    p.val[i] <- smry.lm$coefficients[2,4]
    eff[i] <- smry.lm$coefficients[1,4]
    r.square[i] <- smry.lm$r.squared
}

sig_snp <- which(p.val < 0.05/dim(x)[2])

#--------------------- save result ------------------------
result_dir <- paste0("/ysm-gpfs/pi/zhao/wz284/trans_eqtl/result/linear_regression/",
	chr,"/")
dir.create(result_dir)
if (length(sig_snp) == 0){
	save(x.raw, x, y, p.val, eff, r.square, 
    	file=paste0(result_dir, gene, "_", tissue, "_lm.Rdata")) 	
} else {
	save(x.raw, x, y, p.val, eff, r.square, 
    	file=paste0(result_dir, gene, "_", tissue, "_signal_lm.Rdata")) 	
}







