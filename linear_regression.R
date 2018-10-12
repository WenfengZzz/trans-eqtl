rm(list=ls())
setwd("~/Box Sync/Graduate/Research/Hi-C")

gene <- na.omit(gene_intra)
i <- 3
for (i in 1:length(gene)){
  snp_gene <- gene[i]
  # snp_gene <- "ENSG00000244674"
  result_file <- paste0("result/ENSG00000131697/", snp_gene, "_Liver_elnt.Rdata")
  if(!file.exists(result_file)){
    next()
  }
  load(result_file)
  x <- as.matrix(drop(x.array))
  y <- as.matrix(drop(y.array))
  
  p.val <- numeric()
  eff <- numeric()
  r.square <- numeric()
  
  for (j in 1:dim(x)[2]){
    smry.lm <- summary(lm(y~x[,j]))
    p.val[j] <- smry.lm$coefficients[2,4]
    eff[j] <- smry.lm$coefficients[1,4]
    r.square[j] <- smry.lm$r.squared
  }
  
  snp_dir <- paste0("snp_info/", snp_gene, "/")
  snp_files <- list.files(snp_dir)
  snp_file <- paste0(snp_dir, snp_files[grep("mach.info", snp_files)])
  if(!file.exists(snp_file)){
    next()
  }
  
  machinfo <- read.table(snp_file, header = TRUE, fill = TRUE)
  snp_loc <- machinfo[,1]
  text <- gsub("_[A-Z]_[A-Z]_[a-z][0-9][0-9]", "", snp_loc)
  text <- gsub("[0-9]_", "", text)
  snp_loc <- as.numeric(text)
  # summary(snp_loc)
  # length(snp_loc)
  # dim(x)
  snp_summary <-  cbind(snp_loc, p.val, r.square)
  summary_file <- paste0("result/linear_result/",snp_gene, "_linear_results.Rdata")
  save(snp_summary, file = summary_file)
}

summary(x[,1])
#manhattan plot
mplot <- paste0(files[j], "_man.png")

#png(mplot)
nlogp <- -log10(p.val)
t <- -log10(0.1/3183)
locpval <- cbind(snp_loc, nlogp)
plot(locpval, xlab = "SNP location (chr1)", ylab = "-log10(pvalue)", cex=0.5)
abline(h=t, lty=2)
abline(h=t-1, lty=2)
#dev.off()

## QQ-plot


library(lattice)
qqunif.plot<-function(pvalues, 
                      should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                      xlab=expression(paste("Expected (",-log[10], " p-value)")),
                      ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                      draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }
  
  #draw the plot
  xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points, 
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           panel.xyplot(x,y, ...);
           panel.abline(0,1);
         }, par.settings=par.settings, ...
  )
}


qplot <- paste0(files[j], "_qq.png")
png(qplot)
qqunif.plot(p.val)
dev.off()

##
which(p.val<0.1/3183)
p.val[1083]
log10(0.05/3183)
-log10(p.val[1083])
which(betas.elnt!=0)
r.square[1083]
# nlogp <- -log10(p.val)
# abline(h = t, lty=2)
# t <- -log10(1/3183)
# which(p.val < 1/3183)
# which(r.square>0.1)
# which(p.val<0.005)
# which(betas.elnt!=0)-1
# p.val[which(p.val<0.01)]


# pca <- prcomp(x, center = TRUE, scale = TRUE)
# pc <- pca$x
# for (i in 1:dim(x)[2]){
#   smry.lm <- summary(lm(y~x[,i]+pc[,1:10]))
#   p.val[i] <- smry.lm$coefficients[2,4]
#   eff[i] <- smry.lm$coefficients[1,4]
# }
# p.val[which(p.val<0.01)]



