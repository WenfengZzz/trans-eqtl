# processing results from linear regression
# intra
result_dir <- "/Users/wayne/Box Sync/Graduate/Research/Hi-C/result/linear_result/131697/intra/"
result_files <- list.files(result_dir)
result <- list()
for (i in 1:length(result_files)){
  result_file <- paste0(result_dir, result_files[i])
  load(result_file)
  result[[i]] <- snp_summary
}

result[[1]]
merged_result <- as.matrix(do.call(rbind, result))
dim(merged_result)
head(merged_result)

unique_result <- unique(merged_result)

head(unique_result)
dim(merged_result)
dim(unique_result)

# range
chr1_range <- chr1_intra_region[2:3,]
inrange <- rep(FALSE, dim(unique_result)[1])
i <- 1
# a <- 31200342
test <- FALSE
for (j in 1:dim(unique_result)[1]){
  for (i in 1:dim(chr1_range)[2]) {
    inrange[j] <- inrange[j] | (chr1_range[1,i]<unique_result[j,1] & unique_result[j,1]<chr1_range[2,i])
    #test <- test | (chr1_range[1,i]< a & a <chr1_range[2,i])
  }
}

sum(inrange)
which(inrange == TRUE)
summary(unique_result)

sum(inrange)
inrange_result <- unique_result[inrange, ]
dim(inrange_result)
dim(inrange_result_chr3)
# m-plot
dev.off()
nlogp <- -log10(inrange_result[,2])
t <- -log10(0.05/length(nlogp))
locpval <- cbind(inrange_result[,1], nlogp)
plot(locpval, xlab = "SNP location (chr1)", ylab = "-log10(pvalue)", cex=0.5, xlim = c(0, 2.5e+08))
abline(v=5e+06, col="red", lwd = 2)
#abline(h=t-1, lty=2)

qqunif.plot(na.omit(inrange_result[,2]))

####################################################
####################################################
#inter

# chr9 ENSG00000196873_Liver_elnt.Rdata
# chr3 ENSG00000244674_Liver_elnt.Rdata
# chr16 ENSG00000260518_Liver_elnt.Rdata


result_dir <- "/Users/wayne/Box Sync/Graduate/Research/Hi-C/result/linear_result/131697/inter/"
result_files <- list.files(result_dir)
result <- list()
for (i in 1:length(result_files)){
  i <- 1
  result_file <- paste0(result_dir, result_files[i])
  load(result_file)
  result[[i]] <- snp_summary
}

result[[1]]
merged_result <- as.matrix(do.call(rbind, result))
dim(merged_result)
head(merged_result)

unique_result <- unique(merged_result)

head(unique_result)
dim(merged_result)
dim(unique_result)

# range
inter_chr = 3
chr1_range <- as.matrix(chr1_inter_region[2:3,chr1_inter_region[1,]==inter_chr])
inrange <- rep(FALSE, dim(unique_result)[1])
i <- 1
# a <- 31200342
# test <- FALSE
for (j in 1:dim(unique_result)[1]){
  for (i in 1:dim(chr1_range)[2]) {
    inrange[j] <- inrange[j] | (chr1_range[1,i]<unique_result[j,1] & unique_result[j,1]<chr1_range[2,i])
    #test <- test | (chr1_range[1,i]< a & a <chr1_range[2,i])
  }
}
sum(inrange)
which(inrange == TRUE)
summary(unique_result)

sum(inrange)
inrange_result <- unique_result[inrange, ]
dim(inrange_result)

inrange_result_chr3 <- inrange_result

# m-plot
dev.off()
nlogp <- -log10(inrange_result[,2])
t <- -log10(0.05/length(nlogp))
locpval <- cbind(inrange_result[,1], nlogp)
plot(locpval, xlab = "SNP location (chr3)", ylab = "-log10(pvalue)", cex=0.5, ylim = c(0,4)) #xlim = c(0, 2.5e+08))
abline(v=5e+06, col="red", lwd = 2)
#abline(h=t-1, lty=2)


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


