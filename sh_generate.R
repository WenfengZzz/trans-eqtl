args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
bin.dir <- paste0("/home/wz284/scratch60/trans_eqtl/data/trans_snp/exp_gene_snp_region/",chr,"/")
snp_bin_files <- list.files(bin.dir)
for (i in 1:length(snp_bin_files)){
 	sent <- paste0("module load R; Rscript /home/wz284/scratch60/trans_eqtl/scripts/data_preprocessing/bin_to_position.R ",chr," ", snp_bin_files[i])
	cat(sent,"\n")    
}
