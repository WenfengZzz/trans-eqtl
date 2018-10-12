position_files <- list.files("snp_locations/chr1/")
gene_name <- substr(position_files, 1, 15)


for (i in 1:length(position_files)){
    position.dir <- paste0("/home/wz284/scratch60/trans_eqtl/data/trans_snp/snp_position/")
    position_file <- paste0(position.dir, position_files[i])
    output.dir <- paste0("/home/wz284/scratch60/trans_eqtl/data/trans_snp/snp_info/", gene_name[i], "/")
    output.file <- paste0(output.dir, gene_name[i], ".vcf.gz")
    prefix <- paste0(output.dir, gene_name[i])
    sent <- paste0("mkdir -p ",output.dir,
    "; sh /home/wz284/scratch60/trans_eqtl/scripts/data_preprocessing/extract_geno_position.sh ",
      position_file, " " ,output.file, " ",prefix)
    cat(sent,"\n")    
}

file = paste0("extractsnp.sh")
sink(file)
sh.generate(gene_exp, gene_snp_inter)

