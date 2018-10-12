position_file=$1
#position file format:/ysm-gpfs/pi/zhao/from_louise/yh367/GTEX/cis_snp_by_gene/chr1/ENSG00000273483.1/union
output_file=$2
#ouput_file example:/ysm-gpfs/pi/zhao/wl382/data/GEUV/cis_snp_by_gene/chr17/ENSG00000196544.6/ENSG00000196544.6.vcf.gz
output_prefix=$3
#output_prefix example:ENSG00000196544.6_
vcftools --gzvcf /ysm-gpfs/scratch60/wl382/GTEX/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --positions ${position_file} --recode -c | gzip -c > ${output_file} ; DosageConvertor --vcfDose ${output_file} --prefix ${output_prefix} --type mach --format DS
