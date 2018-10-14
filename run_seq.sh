#!/usr/bin bash
file_nums=(0095 0068 0053 0041 0048 0051 0050 0039 0038 0039 0055)
j=1
for i in {0110..0192}; do
       module load SimpleQueue/3.1
       sqCreateScript -q scavenge -m 10000 -N ${i} -n 20 -w 48:00:00 prediction_4models_vb.annot.sh_${i} > run_${i}.sh
	sed '2 i#SBATCH --constraint=avx2' -i run_${i}.sh
#	module load dSQ
#	dSQ --jobfile /ysm-gpfs/pi/zhao/wl382/snpPred_epi/Annotation/scripts/expr_annot/cross_validation/sh_file/cv_logit_all_wl_chr${j}_sh_${i}  --mem-per-cpu=6g -t 48:00:00 -n 5 -p scavenge -J ld_a -C avx2 --mail-user=123@qq.com > run_chr${j}_${i}.sh
        sbatch run_${i}.sh
	#rm /ysm-gpfs/pi/zhao/wl382/snpPred_epi/Annotation/scripts/results_dealing/sh_files/model_logitchr1_sh_${i}
done
