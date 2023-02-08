#!/bin/bash                                                                                                                                                      

#SBATCH --job-name=test_corr_tissues # Job name                                                                                                                           

#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)                                                                                         
#SBATCH --mail-user=rutendo.sigauke@ucdenver.edu                                                                                             
#SBATCH --nodes=1                                                                                                                                                
#SBATCH --ntasks=1 # Num of CPU                                                                                                  
#SBATCH --time=00:05:00 # Time limit                                                                                                                
#SBATCH --partition short                                                                                                                                        
#SBATCH --mem=2mb # Memory                                                                                                                                 
#SBATCH --output=/scratch/Users/rusi2317/projects/meta_analysis_qc/hg38/scripts/bidir_gene_pairs/test/e_and_o/%x_%j.out                                                                                                           
#SBATCH --error=/scratch/Users/rusi2317/projects/meta_analysis_qc/hg38/scripts/bidir_gene_pairs/test/e_and_o/%x_%j.err                                                                  

echo "BEGIN: Correlations..................."
WD=/scratch/Users/rusi2317/projects/meta_analysis_qc/hg38/scripts/bidir_gene_pairs
WGCNA_CORR=${WD}/R/bidir_gene_correlations_tissues.R
TPMS=${WD}/test/test_tpms.tsv 
CHR=${WD}/test/test_chroms.tsv
OUT=${WD}/test/out_tissues/
META=${WD}/test/230131_db_output.tsv 

#run script
Rscript --vanilla ${WGCNA_CORR} -t ${TPMS} -m ${META} -c ${CHR} -o ${OUT}


echo "DONE: Matrix pearson correlation......"
