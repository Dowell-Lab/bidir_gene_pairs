#!/bin/bash                                                                                                                                                      

#SBATCH --job-name=test_corr_filt_sig # Job name                                                                                                                           

#SBATCH --mail-type=FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)                                                                                         
#SBATCH --mail-user=rutendo.sigauke@ucdenver.edu                                                                                             
#SBATCH --nodes=1                                                                                                                                                
#SBATCH --ntasks=1 # Num of CPU                                                                                                  
#SBATCH --time=00:10:00 # Time limit                                                                                                                
#SBATCH --partition short                                                                                                                                        
#SBATCH --mem=5mb # Memory                                                                                                                                 
#SBATCH --output=/scratch/Users/rusi2317/projects/meta_analysis_qc/hg38/bin/bidir_gene_pairs/test/e_and_o/%x_%j.out                                                                                                           
#SBATCH --error=/scratch/Users/rusi2317/projects/meta_analysis_qc/hg38/bin/bidir_gene_pairs/test/e_and_o/%x_%j.err                                                                  

echo "BEGIN: Correlations..................."

#get sample paths
#working directory
WD=/scratch/Users/rusi2317/projects/meta_analysis_qc/hg38/bin/bidir_gene_pairs

#input to filtering script
FILT_CORR=${WD}/R/filter_significant_pairs.R 
OUT=${WD}/test/out_significant/
CORR=${WD}/test/out_tissues
CORR_ALL=${WD}/test/out_allsamples

ANNOTATIONS=${WD}/test/test_chroms.tsv

##################################
#run script

Rscript --vanilla ${FILT_CORR} -c ${ANNOTATIONS} -r ${CORR} -o ${OUT} -t
Rscript --vanilla ${FILT_CORR} -c ${ANNOTATIONS} -r ${CORR_ALL} -o ${OUT}

echo "DONE: Matrix pearson correlation......"
