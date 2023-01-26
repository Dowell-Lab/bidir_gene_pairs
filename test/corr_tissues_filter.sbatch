#!/bin/bash                                                                                                                                                      

#SBATCH --job-name=test_corr_tissues_filt_sig # Job name                                                                                                                           

#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)                                                                                         
#SBATCH --mail-user=rutendo.sigauke@ucdenver.edu                                                                                             
#SBATCH --nodes=1                                                                                                                                                
#SBATCH --ntasks=1 # Num of CPU                                                                                                  
#SBATCH --time=00:10:00 # Time limit                                                                                                                
#SBATCH --partition short                                                                                                                                        
#SBATCH --mem=1gb # Memory                                                                                                                                 
#SBATCH --output=/scratch/Users/rusi2317/e_and_o/%x_%j.out                                                                                                           
#SBATCH --error=/scratch/Users/rusi2317/e_and_o/%x_%j.err                                                                  

echo "BEGIN: Correlations..................."

#get sample paths
#working directory
WD=/scratch/Users/rusi2317/projects/meta_analysis_qc/hg38/scripts/bidir_gene_pairs

#input to filtering script
FILT_CORR=${WD}/R/filter_pairs_by_tissues.R
TPMS=${WD}/test/test_tpms_samples.tsv
OUT=${WD}/test/
META=${WD}/test/20211026_db_output_tissues.csv
CORR=${WD}/test/out_tissues

ANNOTATIONS=${WD}/test/test_chroms.tsv

##################################
#run script

Rscript --vanilla ${FILT_CORR} -t ${TPMS} -m ${META} -c ${ANNOTATIONS} -r ${CORR} -o ${OUT}

echo "DONE: Matrix pearson correlation......"
