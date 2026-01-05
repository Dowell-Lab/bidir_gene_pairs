# bidir_gene_pairs
Bidirectional transcript and gene pairs derived from nascent RNA data:

## Requirements

- parallel
- optparse version 1.7.3        
- data.table version 1.14.2     
- tidyr version 1.2.1          
- dplyr version 1.0.10         
- WGCNA version 1.70-3          

## Running in the command line

### `nascent_correlations.R`

Calculated correlations beatween gene and bidirectional transcription using all samples.

```
Rscript --vanilla nascent_correlations.R -h
Usage: nascent_correlations.R [options]

Options:
	-t CHARACTER, --tpms=CHARACTER
		path to TPM normalized counts

	-m CHARACTER, --samplemeta=CHARACTER
		path to metadata table for all samples

	-i CHARACTER, --chr_id=CHARACTER
		chromosome to process

	-c INTEGER, --ncores=INTEGER
		number of cores requisted (Note: more will speed up the run time) [default = 1]

	-w INTEGER, --window=INTEGER
		window in bases around TSS for bidirectionals to include [default = 1000000 bp]

	-n INTEGER, --nlimit=INTEGER
		minimum number of transcribed samples to include [default = 3]

	-u CHARACTER, --tissue=CHARACTER
		tissue to process

	-e CHARACTER, --exclude_missing_data=CHARACTER
		exclude observations where one of the samples has missing data (similar to use='pairwise.complete.obs') [default = FALSE]

	-o CHARACTER, --out=CHARACTER
		path to output directory [default = ./]

	-h, --help
		Show this help message and exit

```

### Example Output

The output is a form of a bed12 file where the first 6 columns are gene coodinates and the following 6 are bidirectional coordinates. Remaining columns are the summary statistics for correlation and the relationship between the gene and bidirectional.

- transcript1_chrom	   : Gene chromosome  
- transcript1_start	   : Gene start coordinate
- transcript1_stop         : Gene stop coordinate
- transcript_1	           : Gene id 
- transcript1_score	   : Gene score (. since none was assigned)
- transcript1_strand 	   : Gene strand
- transcript2_chrom 	   : Bidirectional chromosome
- transcript2_start 	   : Bidirectional start coordinate
- transcript2_stop 	   : Bidirectional stol coordinate
- transcript_2             : Bidirectional id
- transcript2_score 	   : Bidirectional score (i.e. the number of papers that support a bidirectional from muMerge)
- transcript2_strand 	   : Bidirectional strand (. since these are not stranded)
- pcc                      : Pearsons correlation coefficient
- pval                     : P-value
- adj_p_BH                 : Adjusted p-value (Benjamini-Hochberg correction)
- nObs                     : Number of observations in correlation analysis
- t                        : T statistic
- distance_tss 		   : Distance between the gene start (TSS) and the bidirectional start coordinate 
- distance_tes 		   : Distance between the gene stop (TES) and the bidirectional start coordinate
- position 		   : Is the bidirectional upstream or downstream of the TSS
- percent_transcribed_both : Percent of the number of observed samples used in analysis
