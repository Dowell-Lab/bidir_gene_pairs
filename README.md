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

See `test` folder to see example input and output.

### `bidir_gene_correlations_allsamples.R`

Calculated correlations beatween gene and bidirectional transcription using all samples.

```
Rscript --vanilla bidir_gene_correlations_allsamples.R -h
Usage: bidir_gene_correlations_allsamples.R [options]


Options:
	-t CHARACTER, --tpms=CHARACTER
		path to TPM normalized counts

	-c CHARACTER, --chr=CHARACTER
		chromosome file for processing

	-o CHARACTER, --out=CHARACTER
		path to output directory [default= ./]

	-h, --help
		Show this help message and exit
```

### `bidir_gene_correlations_tissues.R`

Calculated correlations beatween gene and bidirectional transcription using specific *Tissues*.

The current limit for the number of samples per tissue used is 15 samples.

```
Rscript --vanilla bidir_gene_correlations_tissues.R -h
Usage: bidir_gene_correlations_tissues.R [options]


Options:
	-t CHARACTER, --tpms=CHARACTER
		path to TPM normalized counts

	-m CHARACTER, --samplemeta=CHARACTER
		path to metadata table for all samples

	-c CHARACTER, --chr=CHARACTER
		chromosome file for processing

	-o CHARACTER, --out=CHARACTER
		path to output directory [default= ./]

	-h, --help
		Show this help message and exit
```


### `filter_significant_pairs.R`

Takes input from `bidir_gene_correlations_tissues.R` and `bidir_gene_correlations_allsamples.R`. Filtering is done based on:

1. Adjusted p-value : 0.01
2. Pearson correlation coefficient (PCC) : >/< 0.6
3. Percent of samples with transcription : > 5%

```
Rscript --vanilla filter_significant_pairs.R -h
Usage: filter_significant_pairs.R [options]

Options:
	-r CHARACTER, --correlations=CHARACTER
		path to chromosome level correnations from bidir_gene_correlations_allsamples.R or  bidir_gene_correlations_tissues.R

	-c CHARACTER, --chr=CHARACTER
		chromosome file for processing

	-d INTEGER, --min_dist=INTEGER
		minimum distance between pairs [default = 1000000 bases]

	-e INTEGER, --percent_trans=INTEGER
		percent of sample with transcription for transcript [default = 5]

	-v DOUBLE, --rvalue=DOUBLE
		pearson's R value cut-off [default = 0.6 ]

	-p DOUBLE, --adj_pvalue=DOUBLE
		adjusted p-value filter for called pairs [default = 0.01 ]

	-o CHARACTER, --out=CHARACTER
		path to output directory [default= ./]

	-t, --tissue
		Are the correlations based on a per-tissue and per-chromosome basis? [default = FALSE]

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
- tissue 		   : Tissue id based on metadata for tissue derived correlations (labeled `All_samples` if all samples are used)
- percent_transcribed_both : Percent of the number of observed samples used in analysis