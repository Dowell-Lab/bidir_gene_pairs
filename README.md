# bidir_gene_pairs
Bidirectional transcript and gene pairs derived from nascent RNA data

See `test` folder to see example input and output.

Scripts:

- `bidir_gene_correlations_allsamples.R`

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

- `bidir_gene_correlations_tissues.R`

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


- `filter_pairs_by_tissues.R`

Takes input from `bidir_gene_correlations_tissues.R` abd filters pairs by chromosomes. Filtering is done based on:

1. Adjusted p-value : 0.01
2. Pearson correlation coefficient (PCC) : >/< 0.9
3. Percent of samples with transcription : > 75%

```
Rscript --vanilla filter_pairs_by_tissues.R -h
Usage: filter_pairs_by_tissues.R [options]


Options:
	-t CHARACTER, --tpms=CHARACTER
	   path to TPM normalized counts

	   -m CHARACTER, --samplemeta=CHARACTER
	      path to metadata table for all samples

	      -r CHARACTER, --correlations=CHARACTER
	      	 path to chromosome level correnations from nascent_transcripts_pearsons_corr.R

		 -c CHARACTER, --chr=CHARACTER
		    chromosome file for processing

		    -d INTEGER, --min_dist=INTEGER
		       minimum distance between pairs [default = 1000000 bases]

		       -e INTEGER, --percent_trans=INTEGER
		       	  percent of sample with transcription for transcript [default = 25]

			  -v DOUBLE, --rvalue=DOUBLE
			     pearson's R value cut-off [default = 0.9 ]

			     -p DOUBLE, --adj_pvalue=DOUBLE
			     	adjusted p-value filter for called pairs [default = 0.001 ]

				-o CHARACTER, --out=CHARACTER
				   path to output directory [default= ./]

				   -h, --help
				       Show this help message and exit


```