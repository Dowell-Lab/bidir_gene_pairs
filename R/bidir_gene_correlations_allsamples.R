#! /usr/bin/env Rscript

##################
##load packages ##
##################

suppressMessages(library(WGCNA)) ## faster cor()
suppressMessages(library(dplyr)) ## for the R pipes
suppressMessages(library(tidyr)) ## for tidying the dataframes
suppressMessages(library(data.table)) ## load files into R faster
suppressMessages(library(parallel)) ##running code in parallel
suppressMessages(library(optparse)) ##adding arguments 

#############################################
## Initialize command options for script   ##
#############################################
# define input and output options
option_list = list(
	make_option(c("-t", "--tpms"), type="character", default=NULL, 
	      help="path to TPM normalized counts", metavar="character"),
        make_option(c("-c", "--chr"), type="character", default="NULL",     
              help="chromosome file for processing", metavar="character"),
	make_option(c("-o", "--out"), type="character", default="./", 
              help="path to output directory [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$tpms)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input tpms).n", call.=FALSE)
}

#initalized paths and options for script
tpms_datatable <- data.table::fread(opt$tpms)
chroms <- data.table::fread(opt$chr)
nsamples <- ncol(tpms_datatable)
output_folder <- opt$out

##Start run Log
sink(paste0(output_folder,"log_all_",Sys.Date(),".txt"))

cat("Running co-transcription analyses")
cat("\n")

cat(paste0("Date: ", Sys.Date()))
cat("\n")

print("START")
cat("\n")

print(paste("Normalized Counts :", as.character(opt$tpms)))
cat("\n")

print(paste("Chromosome File   :", as.character(opt$chr)))
cat("\n")

print(paste("Output Directory  :", as.character(opt$out)))
cat("\n")

print(paste("Number of Samples :", as.character(nsamples)))
cat("\n")

##long format of pairwise matrix matrix 
restructure_cor_matrix <- function(matrix){
    
    #' get long format for matrix output from WGCNAs corAndPvalue() function
    #' 
    #' @description Summarize the matrix in 'transcript1, transcript2, value' 
    #' format
    #' 
    #' @param matrix : pairwise comparisons in matrix format
    #'
    #'
    #' @usage restructure_cor_matrix(matrix)
    #' @return A tibble with values for pairs from matrix
    #' @export
    
    #convert to dataframe
    matrix_df <- as.data.frame(matrix)
    matrix_df$transcript_1 <- colnames(matrix_df)
    
    #change the structure of the dataframe to 3 column dataframe
    #> transcript_1, transcript_2, coefficient
    matrix_df_long <- matrix_df %>% tidyr::pivot_longer(!transcript_1, 
                                                        names_to = "transcript_2", 
                                                        values_to = "value",
                                                        values_drop_na = TRUE)
    
    #remove same pair correlations
    matrix_df_unique <- subset(matrix_df_long, 
                               transcript_1 != transcript_2)
    
    #remove NA comparisons (these are transcripts not in all samples)
    matrix_df_noNA <- subset(matrix_df_unique, 
                             !is.na(matrix_df_unique[[3]]))
    
    matrix_df_noNA$pair_id <- paste0(matrix_df_noNA$transcript_1,':',
                                     matrix_df_noNA$transcript_2)
    
    return(matrix_df_noNA)
    
}

##combine correlation outputs 
cor_summary_stats <- function(corAndPvalueOut_list) {
    
    #' get summary from corAndPvalue() output tibble list 
    #' takes output from restructure_cor_matrix() function
    #'
    #' @description All corAndPvalue() statistics in a single tibble 
    #' 
    #' 
    #' @param corAndPvalueOut_list: list of matrix for pairwise comparisons with summary stats
    #'
    #'
    #' @usage cor_summary_stats(list_of_tibbles)
    #' @return A tibble with values for pairs from matrix list
    #' @export
    
    ##summary table with all the statistics
    ##combined as shown below
    corAndPvalueOut_tibble <- purrr::reduce(list(corAndPvalueOut_list$cor,
                                              corAndPvalueOut_list$p[c(3,4)],
                                              corAndPvalueOut_list$t[c(3,4)],
                                              corAndPvalueOut_list$nObs[c(3,4)]),
                                         dplyr::left_join, by = 'pair_id')
    ##renaming column names
    corAndPvalueOut_tibble <- corAndPvalueOut_tibble %>% dplyr::rename("pcc" ="value.x",
                                                             "pval" ="value.y",
                                                              "t"="value.x.x",
                                                              "nObs"="value.y.y")
    
    ##calculating the adjusted p-values
    corAndPvalueOut_tibble$adj_p_BH <- p.adjust(corAndPvalueOut_tibble$pval,
                                                      method = 'BH')
    
    ##rearrange columns
    corAndPvalueOut_tibble <- corAndPvalueOut_tibble %>% dplyr::relocate(pcc,
                                                                         .after = pair_id)
    corAndPvalueOut_tibble <- corAndPvalueOut_tibble %>% dplyr::relocate(adj_p_BH,
                                                                         .after = pcc)
    return(corAndPvalueOut_tibble)   
    
}

##get metadata for pairs 
cor_pair_metadata <- function(tpm_filtered_chrm, corAndPvalueOut_tibble) {

    #' get metadata for all pairs
    #'
    #' @description Take tpm input file with metadata and the single tibble 
    #' file with correlations summary stats and combines to give a unified
    #' table with all information
    #' 
    #' @param tpm_filtered_chrm : TPMs in a BED6 format
    #'
    #' @param corAndPvalueOut_tibble : tibble file with corr and p-values
    #'
    #'
    #' @usage cor_summary_stats(list_of_tibbles)
    #' @return A tibble with values for pairs from matrix list
    #' @export

    transcript_coords <- tpm_filtered_chrm[,1:6]
    colnames(transcript_coords) <- c("chrom","start",
                                     "stop","gene_transcript",
                                     "score","strand")
       
    transcript_coords$transcript_type <- ifelse(grepl("chr*",
                                                      transcript_coords$gene_transcript),
                                                      "Bidirectional", "Gene")

    # subset transcript 1 coordinate data
    transcript1_choords <- transcript_coords[transcript_coords$gene_transcript %in% corAndPvalueOut_tibble$transcript_1,]
    colnames(transcript1_choords) <- c("transcript1_chrom","transcript1_start",
                                       "transcript1_stop",
                                       "gene_transcript",
                                       "transcript1_score",
                                       "transcript1_strand",
                                       "transcript1_type")

    # subset transcript 2 coordinate data
    transcript2_choords <- transcript_coords[transcript_coords$gene_transcript %in% corAndPvalueOut_tibble$transcript_2,]
    colnames(transcript2_choords) <- c("transcript2_chrom","transcript2_start",
                                       "transcript2_stop",
                                       "gene_transcript",
                                       "transcript2_score",
                                       "transcript2_strand",
                                       "transcript2_type")

    # combine transcript coordinates with correlations
    corAndPvalueOut_transcript1 <- dplyr::left_join(corAndPvalueOut_tibble,
                                             transcript1_choords,
                                             by = c("transcript_1"="gene_transcript"))

    corAndPvalueOut_transcript1and2 <- dplyr::left_join(corAndPvalueOut_transcript1,
                                                transcript2_choords,
                                                by = c("transcript_2"="gene_transcript"))

    # remove redundant pairs and only report unique Gene-Bidirectional correlations
    # now transcript_1 are Genes and transcript_2 are Bidirectionals
    corAndPvalueOut_gene_bidirs <- subset(corAndPvalueOut_transcript1and2,
                                      transcript1_type == 'Gene' &
                                      transcript2_type != 'Gene')

    # Now calculating the distance between genes and bidirectionals
    # relative center position of the bidirectional transcript
    bidir_center_pos <- (corAndPvalueOut_gene_bidirs$transcript2_stop - corAndPvalueOut_gene_bidirs$transcript2_start)/2

    # get genomic center position
    bidir_center <- round(bidir_center_pos, digits = 0) + corAndPvalueOut_gene_bidirs$transcript2_start

    # the distance calculation is from the center position of bidirectional to the start of the gene
    corAndPvalueOut_gene_bidirs$distance_tss <- bidir_center - corAndPvalueOut_gene_bidirs$transcript1_start
    corAndPvalueOut_gene_bidirs$distance_tes <- bidir_center - corAndPvalueOut_gene_bidirs$transcript1_stop
    
    #reordering data.frame as a bed 12 plus PCC output and distances
    column_order <- c("transcript1_chrom","transcript1_start","transcript1_stop",
                     "transcript_1", "transcript1_score", "transcript1_strand",
                     "transcript2_chrom","transcript2_start","transcript2_stop",
                     "transcript_2", "transcript2_score", "transcript2_strand",
                     "pcc","pval","adj_p_BH","nObs","t","distance_tss","distance_tes")
    
    corAndPvalueOut_gene_bidirs_bedformat <- corAndPvalueOut_gene_bidirs[, column_order]
    
    #annotate position of bidirectional relative to gene 
    corAndPvalueOut_gene_bidirs_bedformat$position <- ifelse(corAndPvalueOut_gene_bidirs_bedformat$transcript1_strand=='+' & 
                                                             corAndPvalueOut_gene_bidirs_bedformat$distance_tss<0 |
                                                             corAndPvalueOut_gene_bidirs_bedformat$transcript1_strand=='-' & 
                                                             corAndPvalueOut_gene_bidirs_bedformat$distance_tss>0, "upstream","downstream")

    return(corAndPvalueOut_gene_bidirs_bedformat)

}

##############################################
##matrix correlation by chromosome function ##
##############################################
transcript_pearsons_by_chromosome <- function(tpms_datatable, chromosome, nsamples){
    
    #' calculate pearson's correlations for all transcripts in input
    #' 
    #' @description This function will calculate person's R and significance for 
    #' input normalized counts  
    #' 
    #' @param tpms_datatable path i.e. path to normalized counts
    #'
    #' @param chromosome id based on the input list of chromosomes 
    #'
    #' @param output_folder output directory
    #'
    #' @usage transcript_pearsons_by_chromosome(tpms_datatable, chromosome, output_folder)
    #' @return A data.frame with all pairwise correlations and significance
    #' @export
    
    #get a subset of genes and bidirs by chromosome id
    tpms_chrm <- as.data.frame(subset(tpms_datatable, chrom == chromosome))
    
    #transpose the dataframe (ie. rows are samples and columns are transcripts)
    #this will also convert dataframe to matrix
    tpms_chrm_t <- t(tpms_chrm[c(7:nsamples)])
    colnames(tpms_chrm_t) <- tpms_chrm$gene_transcript
    
    #log transform the normalized tpm counts(base 10)
    ##NOTE: running with adding 1s and converting 0s to NA
    #tpms_chrm_log10 <- log(tpms_chrm_t, base=10) ##original running process
    tpms_chrm_log10 <- log(tpms_chrm_t+1, base=10)

    ##NOTE: running correlations on pairs BUT ignoring the 
    ##samples with missing data
    tpms_chrm_log10_NAs <- dplyr::na_if(tpms_chrm_log10, 0)

    ########################################################
    ##Using WGCNA calculate correlations and p-values for ##
    ##relavant samples with trancriptio                   ##
    ########################################################
    ##calculated pcc, pvalue, z stat, t stat and number of observations (i.e.)
    ##samples with both genes and bidirectionals transcribed         
    corAndPvalueOut <- WGCNA::corAndPvalue(tpms_chrm_log10_NAs, 
                                           use="pairwise.complete.obs")
    
    ##restructure all the matrix outputs long formats 
    corAndPvalueOut_matrix_list <- lapply(corAndPvalueOut, restructure_cor_matrix)
    
    ##combine all summary stats all in one tibble
    corAndPvalueOut_all_tibble <- cor_summary_stats(corAndPvalueOut_matrix_list)
    
    ##remove redundant pairs and add metadata
    ##bidir and gene pair ids, distances 
    corAndPvalueOut_pairs <- cor_pair_metadata(tpms_chrm, corAndPvalueOut_all_tibble)

    ##add sample summary of counts
    corAndPvalueOut_pairs$tissue <- "All_samples"
    corAndPvalueOut_pairs$percent_transcribed_both <- (corAndPvalueOut_pairs$nObs/(nsamples-6))*100
    
    
    #saving the final dataframes
    #including date and chromosome ids to the file names
    outdir <- output_folder 
    
    final_path <- paste0(outdir,
                         'pearson_correlation_', 
                         chromosome,
                         '.tsv.gz' )
    
    data.table::fwrite(corAndPvalueOut_pairs,
                       final_path,
                       sep='\t')
}

chromosome_list <- as.character(chroms$V1) 

parallel::mclapply(chromosome_list, transcript_pearsons_by_chromosome,
		   tpms_datatable = tpms_datatable, nsamples = nsamples, 
		   mc.cores = length(chromosome_list))

cat("Session Summary")
cat("\n")

print(sessionInfo())
cat("\n")

cat("DONE!")
cat("\n")

sink()
