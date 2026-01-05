#! /usr/bin/env Rscript

##################
##load packages ##
##################

suppressMessages(library(WGCNA)) ## faster cor()
suppressMessages(library(dplyr)) ## for the R pipes
suppressMessages(library(tidyr)) ## for tidying the dataframes
suppressMessages(library(data.table)) ## load files into R faster
suppressMessages(library(parallel)) ## running code in parallel
suppressMessages(library(optparse)) ## adding arguments 
suppressMessages(library(reshape2)) ## restructure matrix

#############################################
## Initialize command options for script   ##
#############################################
# define input and output options
option_list = list(
    make_option(c("-t", "--tpms"), type="character", default=NULL,
                help="path to TPM normalized counts", metavar="character"),
    make_option(c("-m", "--samplemeta"), type="character", default=FALSE,
                help="path to metadata table for all samples", metavar="character"),
    make_option(c("-i", "--chr_id"), type="character", default="chrY",
                help="chromosome to process", metavar="character"),
    make_option(c("-c", "--ncores"), type="integer", default=1,
                help="number of cores requisted (Note: more will speed up the run time) [default = %default]",
                metavar="integer"),
    make_option(c("-w", "--window"), type="integer", default=1000000,
                help="window in bases around TSS for bidirectionals to include [default = %default bp]",
                metavar="integer"),
    make_option(c("-n", "--nlimit"), type="integer", default=3,
                help="minimum number of transcribed samples to include [default = %default]",
                metavar="integer"),
    make_option(c("-u", "--tissue"), type="character", default="blood",
                help="tissue to process", metavar="character"),
    make_option(c("-e", "--exclude_missing_data"), type="character", default="FALSE",
                help="exclude observations where one of the samples has missing data (similar to use='pairwise.complete.obs') [default = %default]", 
                metavar="character"),
    make_option(c("-o", "--out"), type="character", default="./",
                help="path to output directory [default = %default]", metavar="character")
);
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$tpms)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input tpms).n", call.=FALSE)
}

###########################################
##Initialize variables and paths         ##
###########################################
#variables
ncores <- opt$ncores
output_folder <- opt$out
chromosome_id <- opt$chr_id
window <- opt$window
tissue <- opt$tissue
nlimit <- opt$nlimit

#files and paths
tpms_datatable <- data.table::fread(opt$tpms, nThread=1)

#update names for the bed6 counts 
bed6_colnames <- c("chrom","start","stop","gene_transcript","score","strand",
                     colnames(tpms_datatable)[7:ncol(tpms_datatable)])
colnames(tpms_datatable) <- bed6_colnames


if (opt$samplemeta == FALSE){
	metadata <- FALSE 
	}else{	
	metadata <- data.table::fread(opt$samplemeta, nThread=1)
 }

###########################################
##Processing functions                   ##
###########################################

##-----------------------------------------
## Get transcripts in a window
##-----------------------------------------

get_transcripts_in_window <- function(gene_name, gene_tpms_df, window = window){
    
    #' get bidirectional transcripts with a genes window
    #'
    #' @description takes in the gene name as character and data.table or 
    #' data.frame with counts and takes a specified window 
    #' 
    #' @param gene_name : character
    #'
    #' @param gene_tpms_df : data.frame or data.table with counts in BED6 formats
    #'
    #' @param window : interger with window size (default = 1000000)
    #'
    #' @usage get_transcripts_in_window(gene_name, gene_tpms_df, window = 1000000)
    #' @return A counts data.table with specified gene and bidirectional regions
    #' within the window
    #' @export
    
    ##get transcripts including genes within the 1MB window gene plus a 5kb pad
    window <- window + 5000
    
    ##filter gene of interest
    #gene_counts <- gene_tpms_df[gene_tpms_df$gene_transcript %in% gene_name,] #gene_tpms_df[grepl(gene_name, gene_tpms_df$gene_transcript),]
    gene_counts <- subset(gene_tpms_df, gene_transcript == gene_name)
    gene_chrom <- gene_counts$chrom
    gene_start <- gene_counts$start
    gene_stop <- gene_counts$stop
    gene_strand <- gene_counts$strand
    
    if (gene_strand == "+"){
    
        #get bidirectional transcripts in the specified window
        gene_bidir_window_tpms_df <- subset(gene_tpms_df,
                                     chrom == gene_chrom &
                                     start > gene_start - window & 
                                     start < gene_start + window)
        } else {
        
        #get bidirectional transcripts in the specified window
        gene_bidir_window_tpms_df <- subset(gene_tpms_df,
                                     chrom == gene_chrom &
                                     stop > gene_stop - window & 
                                     stop < gene_stop + window)
    }

    bidir_window_tpms_df <- gene_bidir_window_tpms_df[grepl("chr",gene_bidir_window_tpms_df$gene_transcript),]

    gene_bidir_tpms_df <- rbind(gene_counts, bidir_window_tpms_df)
    
    return(gene_bidir_tpms_df)
}

##-----------------------------------------
## Restructure correlation matrix
##-----------------------------------------

##get long format for matrix
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

##-----------------------------------------
## Get correlations summary stats
##-----------------------------------------

##get and combine summary statistics
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

##-----------------------------------------
## Get pair metadata
##-----------------------------------------

##add trascript metadata to pairs
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
    transcript1_choords <- transcript_coords[transcript_coords$gene_transcript %in% 
                                             corAndPvalueOut_tibble$transcript_1,]
    colnames(transcript1_choords) <- c("transcript1_chrom","transcript1_start",
                                       "transcript1_stop",
                                       "gene_transcript",
                                       "transcript1_score",
                                       "transcript1_strand",
                                       "transcript1_type")

    # subset transcript 2 coordinate data
    transcript2_choords <- transcript_coords[transcript_coords$gene_transcript %in% 
                                             corAndPvalueOut_tibble$transcript_2,]
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
    # here fixing the distance so that it is strand specific and relative to the genes start    
    corAndPvalueOut_gene_bidirs$distance_tss <- ifelse(corAndPvalueOut_gene_bidirs$transcript1_strand=='+',
                                                  bidir_center - corAndPvalueOut_gene_bidirs$transcript1_start,
                                                  corAndPvalueOut_gene_bidirs$transcript1_stop - bidir_center)

    corAndPvalueOut_gene_bidirs$distance_tes <- ifelse(corAndPvalueOut_gene_bidirs$transcript1_strand=='+',
                                                  bidir_center - corAndPvalueOut_gene_bidirs$transcript1_stop,
                                                  corAndPvalueOut_gene_bidirs$transcript1_start - bidir_center)


    #reordering data.frame as a bed 12 plus PCC output and distances
    column_order <- c("transcript1_chrom","transcript1_start","transcript1_stop",
                     "transcript_1", "transcript1_score", "transcript1_strand",
                     "transcript2_chrom","transcript2_start","transcript2_stop",
                     "transcript_2", "transcript2_score", "transcript2_strand",
                     "pcc","pval","adj_p_BH","nObs","t","distance_tss","distance_tes")
    
    corAndPvalueOut_gene_bidirs_bedformat <- corAndPvalueOut_gene_bidirs[, column_order]
    
    #annotate position of bidirectional relative to gene
    corAndPvalueOut_gene_bidirs_bedformat$position <- ifelse(corAndPvalueOut_gene_bidirs_bedformat$distance_tss < 0,
                                                             "upstream","downstream")

    return(corAndPvalueOut_gene_bidirs_bedformat)

}

##-----------------------------------------
## Process correlations
##-----------------------------------------

##put it all together and calculate correlations for tissue specific
transcript_pearsons_by_chromosome_tissue <- function(tpms_datatable, metadata, chromosome_id, tissue_name){
    
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

    nsamples <- ncol(tpms_datatable)
    ##get metadata for the samples analyzed
    sample_ids <- colnames(tpms_datatable[,7:nsamples])
    metadata_analyzed <- metadata[metadata$sample_name %in% sample_ids,]

    # get metadata for specific tissue of interest
    metadata_tissue <- subset(metadata_analyzed, tissue == tissue_name)
    #print(paste0("Tissue metadata ",tissue_name," : ", as.character(nrow(metadata_tissue))))
    #get a subset of genes and bidirs by chromosome id
    tpms_chrm <- subset(tpms_datatable, chrom == chromosome_id)
    
    #filter samples that match the tissue of interest
    tpms_chrms_tissue <- t(tpms_chrm[ ,colnames(tpms_chrm) %in% metadata_tissue$sample_name, with=FALSE]) 
    colnames(tpms_chrms_tissue) <- tpms_chrm$gene_transcript

    # log transform the matrix of tpms
    tpms_chrms_tissue_log10 <- log(tpms_chrms_tissue+1, base=10)

    # make sure that samples with 0 counts are excluded from the log() transformation
    #log transform the normalized tpm counts(base 10)
    ##NOTE: running with adding 1s and converting 0s to NA
     
    tpms_chrms_tissue_log10_NAs <- tpms_chrms_tissue_log10

    tpms_chrms_tissue_log10_NAs[tpms_chrms_tissue_log10_NAs == 0] <- NA
    #tpms_chrms_tissue_log10_NAs <- dplyr::na_if(tpms_chrms_tissue_log10, 0)  

    ########################################################
    ##Using WGCNA calculate correlations and p-values for ##
    ##relavant samples with trancriptio                   ##
    ########################################################
    ##calculated pcc, pvalue, z stat, t stat and number of observations (i.e.)
    ##samples with both genes and bidirectionals transcribed         
    corAndPvalueOut <- WGCNA::corAndPvalue(tpms_chrms_tissue_log10_NAs, 
                                           use="pairwise.complete.obs")
    
    ##restructure all the matrix outputs long formats 
    corAndPvalueOut_matrix_list <- lapply(corAndPvalueOut, restructure_cor_matrix)
    
    ##combine all summary stats all in one tibble
    corAndPvalueOut_all_tibble <- cor_summary_stats(corAndPvalueOut_matrix_list)
    #print(dim(corAndPvalueOut_all_tibble))
    #print(head(corAndPvalueOut_all_tibble))
    
    ##remove redundant pairs and add metadata
    ##bidir and gene pair ids, distances 
    corAndPvalueOut_pairs <- cor_pair_metadata(tpms_chrm, 
                                               corAndPvalueOut_all_tibble)
    
    ##add tissue summary of counts
    corAndPvalueOut_pairs$tissue <- tissue_name
    corAndPvalueOut_pairs$percent_transcribed_both <- (corAndPvalueOut_pairs$nObs/nrow(metadata_tissue))*100
    
    return(corAndPvalueOut_pairs)

    }

## calculate correlations for all samples
transcript_pearsons_by_chromosome <- function(tpms_datatable, chromosome_id){
    
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

    nsamples <- ncol(tpms_datatable)
    ##get metadata for the samples analyzed
    sample_ids <- colnames(tpms_datatable[,7:nsamples])
    #metadata_analyzed <- metadata[metadata$sample_name %in% sample_ids,]

    # get metadata for specific tissue of interest
    #metadata_tissue <- subset(metadata_analyzed, tissue == tissue_name)
    #print(paste0("Tissue metadata ",tissue_name," : ", as.character(nrow(metadata_tissue))))
    #get a subset of genes and bidirs by chromosome id
    tpms_chrm <- subset(tpms_datatable, chrom == chromosome_id)
    
    #filter samples that match the tissue of interest
    #tpms_chrms_tissue <- t(tpms_chrm[ ,colnames(tpms_chrm) %in% metadata_tissue$sample_name, with=FALSE]) 
    #colnames(tpms_chrms_tissue) <- tpms_chrm$gene_transcript
    
    # transpose the counts table 
    tpms_chrm_t <- t(tpms_chrm[, colnames(tpms_chrm) %in% sample_ids, with = FALSE])
    colnames(tpms_chrm_t) <- tpms_chrm$gene_transcript

    # log transform the matrix of tpms
    tpms_chrms_log10 <- log(tpms_chrm_t+1, base=10)

    # make sure that samples with 0 counts are excluded from the log() transformation
    #log transform the normalized tpm counts(base 10)
    ##NOTE: running with adding 1s and converting 0s to NA
     
    tpms_chrms_log10_NAs <- tpms_chrms_log10

    tpms_chrms_log10_NAs[tpms_chrms_log10_NAs == 0] <- NA
    #tpms_chrms_tissue_log10_NAs <- dplyr::na_if(tpms_chrms_tissue_log10, 0)  

    ########################################################
    ##Using WGCNA calculate correlations and p-values for ##
    ##relavant samples with trancriptio                   ##
    ########################################################
    ##calculated pcc, pvalue, z stat, t stat and number of observations (i.e.)
    ##samples with both genes and bidirectionals transcribed         
    corAndPvalueOut <- WGCNA::corAndPvalue(tpms_chrms_log10_NAs, 
                                           use="pairwise.complete.obs")
    
    ##restructure all the matrix outputs long formats 
    corAndPvalueOut_matrix_list <- lapply(corAndPvalueOut, restructure_cor_matrix)
    
    ##combine all summary stats all in one tibble
    corAndPvalueOut_all_tibble <- cor_summary_stats(corAndPvalueOut_matrix_list)
    #print(dim(corAndPvalueOut_all_tibble))
    #print(head(corAndPvalueOut_all_tibble))
    
    ##remove redundant pairs and add metadata
    ##bidir and gene pair ids, distances 
    corAndPvalueOut_pairs <- cor_pair_metadata(tpms_chrm, 
                                               corAndPvalueOut_all_tibble)
    
    ##add tissue summary of counts
    #corAndPvalueOut_pairs$tissue <- tissue_name
    corAndPvalueOut_pairs$percent_transcribed_both <- (corAndPvalueOut_pairs$nObs/length(sample_ids))*100
    
    return(corAndPvalueOut_pairs)

    }


##-----------------------------------------
## Process by gene
##-----------------------------------------

process_by_gene <- function(gene_id, gene_bidir_tpm_df, metadata_celltype, tissue, window){
    
    #' Run correlations pipeline per gene
    #' 
    #' @description This function will run the correlations between genes and
    #' bidirections  
    #' 
    #' @param gene_id : gene id as characters
    #'
    #' @param gene_bidir_tpm_df : normalized count data.table
    #'    
    #' @param metadata_celltype :  metadata data.table 
    #'
    #' @param tissue : tissue id to process
    #'
    #' @param window : window around TSS for gene
    #'
    #' @usage process_by_gene(gene_id, gene_bidir_tpm_df, metadata_celltype, tissue)
    #' @return A data.frame with all pairwise correlations and significance
    #' @export
    
    # get gene to process
    # using the tryCatch for exceptions
    tryCatch(expr = {gene_a <- gene_bidir_tpm_df[grepl(gsub("([()])","\\\\\\1", gene_id),
                                                       gene_bidir_tpm_df$gene_transcript),]

    
    get_transcripts_a <- get_transcripts_in_window(gene_name = gene_a$gene_transcript, 
                                          gene_tpms_df = gene_bidir_tpm_df, 
                                          window = window)
    if (metadata_celltype == FALSE){
    gene_a_pcc <- transcript_pearsons_by_chromosome(tpms_datatable = get_transcripts_a,
                                                    chromosome_id = gene_a$chrom)
    pairs_pcc <- gene_a_pcc
	}
	else {
    gene_a_pcc <- transcript_pearsons_by_chromosome_tissue(tpms_datatable = get_transcripts_a,
                                                           metadata = metadata_celltype, 
                                                           chromosome_id = gene_a$chrom, 
                                                           tissue_name=tissue)
    pairs_pcc <- gene_a_pcc
}
    return(pairs_pcc)

                     },
             error = function(e){
                 message('Gene pairs in this tissue do not exist!')
                 print(e)
             },
             finally = {
                 message('All done, quitting.')
             }
            )
    
}

##-----------------------------------------
## Save files 
##-----------------------------------------

save_pairs <- function(pairs_dt, 
                       output_folder,
                       tissue_name,
                       chromosome_id){
    
    #saving the final dataframes
    #including chromosome ids to the file names
    tissue_name_nospace <- gsub(" ", "-", tissue_name)

    final_path <- paste0(output_folder,
                         'pearson_correlation_',
                         chromosome_id,
                         '_',
                         tissue_name_nospace,
                         '.tsv.gz' )

    data.table::fwrite(pairs_dt,
                       final_path,
                       sep='\t')
    
}

##########################################
##Run analysis                          ##
##########################################

##----------------------------------------
##0 : start run Log
##----------------------------------------

date_time <- format(Sys.time(), "%Y_%B_%d_%H_%M_%S")
tissue_name_nospace <- gsub(" ", "-", tissue)
sink(paste0(output_folder,"log_pearson_correlations_",chromosome_id,"_",tissue_name_nospace,"_",date_time,".txt"))
cat("Running pearson's correlation analysis")
cat("\n")

cat(paste0("Date: ", Sys.Date()))
cat("\n")

print("START")
cat("\n")

##----------------------------------------
##1: get chromosome gene
##----------------------------------------

#get chromosome transcripts
tpms_datatable_chr <- subset(tpms_datatable, chrom == chromosome_id)

#get list of genes on chromosome
gene_name <- tpms_datatable_chr[grepl("_", tpms_datatable_chr$gene_transcript),]$gene_transcript

##----------------------------------------
##2: process and run correlations on genes 
##----------------------------------------

pearson_pairs <- mclapply(gene_name,
                       process_by_gene, 
                       gene_bidir_tpm_df=tpms_datatable_chr,
                       metadata_celltype=metadata, 
                       tissue=tissue,
                           window = window,
                          mc.cores = ncores)

##----------------------------------------
##3: merge and filter pairs for all genes
##----------------------------------------

#merge the list of gene pairs
pearson_pairs_dt <- data.table::rbindlist(pearson_pairs, use.names=TRUE, fill=TRUE)

#filter to keep thos within specified window
pearson_pairs_filt_dt <- subset(pearson_pairs_dt, nObs >= nlimit & abs(distance_tss) <= window)

#recalculate the adjusted p-values
pearson_pairs_filt_dt$adj_p_BH <- p.adjust(pearson_pairs_filt_dt$pval, method = 'BH')

##----------------------------------------
##4: save the pairs as a table
##----------------------------------------

save_pairs(pairs_dt = pearson_pairs_filt_dt,
           output_folder = output_folder,
           tissue_name = tissue,
           chromosome_id = chromosome_id)

cat("Session Summary")
cat("\n")

print(sessionInfo())
cat("\n")

cat("DONE!")
cat("\n")

sink()
