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
        make_option(c("-r", "--correlations"), type="character", default=NULL,
              help="path to chromosome level correnations from bidir_gene_correlations_allsamples.R or  bidir_gene_correlations_tissues.R", metavar="character"),
        make_option(c("-c", "--chr"), type="character", default="NULL",     
              help="chromosome file for processing", metavar="character"),
	make_option(c("-d", "--min_dist"), type="integer", default=1000000, 
	      help="minimum distance between pairs [default = %default bases]", metavar="integer"),
	make_option(c("-e", "--percent_trans"), type="integer", default=5, 
	      help="percent of sample with transcription for transcript [default = %default]", metavar="integer"),
	make_option(c("-v", "--rvalue"), type="double", default=0.6, 
	      help="pearson's R value cut-off [default = %default ]", metavar="double"),
	make_option(c("-p", "--adj_pvalue"), type="double", default=0.01, 
	      help="adjusted p-value filter for called pairs [default = %default ]", metavar="double"),
	make_option(c("-o", "--out"), type="character", default="./", 
              help="path to output directory [default= %default]", metavar="character"),
        make_option(c("-t", "--tissue"), action="store_true", default=FALSE,
              help="Are the correlations based on a per-tissue and per-chromosome basis? [default = FALSE]")	      
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$correlations)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input correlations).n", call.=FALSE)
}

#initalized paths and options for script
chroms <- data.table::fread(opt$chr)
output_folder <- opt$out
correlation_files <- opt$correlations

#initialize parameter variables inputed via arg parse
minimum_distance <- opt$min_dist
percent_transcribed <- opt$percent_trans
rvalue_cutoff <- opt$rvalue
adjusted_pvalue <- opt$adj_pvalue 

#add log file
date_time <- format(Sys.time(), "%Y_%B_%d_%H_%M_%S")
sink(paste0(output_folder,"log_filter_",date_time,".txt"))
cat("Filtering correlated bidirectional gene pairs")
cat("\n")

cat(paste0("Date: ", Sys.Date()))
cat("\n")

print("START")
cat("\n")

print(paste("Minimum distance (bp)            :", as.character(opt$min_dist)))
cat("\n")

print(paste("|PCC| >                          :", as.character(opt$rvalue)))
cat("\n")

print(paste("Adjusted p-value <               :", as.character(opt$adj_pvalue)))
cat("\n")

print(paste("Minimun % samples in correlation :", as.character(opt$percent_trans)))
cat("\n")
####################################################################################
##Given an input correlation file, filter significant pairs based on distance etc ##
####################################################################################

get_sig_pairs <- function(gene_bidir_pairs, min_distance = minimum_distance, r_value = rvalue_cutoff,
                            adjusted_p = adjusted_pvalue, perc_trans = percent_transcribed){
    
    #' Get bidirectionals associated with a gene by tissue correlations
    #' 
    #' @description This function will filter % of samples supporting the 
    #' correlations based on number of observations and return pairs supported by % of samples 
    #' 
    #' @param gene_bidir_pairs data.table with gene and bidirectional pair correlations
    #'
    #' @param min_distance
    #'
    #' @param r_value
    #'
    #' @param adjusted_p
    #'
    #' @param perc_trans
    #'
    #' @usage get_sig_pairs(gene_bidir_pairs)
    #' @return A data.frame with significant pairwise correlations and significance
    #' @export
    
    #siginificant pairs
    gene_bidirs_sig <- subset(gene_bidir_pairs, 
                         abs(distance_tss) < min_distance &
                         abs(pcc) > r_value &
                         adj_p_BH < adjusted_p &
			 percent_transcribed_both > perc_trans) 
    
    return(gene_bidirs_sig)
    
}

##########################################
###Processing pipeline by chromosome   ###
##########################################

processing_filters_by_chrNtissue <- function(chromosome_id, corr_path){
    
    ##get paths for the counts tables
    corr_files <- list.files(path=corr_path, 
                              pattern=paste0("^pearson_correlation_",chromosome_id,"_"),
                             full.names=TRUE)
    
    ##load the data
    corr_DT_list <- lapply(corr_files, 
                           data.table::fread) 


    ##filter significant pairs
    sig_corr_DT_list <- lapply(corr_DT_list, 
                               get_sig_pairs)
    
    #combine all chromosome pairs
    sig_corr_DT <- do.call(rbind, sig_corr_DT_list)

    #saving the final dataframes
    #including date and chromosome ids to the file names
    final_path <- paste0(output_folder,
                         chromosome_id,
                         '_significant_pairs_by_tissue.tsv.gz')

    data.table::fwrite(sig_corr_DT,
			final_path,
			sep='\t')

}


processing_filters_by_chr <- function(chromosome_id, corr_path){
    
    ##get paths for the counts tables
    corr_files <- list.files(path=corr_path, 
                              pattern=paste0("^pearson_correlation_",chromosome_id),
                             full.names=TRUE)
    
    #load file
    corr_DT <- data.table::fread(corr_files[1])	

    #filter pairs
    sig_corr_DT <- get_sig_pairs(corr_DT)	
    
    #saving the final dataframes
    #including chromosome ids to the file names
    final_path <- paste0(output_folder,
                         chromosome_id,
                         '_significant_pairs.tsv.gz' )

    data.table::fwrite(sig_corr_DT,
			final_path,
			sep='\t')

}

chromosome_list <- as.character(chroms$V1) 
print(paste("Processing # chromosomes : ",as.character(length(chromosome_list)) ))
cat("\n")

# print some progress messages to stderr if "quietly" wasn't requested
if ( opt$tissue ) { 
    print("Getting significant pairs from tissue correlations") 
    parallel::mclapply(chromosome_list, 
	           processing_filters_by_chrNtissue, 
		   corr_path = correlation_files,
		   mc.cores = length(chromosome_list))

} else {
    print("Getting significant pairs from all sample correlations")
    parallel::mclapply(chromosome_list, 
	           processing_filters_by_chr, 
		   corr_path = correlation_files,
		   mc.cores = length(chromosome_list))


}

cat("\n")
print("Session Summary")
cat("\n")

print(sessionInfo())
cat("\n")

print("DONE!")
cat("\n")

sink()


