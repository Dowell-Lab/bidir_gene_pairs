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
        make_option(c("-m", "--samplemeta"), type="character", default=NULL,
              help="path to metadata table for all samples", metavar="character"),
        make_option(c("-r", "--correlations"), type="character", default=NULL,
              help="path to chromosome level correnations from nascent_transcripts_pearsons_corr.R", metavar="character"),
        make_option(c("-c", "--chr"), type="character", default="NULL",     
              help="chromosome file for processing", metavar="character"),
	make_option(c("-d", "--min_dist"), type="integer", default=1000000, 
	      help="minimum distance between pairs [default = %default bases]", metavar="integer"),
	make_option(c("-e", "--percent_trans"), type="integer", default=25, 
	      help="percent of sample with transcription for transcript [default = %default]", metavar="integer"),
	make_option(c("-v", "--rvalue"), type="double", default=0.9, 
	      help="pearson's R value cut-off [default = %default ]", metavar="double"),
	make_option(c("-p", "--adj_pvalue"), type="double", default=0.001, 
	      help="adjusted p-value filter for called pairs [default = %default ]", metavar="double"),
	make_option(c("-o", "--out"), type="character", default="./", 
              help="path to output directory [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$tpms)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input tpms).n", call.=FALSE)
}

print("START")
#initalized paths and options for script
tpms_datatable <- data.table::fread(opt$tpms)
metadata <- read.csv(opt$samplemeta)
chroms <- data.table::fread(opt$chr)
nsamples <- ncol(tpms_datatable)
output_folder <- opt$out
correlation_files <- opt$correlations

#initialize parameter variables inputed via arg parse
minimum_distance <- opt$min_dist
percent_transcribed <- opt$percent_trans
rvalue_cutoff <- opt$rvalue
adjusted_pvalue <- opt$adj_pvalue 

print(paste("Normalized Counts :", as.character(opt$tpms)))
print(paste("Sample metadata   :", as.character(opt$samplemeta)))
print(paste("Chromosome File   :", as.character(opt$chr)))
print(paste("Output Directory  :", as.character(opt$out)))
print(paste("Number of Samples :", as.character(nsamples)))


####################################################################################
##Given an input correlation file, filter significant pairs based on distance etc ##
####################################################################################

get_sig_pairs <- function(erna_gene_pairs_tissue, metadata, tpms, 
                            min_distance = minimum_distance, r_value = rvalue_cutoff,
                            adjusted_p = adjusted_pvalue, perc_trans = percent_transcribed){
    
    #' Get bidirectionals associated with a gene by tissue correlations
    #' 
    #' @description This function will calculate % of samples supporting the 
    #' correlations by tissue and return pairs sipported by majority of samples 
    #' 
    #' @param metadata data.table with meta data for samples for the specific
    #' tissue
    #'
    #'
    #' @param tpms data.table with normalized counts for all samples
    #'
    #' @param tissue_id name of tissue to filter on
    #'
    #' @usage transcript_pearsons_by_chromosome(tpms_datatable, chromosome, output_folder)
    #' @return A data.frame with all pairwise correlations and significance
    #' @export
    
    #get tissue id 
    tissue_id <- unique(erna_gene_pairs_tissue$tissue)
    
    #filter sample meta data by tissue of interest
    tissue <- subset(metadata, tissue == tissue_id)

    #get the number of samples in tissues 
    num_tissue <- nrow(tissue)

    #filter the tpms for the samples within specific tissue id
    tpms_tissue <- tpms[ ,
                        colnames(tpms) %in%
                        tissue$srz, with=FALSE]

    #add a gene id amd calculate the percent of samples with transcription
    tpms_tissue$GeneID <- tpms$GeneID
    tpms_tissue$num_zeros <- rowSums(tpms_tissue[,1:num_tissue] == 0)
    tpms_tissue$num_transcribed <- num_tissue-tpms_tissue$num_zeros
    tpms_tissue$percent_transcribed <- ((num_tissue-tpms_tissue$num_zeros)/num_tissue)*100
  
    #siginificant pairs
    gene_bidirs_sig <- subset(erna_gene_pairs_tissue, 
                         abs(distance) < min_distance &
                         abs(coefficient) > r_value &
                         adj_p_BH < adjusted_p ) 
    
    #merge with sample counts meta data
    #1: get the % transcription rate for BIDIRECTIONALS
    gene_bidirs_meta_txp2 <- merge(gene_bidirs_sig,
                              tpms_tissue[,c("GeneID",
                                             "percent_transcribed")],
                              by.x="transcript_2",
                              by.y="GeneID")
    
    names(gene_bidirs_meta_txp2)[names(gene_bidirs_meta_txp2) == 'percent_transcribed'] <- 'transcript2_percent_transcribed'

    #2: get the % transcription rate for GENES
    gene_bidirs_meta <- merge(gene_bidirs_meta_txp2,
                              tpms_tissue[,c("GeneID",
                                             "percent_transcribed")],
                              by.x="transcript_1",
                              by.y="GeneID")
    
    names(gene_bidirs_meta)[names(gene_bidirs_meta) == 'percent_transcribed'] <- 'transcript1_percent_transcribed'

    gene_bidirs_meta$confidence <- (gene_bidirs_meta$coefficient*((gene_bidirs_meta$transcript1_percent_transcribed+gene_bidirs_meta$transcript2_percent_transcribed)/2))/100

    #transcribed by most samples
    gene_bidirs_transcribed <- subset(gene_bidirs_meta,
                                      transcript1_percent_transcribed > perc_trans &
                                      transcript2_percent_transcribed > perc_trans)

    return(gene_bidirs_transcribed)
    
}

##########################################
###Processing pipeline by chromosome   ###
##########################################

processing_filters_by_chr <- function(chromosome_id, corr_path, 
                                      metadata_counted, norm_counts){
    
    ##get paths for the counts tables
    corr_files <- list.files(path=corr_path, 
                              pattern=paste0("^pearson_correlation_",chromosome_id,"_"),
                             full.names=TRUE)
    
    ##load the data
    corr_DT_list <- lapply(corr_files, 
                           data.table::fread) 

    ##add tissue names
    file_names <- as.character(tools::file_path_sans_ext(basename(corr_files)))
    tissue_names <- gsub("-", " ", as.character(lapply(strsplit(file_names, '_'), `[`, 4)))

    for (i in 1:length(corr_DT_list)){

        corr_DT_list[[i]]$tissue <- tissue_names[i]

    }

    ##filter significant pairs
    sig_corr_DT_list <- lapply(corr_DT_list, 
                               get_sig_pairs, 
                               metadata=metadata_counted,
                               tpms=norm_counts)
    
    #combine all chromosome pairs
    sig_corr_DT <- do.call(rbind, sig_corr_DT_list)

    #saving the final dataframes
    #including date and chromosome ids to the file names
    final_path <- paste0(output_folder,
                         chromosome_id,
                         '_significant',
			 rvalue_cutoff,
			 '_pairs_by_tissue_',
                         Sys.Date(),
                         '.txt' )

    write.table(sig_corr_DT,
                final_path,
                sep='\t',
                quote = FALSE,
                row.names=FALSE)

    #compress file		
    system(sprintf("gzip %s", final_path))

}


##get metadata for the samples analyzed
sample_ids <- colnames(tpms_datatable[,6:885])
metadata_analyzed <- metadata[metadata$srz %in% sample_ids,]

chromosome_list <- as.character(chroms$V1) 
print(paste("Processing # chromosomes : ",as.character(length(chromosome_list)) ))
parallel::mclapply(chromosome_list, 
	           processing_filters_by_chr, 
		   corr_path = correlation_files,
		   metadata_counted = metadata_analyzed,
		   norm_counts = tpms_datatable,
		   mc.cores = length(chromosome_list))

print("Session Summary")

print(sessionInfo())

print("DONE!")
