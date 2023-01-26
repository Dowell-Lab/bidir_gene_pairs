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

print("START")
#initalized paths and options for script
tpms_datatable <- data.table::fread(opt$tpms)
chroms <- data.table::fread(opt$chr)
nsamples <- ncol(tpms_datatable)
output_folder <- opt$out

print(paste("Normalized Counts :", as.character(opt$tpms)))
print(paste("Chromosome File   :", as.character(opt$chr)))
print(paste("Output Directory  :", as.character(opt$out)))
print(paste("Number of Samples :", as.character(nsamples)))

##############################################
##matrix correlation by chromosome function ##
##############################################

transcript_pearsons_by_chromosome <- function(chromosome){
    
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
    tpms_chrm <- as.data.frame(subset(tpms_datatable, Chr == chromosome))
    
    #transpose the dataframe (ie. rows are samples and columns are transcripts)
    #this will also convert dataframe to matrix
    tpms_chrm_t <- t(tpms_chrm[c(7:nsamples)])
    colnames(tpms_chrm_t) <- tpms_chrm$Geneid
    
    #log transform the normalized tpm counts(base 10)
    ##NOTE: running with adding 1s and converting 0s to NA
    #tpms_chrm_log10 <- log(tpms_chrm_t, base=10)
    tpms_chrm_log10 <- log(tpms_chrm_t+1, base=10)

    ##NOTE: to run correlations on all pairs ignoring the 
    ##just samples with missing data, convert -Inf from 
    ##the above log to NAs
    #tpms_chrm_log10_NAs <- dplyr::na_if(tpms_chrm_log10, -Inf)
    tpms_chrm_log10_NAs <- dplyr::na_if(tpms_chrm_log10, 0)

    ####################################################
    ##R values                                        ##
    ####################################################
    #calculate the R value from pearsons correlation
    #also constrict the that the samples used have all transcripts transcribed 
    ## To be revisited, as there could be cell type specific correlations ...
    pearsons_R <- WGCNA::cor(tpms_chrm_log10_NAs, use = "pairwise.complete.obs")
    
    #convert the correlation matrix to dataframe
    #this makes it easy to restructure the output
    pearsons_R_df <- as.data.frame(pearsons_R)
    pearsons_R_df$transcript_1 <- rownames(pearsons_R_df) #add rownames (same as colnames)
    
    #change the structure of the dataframe to 3 column dataframe
    #> transcript_1, transcript_2, coefficient
    pearsons_R_long <- pearsons_R_df %>% tidyr::pivot_longer(!transcript_1, 
                                                             names_to = "transcript_2", 
                                                             values_to = "coefficient",
							     values_drop_na = TRUE)
    
    #remove same pair correlations
    pearsons_R_long_unique <- subset(pearsons_R_long, 
                                     transcript_1 != transcript_2)
    
    #remove NA comparisons (these are transcripts not in all samples)
    pearsons_R_long_unique_noNA <- subset(pearsons_R_long_unique, 
                                          !is.na(pearsons_R_long_unique[[3]]))
    
    pearsons_R_long_unique_noNA$pair_id <- paste0(pearsons_R_long_unique_noNA$transcript_1,
                                              ':',
                                              pearsons_R_long_unique_noNA$transcript_2)
    
    ####################################################
    ##P-values                                        ##
    ####################################################
    #Calculate p-values
    pearsons_pval <- WGCNA::corPvalueStudent(pearsons_R, nsamples)
    
    #converting matrix to dataframe for easy wrangling
    pearsons_pval_df <- as.data.frame(pearsons_pval)
    pearsons_pval_df$transcript_1 <- rownames(pearsons_pval_df)

    #pivot the dataframe structure
    pearsons_pval_long <- pearsons_pval_df %>% tidyr::pivot_longer(!transcript_1, 
                                 names_to = "transcript_2", 
                                 values_to = "p_value",
				 values_drop_na = TRUE)

    #remove same pair correlations
    pearsons_pval_long_unique <- subset(pearsons_pval_long, 
                                        transcript_1 != transcript_2)

    #remove NA comparisons (these are transcripts not in all samples)
    pearsons_pval_long_unique_noNA <- subset(pearsons_pval_long_unique, 
                                          !is.na(pearsons_pval_long_unique[[3]]))

    pearsons_pval_long_unique_noNA$pair_id <- paste0(pearsons_pval_long_unique_noNA$transcript_1,
                                              ':',
                                              pearsons_pval_long_unique_noNA$transcript_2)
    
    #####################################################
    #Final dataframe with both R values and p-values   ##
    #####################################################
    pearsons_Rpval_long_unique_noNA <- dplyr::left_join(pearsons_R_long_unique_noNA, 
                                                    pearsons_pval_long_unique_noNA[c(3,4)],
                                                    by = "pair_id")

    pearsons_Rpval_long_unique_noNA$adj_p_BH <- p.adjust(pearsons_Rpval_long_unique_noNA$p_value,
                                                      method = 'BH')

    #add distance calculations
    
    transcript_coords <- tpms_chrm[c(1,3,4,6)]

    # subset transcript 1 coordinate data
    transcript1_choords <- transcript_coords[transcript_coords$Geneid %in% pearsons_Rpval_long_unique_noNA$transcript_1,]
    colnames(transcript1_choords)[2:4] <- c("transcript1_start","transcript1_end","transcript1_biotype")

    # subset transcript 2 coordinate data
    transcript2_choords <- transcript_coords[transcript_coords$Geneid %in% pearsons_Rpval_long_unique_noNA$transcript_2,]
    colnames(transcript2_choords)[2:4] <- c("transcript2_start","transcript2_end","transcript2_biotype")

    # combine transcript coordinates with correlations
    pearsons_Rpval_txpt1_coord <- dplyr::left_join(pearsons_Rpval_long_unique_noNA,
                                                    transcript1_choords,
                                                    by = c("transcript_1"="Geneid"))

    pearsons_Rpval_txpt_coord <- dplyr::left_join(pearsons_Rpval_txpt1_coord,
                                                    transcript2_choords,
                                                    by = c("transcript_2"="Geneid"))

    pearsons_Rpval_txpt_coord$distance <- pearsons_Rpval_txpt_coord$transcript2_start - pearsons_Rpval_txpt_coord$transcript1_start

    # remove redundant pairs and only report unique Gene-Bidirectional correlations
    pearsons_Rpval_txpt_coord_pairs <- subset(pearsons_Rpval_txpt_coord, 
						transcript1_biotype == 'Gene' & 
    			      	 		transcript2_biotype != 'Gene')   
 
    #saving the final dataframes
    #including date and chromosome ids to the file names
    outdir <- output_folder 
    
    final_path <- paste0(outdir,
                         'pearson_correlation_', 
                         chromosome,
                         '_',
                         Sys.Date(),
                         '.txt' )
    
    write.table(pearsons_Rpval_txpt_coord_pairs, #pearsons_Rpval_long_unique_noNA, 
                final_path,
                sep='\t',
                quote = FALSE,
                row.names=FALSE)

    #system("gzip final_path")
    system(sprintf("gzip %s", final_path))
    
}

chromosome_list <- as.character(chroms$V1) 

parallel::mclapply(chromosome_list, transcript_pearsons_by_chromosome, mc.cores = length(chromosome_list))

print("Session Summary")

print(sessionInfo())

print("DONE!")
