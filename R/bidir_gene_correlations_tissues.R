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
metadata <- read.csv(opt$samplemeta)
chroms <- data.table::fread(opt$chr)
nsamples <- ncol(tpms_datatable)
output_folder <- opt$out

print(paste("Normalized Counts :", as.character(opt$tpms)))
print(paste("Sample metadata   :", as.character(opt$samplemeta)))
print(paste("Chromosome File   :", as.character(opt$chr)))
print(paste("Output Directory  :", as.character(opt$out)))
print(paste("Number of Samples :", as.character(nsamples)))

##############################################
##matrix correlation by chromosome function ##
##############################################

##get metadata for the samples analyzed
sample_ids <- colnames(tpms_datatable[,6:885])
metadata_analyzed <- metadata[metadata$srz %in% sample_ids,]
print(dim(metadata_analyzed))

transcript_pearsons_by_chromosome <- function(chromosome_id, tissue_name){
    
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

    ##load data
    #tpms_datatable <- data.table::fread(opt$tpms)
    #metadata <- read.csv(opt$samplemeta)
    
    ##get metadata for the samples analyzed
    sample_ids <- colnames(tpms_datatable[,6:885])
    metadata_analyzed <- metadata[metadata$srz %in% sample_ids,]

    # get metadata for specific tissue of interest
    metadata_tissue <- subset(metadata_analyzed, tissue == tissue_name)
    print(paste("Tissue metadata : ", as.character(nrow(metadata_tissue))))
   
    #get a subset of genes and bidirs by chromosome id
    tpms_chrm <- subset(tpms_datatable, chromosome == chromosome_id)
    
    #filter samples that match the tissue of interest
    tpms_chrms_tissue <- t(tpms_chrm[ ,colnames(tpms_chrm) %in% metadata_tissue$srz, with=FALSE]) 
    colnames(tpms_chrms_tissue) <- tpms_chrm$GeneID

    # log transform the matrix of tpms
    tpms_chrms_tissue_log10 <- log(tpms_chrms_tissue+1, base=10)

    # make sure that samples with 0 counts are excluded from the log() transformation
    #log transform the normalized tpm counts(base 10)
    ##NOTE: running with adding 1s and converting 0s to NA
    tpms_chrms_tissue_log10_NAs <- dplyr::na_if(tpms_chrms_tissue_log10, 0)


    #########################
    ##Pearsons correlation ##
    #########################
    # Using WGCAs cor function
    pearsons_R_chrms_tissue <- WGCNA::cor(tpms_chrms_tissue_log10_NAs, 
                                     use = "pairwise.complete.obs")

    #convert the correlation matrix to dataframe
    #this makes it easy to restructure the output
    pearsons_R_df_chrms_tissue <- as.data.frame(pearsons_R_chrms_tissue)
    pearsons_R_df_chrms_tissue$transcript_1 <- rownames(pearsons_R_df_chrms_tissue) #add rownames (same as colnames)

    #change the structure of the dataframe to 3 column dataframe
    #> transcript_1, transcript_2, coefficient
    pearsons_R_long_chrms_tissue <- pearsons_R_df_chrms_tissue %>% tidyr::pivot_longer(!transcript_1, 
                                                             names_to = "transcript_2", 
                                                             values_to = "coefficient",
                                 values_drop_na = TRUE)

    #remove same pair correlations
    pearsons_R_long_unique_chrms_tissue <- subset(pearsons_R_long_chrms_tissue, 
                                             transcript_1 != transcript_2)

    #remove NA comparisons (these are transcripts not in all samples)
    pearsons_R_long_unique_noNA_chrms_tissue <- subset(pearsons_R_long_unique_chrms_tissue, 
                                          !is.na(pearsons_R_long_unique_chrms_tissue[[3]]))

    pearsons_R_long_unique_noNA_chrms_tissue$pair_id <- paste0(pearsons_R_long_unique_noNA_chrms_tissue$transcript_1,
                                              ':',
                                              pearsons_R_long_unique_noNA_chrms_tissue$transcript_2)

    #########################
    ## P-value calculation ##
    #########################
    pearsons_pval_chrms_tissue <- WGCNA::corPvalueStudent(pearsons_R_chrms_tissue,
                                                 nrow(tpms_chrms_tissue))

    #converting matrix to dataframe for easy wrangling
    pearsons_pval_df_chrms_tissue <- as.data.frame(pearsons_pval_chrms_tissue)

    #converting matrix to dataframe for easy wrangling
    pearsons_pval_df_chrms_tissue <- as.data.frame(pearsons_pval_chrms_tissue)
    pearsons_pval_df_chrms_tissue$transcript_1 <- rownames(pearsons_pval_df_chrms_tissue)

    #pivot the dataframe structure
    pearsons_pval_long_chrms_tissue <- pearsons_pval_df_chrms_tissue %>% tidyr::pivot_longer(!transcript_1, 
                                 names_to = "transcript_2", 
                                 values_to = "p_value",
                 values_drop_na = TRUE)

    #remove same pair correlations
    pearsons_pval_long_unique_chrms_tissue <- subset(pearsons_pval_long_chrms_tissue, 
                                        transcript_1 != transcript_2)

    #remove NA comparisons (these are transcripts not in all samples)
    pearsons_pval_long_unique_noNA_chrms_tissue <- subset(pearsons_pval_long_unique_chrms_tissue, 
                                          !is.na(pearsons_pval_long_unique_chrms_tissue[[3]]))

    pearsons_pval_long_unique_noNA_chrms_tissue$pair_id <- paste0(pearsons_pval_long_unique_noNA_chrms_tissue$transcript_1,
                                              ':',
                                              pearsons_pval_long_unique_noNA_chrms_tissue$transcript_2)
    
    #####################################################
    #Final dataframe with both R values and p-values   ##
    #####################################################
    pearsons_Rpval_long_unique_noNA_chrms_tissue <- dplyr::left_join(pearsons_R_long_unique_noNA_chrms_tissue, 
                                                    pearsons_pval_long_unique_noNA_chrms_tissue[c(3,4)],
                                                    by = "pair_id")

    pearsons_Rpval_long_unique_noNA_chrms_tissue$adj_p_BH <- p.adjust(pearsons_Rpval_long_unique_noNA_chrms_tissue$p_value,
                                                      method = 'BH')
    
    #add distance calculations
    transcript_coords_chrms_tissue <- tpms_chrm[ ,colnames(tpms_chrm) %in% c('GeneID','chromosome','start','end','bidir_id'), with=FALSE]
    
    # subset transcript 1 coordinate data
    transcript1_choords_chrms_tissue <- transcript_coords_chrms_tissue[transcript_coords_chrms_tissue$GeneID %in% pearsons_Rpval_long_unique_noNA_chrms_tissue$transcript_1,]
    colnames(transcript1_choords_chrms_tissue) <- c("chromosome","transcript1_start","transcript1_end","GeneID","transcript1_biotype")
    
    # subset transcript 2 coordinate data
    transcript2_choords_chrms_tissue <- transcript_coords_chrms_tissue[transcript_coords_chrms_tissue$GeneID %in% pearsons_Rpval_long_unique_noNA_chrms_tissue$transcript_2,]
    colnames(transcript2_choords_chrms_tissue) <- c("chromosome","transcript2_start","transcript2_end","GeneID","transcript2_biotype")
    
    # combine transcript coordinates with correlations
    pearsons_Rpval_txpt1_coord_chrms_tissue <- dplyr::left_join(pearsons_Rpval_long_unique_noNA_chrms_tissue,
                                                    transcript1_choords_chrms_tissue,
                                                    by = c("transcript_1"="GeneID"))

    pearsons_Rpval_txpt_coord_chrms_tissue <- dplyr::left_join(pearsons_Rpval_txpt1_coord_chrms_tissue,
                                                    transcript2_choords_chrms_tissue,
                                                    by = c("transcript_2"="GeneID"))

    pearsons_Rpval_txpt_coord_chrms_tissue$distance <- pearsons_Rpval_txpt_coord_chrms_tissue$transcript2_start - pearsons_Rpval_txpt_coord_chrms_tissue$transcript1_start
    
    # remove redundant pairs and only report unique Gene-Bidirectional correlations
    pearsons_Rpval_txpt_coord_pairs_chrms_tissue <- subset(pearsons_Rpval_txpt_coord_chrms_tissue, 
                        transcript1_biotype == 'Gene' & 
                                transcript2_biotype != 'Gene') 
 
    #saving the final dataframes
    #including date and chromosome ids to the file names
    outdir <- output_folder 

    tissue_name_nospace <- gsub(" ", "-", tissue_name)
    
    final_path <- paste0(outdir,
                         'pearson_correlation_', 
                         chromosome_id,
                         '_', 
			 tissue_name_nospace,
			 '_',
                         Sys.Date(),
                         'noMissingSamples.txt' )
    
    write.table(pearsons_Rpval_txpt_coord_pairs_chrms_tissue, #pearsons_Rpval_long_unique_noNA, 
                final_path,
                sep='\t',
                quote = FALSE,
                row.names=FALSE)

    #system("gzip final_path")
    system(sprintf("gzip %s", final_path))
    
}

chromosome_list <- as.character(chroms$V1) 
print(paste("Processing # chromosomes : ",as.character(length(chromosome_list)) ))
#parallel::mclapply(chromosome_list, transcript_pearsons_by_chromosome, mc.cores = length(chromosome_list))

#get samples used in analyses
tissue_counts <- as.data.frame(table(as.character(metadata_analyzed$tissue)))
tissue_list <- as.character(subset(tissue_counts, Freq>15)$Var1)
print(paste("Processing # tissues : ",as.character(length(tissue_list))))

cores_requested <- makeCluster(1, type="FORK")
parallel::clusterMap(cores_requested, transcript_pearsons_by_chromosome, chromosome_list, tissue_list)
stopCluster(cores_requested) 

print("Session Summary")

print(sessionInfo())

print("DONE!")
