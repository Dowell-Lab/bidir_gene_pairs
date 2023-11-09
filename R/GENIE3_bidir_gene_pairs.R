#! /usr/bin/env Rscript
##################
##load packages ##
##################

suppressMessages(library(GENIE3)) ## random forest feature selection
suppressMessages(library(data.table)) ## load files into R faster
suppressMessages(library(parallel)) ##running code in parallel
suppressMessages(library(optparse)) ##adding arguments 

#############################################
## Initialize command options for script   ##
#############################################
# define input and output options
option_list = list(
        make_option(c("-c", "--correlations"), type="character", default=NULL,
              help="path to significant pairs", metavar="character"),
        make_option(c("-g", "--gtex"), type="character", default="NULL",     
              help="significant pairs supported by gtex", metavar="character"),
        make_option(c("-t", "--normalized_counts"), type="character", default="NULL",     
              help="normalized counts table", metavar="character"),
        make_option(c("-m", "--metadata"), type="character", default="NULL",     
              help="sample metadata used to filter counts for analysis", metavar="character"),
	make_option(c("-d", "--dopas"), type="integer", default=15000, 
	      help="distance of bidirectionals to be removed downstream on the polyadenylation site (DoPAS) [default DoPAS = %default bp]", metavar="integer"),
	make_option(c("-n", "--nobs"), type="integer", default=10, 
	      help="minimum number of observations for a correlation to be considered [default samples = %default]", metavar="integer"),
	make_option(c("-s", "--seed"), type="integer", default=123, 
	      help="default seed for the random forest GENIE3 analysis [default seed = %default]", metavar="integer"),
	make_option(c("-o", "--out"), type="character", default="./", 
              help="path to output directory [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$correlations)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input correlations).n", call.=FALSE)
}

#initalized paths and options for script
output_folder <- opt$out
correlation_files <- opt$correlations
metadata_celltype <- opt$metadata
gtex_pairs <- opt$gtex
counts <- opt$normalized_counts

#initialize parameter variables inputed via arg parse
minimum_distance <- opt$dopas
seed_val <- opt$seed
nObs <- opt$nobs

sink_file_name <- base::basename(correlation_files)

sink(paste0(output_folder,"log_filter_",sink_file_name,"_",Sys.Date(),".txt"))
cat("Filtering correlated bidirectional gene pairs")
cat("\n")

cat(paste0("Date: ", Sys.Date()))
cat("\n")

print("START")
cat("\n")

print(paste("Correlations file                :", as.character(correlation_files)))
cat("\n")

print(paste("Normalized counts                :", as.character(counts)))
cat("\n")

print(paste("Metadata                         :", as.character(metadata_celltype)))
cat("\n")

print(paste("GTEx pairs                       :", as.character(gtex_pairs)))
cat("\n")

print(paste("Output folder                    :", as.character(output_folder)))
cat("\n")

print(paste("Seed                             :", as.character(seed_val)))
cat("\n")

####################################################################################
## Set the seed!                                                                  ##
####################################################################################
set.seed(seed_val)

####################################################################################
##Further filter correlations by number of samples they are found in and DoPAS    ##
####################################################################################
removing_DoPAD_bidirs <- function(pairs, dist_DoPAD=minimum_distance, number_observed=nObs) {
    #' Remove pairs where the bidirectionals transcript is Downstream of Polyadenylation Site (DoPAD)
    #1: Subset by distance from TES
    pairs_dopad <- subset(pairs, abs(distance_tes) < dist_DoPAD)
    
    #2: Label bidirectionals that fall closer to the TES than to the TSS
    pairs_dopad$filter <- ifelse(abs(pairs_dopad$distance_tss) < abs(pairs_dopad$distance_tes),
                                        FALSE, TRUE)
    
    pairs_dopad_removing <- pairs_dopad[pairs_dopad$filter=='TRUE',]

    #This removes any pairs with bidirectionals in DoPAD as opposed to only removing the pair
    pairs_keeping <- pairs[!pairs$transcript_2 %in% pairs_dopad_removing$transcript_2,]
    
    pairs_keeping_nObs <- subset(pairs_keeping, nObs > number_observed)
    return(pairs_keeping_nObs)
}


####################################################################################
##Running GENIE3 for each gene in the database                                    ##
####################################################################################
genie3_feature_ranks <- function(sig_pairs_dt, counts_tpm, metadata, gtex_table, gene_id, tissue_id, seed=seed_val){

    set.seed(seed)
    
    #1: get pairs by tissue for given gene
    gene_sig_pairs <- subset(sig_pairs_dt, 
                             transcript_1 == gene_id & 
                             tissue == tissue_id)

    #2: subset metadata for tissue
    metadata_celltype <- subset(metadata, tissue == tissue_id)
    
    #3: get normalized counts for the tissues
    gene_bidirs_tpms <- counts_tpm[counts_tpm$gene_transcript %in% c(gene_id, gene_sig_pairs$transcript_2),
                                   colnames(counts_tpm) %in%
                                   c("gene_transcript",
                                     metadata_celltype$sample_name), with=FALSE]
    
    #4: convert counts data.table to matrix
    gene_bidirs_matrix <- as.matrix(gene_bidirs_tpms[,2:ncol(gene_bidirs_tpms)])
    rownames(gene_bidirs_matrix) <- gene_bidirs_tpms$gene_transcript

    #5: list all unique bidirectionals
    gene_bidirs <- unique(gene_sig_pairs$transcript_2)
    
    #6: run GENIE3
    gene_weightMat <- GENIE3(gene_bidirs_matrix, regulators=gene_bidirs)
    
    #7: get all links
    gene_linkList <- getLinkList(gene_weightMat)
    
    #8: get final ranks
    gene_linkList_genes <- subset(gene_linkList, targetGene==gene_id)
    gene_linkList_genes$pair_id <- paste0(gene_linkList_genes$targetGene, '~',
                                          gene_linkList_genes$regulatoryGene)
    gene_linkList_genes$gtex <- ifelse(gene_linkList_genes$pair_id %in% gtex_table$pair_id,
                                        1, 0)
    gene_linkList_genes$tissue <- tissue_id
    gene_linkList_genes$ranks <- as.numeric(seq(1,nrow(gene_linkList_genes), 1))
    
    return(gene_linkList_genes)
} 

#######################################################################################
##Run GENIE3 per gene and across tissues                                             ##
#######################################################################################
print("Processing Gene~Bidirectional pairs with GENIE3  ")
cat("\n")

#load counts
counts_dt <- data.table::fread(counts)

#Import gtex pairs
gtex_dt <- data.table::fread(gtex_pairs)

#load significant pairs
sig_pairs_dt <- data.table::fread(correlation_files)
sig_pairs_dt_filt <- removing_DoPAD_bidirs(sig_pairs_dt)

#load metadata
metadata_dt <- data.table::fread(metadata_celltype)

genes <- as.character(unique(sig_pairs_dt_filt$transcript_1))
tissues <- as.character(unique(sig_pairs_dt_filt$tissue))

#create a list for the results data.tables
genie3_out <- list()

for (i in genes) {

     for(j in tissues) {
         
         #get gene from - tissue
         gene_tissue <-  paste(i,'-->',j) 
         print("Gene:Transcript --> Tissue")
         print(gene_tissue)
         tryCatch(expr = {genie3_out[[gene_tissue]] <- genie3_feature_ranks(sig_pairs_dt=sig_pairs_dt_filt, 
	 	       	 			       			    counts_tpm=counts_dt, 
                                                      			    metadata=metadata_dt, 
                                                      			    gtex_table=gtex_dt,				 
                                                      			    gene_id=i, 
                                                      			    tissue_id=j)
                         						    },
                  error = function(e){
                      message('Gene pairs in this tissue do not exist!')
                      print(e)
                  },
                  warning = function(w){
                      message('Pair does not exist??')
                      print(w)
                  },
                  finally = {
                      message('All done, quitting.')
                  }
                 )
    
     }
}

genie3_out_table <- do.call(rbind, genie3_out)

correlation_files_basename <- base::basename(correlation_files)
correlation_files_name <- as.character(lapply(strsplit(correlation_files_basename, '\\.'), `[`, 1))

##new file name
genie3_out_file <- paste0(output_folder, correlation_files_name, "_GENIE3_ranks.tsv.gz")
data.table::fwrite(genie3_out_table,
                  genie3_out_file,
                  sep='\t')

cat("\n")
print("Session Summary")
cat("\n")

print(sessionInfo())
cat("\n")

print("DONE!")
cat("\n")

sink()


