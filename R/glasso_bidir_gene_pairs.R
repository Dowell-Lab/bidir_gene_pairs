#! /usr/bin/env Rscript

##################
##load packages ##
##################

#suppressMessages(library(WGCNA)) ## faster cor()
suppressMessages(library(dplyr)) ## for the R pipes
suppressMessages(library(tidyr)) ## for tidying the dataframes
suppressMessages(library(data.table)) ## load files into R faster
suppressMessages(library(parallel)) ## running code in parallel
suppressMessages(library(optparse)) ## adding arguments 
suppressMessages(library(glasso)) ## running glasso on pairs
suppressMessages(library(reshape2)) ## restructure matrix

#############################################
## Initialize command options for script   ##
#############################################
# define input and output options
option_list = list(
    make_option(c("-t", "--tpms"), type="character", default=NULL,
                help="path to TPM normalized counts", metavar="character"),
    make_option(c("-m", "--samplemeta"), type="character", default=NULL,
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
    make_option(c("-d", "--dist"), type="integer", default=250,
                help="minimum distance in kb to include in penalty [default = %default kb]",
                metavar="integer"),
    make_option(c("-p", "--param_scale"), type="integer", default=2,
                help="scaling parameter controlling for when we expect meaningful contacts [default = %default]",
                metavar="integer"),
    make_option(c("-s", "--constant"), type="double", default=0.75,
                help="constant following power-law distribution [default = %default]",
                metavar="double"),
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
##Initialize variables                   ##
###########################################
#files and paths
tpms_datatable <- data.table::fread(opt$tpms)
metadata <- data.table::fread(opt$samplemeta)
output_folder <- opt$out
chromosome_id <- opt$chr_id

###########################################
##Processing functions                   ##
###########################################
get_transcripts_in_window <- function(gene_name, gene_tpms_df, window){
    
    ##get transcripts including genes within the 1MB window gene
    window <- window
    
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

restructure_matrix <- function(mat_in){
    
    # restructure matrix with melt
    # use the data.table::melt function
    melted_dt <- data.table::melt(as.data.table(mat_in, 
                                              keep.rownames=TRUE), 
                                measure=patterns("[0-9]"))
    
    # add column names
    colnames(melted_dt) <- c("transcript_1","transcript_2", "value")
    melted_dt$pair_id <- paste0(melted_dt$transcript_1,
                                "~",
                                melted_dt$transcript_2)
    
    return(melted_dt)
    
}


get_rho_mat <- function(dist_matrix, distmin=opt$dist, distance_parameter=opt$param_scale, s=opt$constant) {
    xmin <- distmin #kilobases as a minimum to include in the penalty

    out <- (1-(xmin/dist_matrix)^s) * distance_parameter
    out[!is.finite(out)] <- 0
    out[out < 0] <- 0
    return(out)
}

filter_low_transcribed <- function(tpm_matrix, nlimit=opt$nlimit){
    
    #binary for trancribed or not
    tpm_matrix_binary <- ifelse(tpm_matrix != 0, 1, 0)

    #count the number of samples with counts
    num_samples <- as.data.frame(colSums(tpm_matrix_binary))
    colnames(num_samples) <- "N"
    num_samples$transcripts <- rownames(num_samples)
    
    #filter only transcripts with more than nlimit samples
    num_samples_transcribed <- subset(num_samples, N>=nlimit)
    
    #filter samples with greater than or equal to nlimits samples
    tpm_transcribed_matrix <- tpm_matrix[, 
                                         colnames(tpm_matrix) %in%
                                         num_samples_transcribed$transcripts]
    return(tpm_transcribed_matrix)
    
}


get_distance_matrix <- function(selected_transcripts, nsamples = nsamples){
    
    #get number of transcripts in data set
    n_txpt <- nsamples #nrow(selected_transcripts)
    totals <- n_txpt*n_txpt
    
    cat("Total pairs       = ", totals, "\n")
    cat("number of samples = ", n_txpt)
    
    #create a matrix with 0s
    dist_txpt <- matrix(rep(0, totals),
                       nrow=n_txpt,
                       ncol=n_txpt)

    #loop through all the transcripts and calculate distance in kbs
    for (i in 1:n_txpt){
        for (j in 1:n_txpt){
            
            #printing pairs
            #print(paste0(i,"~",j))
            
            txpt_start_i <- selected_transcripts$start[i] 
            txpt_start_j <- selected_transcripts$start[j] 

            dist_txpt[i,j] <- abs(txpt_start_i-txpt_start_j)/1000
           
        }
    
    }
    
    return(dist_txpt)
}


process_glasso <- function(gene_name, 
                           gene_bidir_tpms,
                           metadata_tissue, 
                           window, 
                           exclude_missing_data){
    
    ##-----------------------------------------------------------
    ##get the gene of interest and bidirectionals in window
    ##-----------------------------------------------------------
    gene_sel <- gene_bidir_tpms[grepl(gene_name,
                                      gene_bidir_tpms$gene_transcript),]
    
    gene_bidir_sel <- get_transcripts_in_window(gene_name = gene_sel$gene_transcript, 
                                          gene_tpms_df = gene_bidir_tpms, 
                                          window = window)
    
    ##-----------------------------------------------------------
    # filter samples that match the tissue of interest
    ##-----------------------------------------------------------
    gene_bidir_sel_matrix  <- t(gene_bidir_sel[ ,colnames(gene_bidir_sel) %in% 
                                               metadata_tissue$sample_name, 
                                               with=FALSE]) 
    colnames(gene_bidir_sel_matrix) <- gene_bidir_sel$gene_transcript
    
    ##-----------------------------------------------------------
    # log transform the matrix of tpms and run cov()
    ##-----------------------------------------------------------
    if (exclude_missing_data == TRUE) {
        
        #log transform
        gene_bidir_sel_matrix_log10 <- log(gene_bidir_sel_matrix+1, base=10)

        #assigned missing values as NAs
        gene_bidir_sel_matrix_log10_NAs <- gene_bidir_sel_matrix_log10
        gene_bidir_sel_matrix_log10_NAs[gene_bidir_sel_matrix_log10_NAs == 0] <- NA

        #run covariance matrix with use = "pairwise.complete.obs"
        gene_bidir_cov_mat <- cov(gene_bidir_sel_matrix_log10_NAs,
                                  use = "pairwise.complete.obs")
        
        #assigned 0 to NAs
        gene_bidir_cov_mat[is.na(gene_bidir_cov_mat)] <- 0
        
        #assigne a diagonal of 1e-4 to ensure stability of GLASSO
        #diag(gene_bidir_cov_mat) <- diag(gene_bidir_cov_mat) + 1e-4
    
        } else {
        
        gene_bidir_sel_matrix_log10 <- log(gene_bidir_sel_matrix+1, base=10)
        gene_bidir_cov_mat <- cov(gene_bidir_sel_matrix_log10)
        diag(gene_bidir_cov_mat) <- diag(gene_bidir_cov_mat) + 1e-4
    }
    
    ##-----------------------------------------------------------
    # get distance matrix
    ##-----------------------------------------------------------
    sig_dist_matrix <- get_distance_matrix(gene_bidir_sel, 
                                           nsamples = nrow(gene_bidir_sel))
    row.names(sig_dist_matrix) <- colnames(sig_dist_matrix) <- rownames(gene_bidir_cov_mat)
    
    ##-----------------------------------------------------------
    # calculate the distance penalty
    ##-----------------------------------------------------------
    gene_bidir_rho_mat <- get_rho_mat(sig_dist_matrix)
    
    ##-----------------------------------------------------------
    # Run Graphical lasso with distance penalty
    ##-----------------------------------------------------------
    GL_gene_bidir <- glasso::glasso(gene_bidir_cov_mat,
                                    gene_bidir_rho_mat)

    # rename row and column names
    colnames(GL_gene_bidir$w) <- row.names(GL_gene_bidir$w) <- row.names(gene_bidir_cov_mat)

    # convert covariance to correlations coefficients
    cors_gene_bidir <- stats::cov2cor(GL_gene_bidir$w)
    
    ##-----------------------------------------------------------
    # restructure matrices
    ##-----------------------------------------------------------
    # correlation matrix
    cors_gene_bidir_dt <- restructure_matrix(cors_gene_bidir)
    colnames(cors_gene_bidir_dt) <- c("transcript_1","transcript_2", "pcc_penalty", "pair_id")
    
    # distance matrix
    dist_gene_bidir_dt <- restructure_matrix(sig_dist_matrix)
    colnames(dist_gene_bidir_dt) <- c("transcript_1","transcript_2", "distance", "pair_id")
    
    # merge the distance and correlation matrix
    cors_dist_gene_bidir_dt <- merge(cors_gene_bidir_dt,
                                    dist_gene_bidir_dt[,c("distance","pair_id")],
                                    by="pair_id")

    cors_dist_pair_dt <- subset(cors_dist_gene_bidir_dt,
                                transcript_1 == gene_name &
                                transcript_2 != gene_name)
    
    return(cors_dist_pair_dt)
    
}

save_pairs <- function(pairs_dt, 
                       output_folder,
                       tissue_name,
                       chromosome_id){
    
    #saving the final dataframes
    #including chromosome ids to the file names
    tissue_name_nospace <- gsub(" ", "-", tissue_name)

    final_path <- paste0(output_folder,
                         'glasso_distance_penalty_',
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

##0 : start run Log
date_time <- format(Sys.time(), "%Y_%B_%d_%H_%M_%S")
sink(paste0(output_folder,"log_glasso_",date_time,".txt"))
cat("Running glasso with distance penalty")
cat("\n")

cat(paste0("Date: ", Sys.Date()))
cat("\n")

print("START")
cat("\n")

##1: get chromosome genes
genes_bidir_chr_df <- subset(tpms_datatable, chrom==opt$chr_id) #get the chromosomes genes and bidirectionals
genes_chr_list <- genes_bidir_chr_df[grepl("_", genes_bidir_chr_df$gene_transcript),]$gene_transcript #get list of genes on chromosome

##2: subset the tissues metadata
metadata_celltype <- subset(metadata, tissue==opt$tissue) 

##3: process and run glasso on all genes in chromosome
glasso_pairs_list <- mclapply(genes_chr_list,
                            process_glasso, 
                            gene_bidir_tpms = genes_bidir_chr_df,
                            metadata_tissue = metadata_celltype,
                            window = opt$window,
			    exclude_missing_data = opt$exclude_missing_data, 	
                            mc.cores = opt$ncores)

##4: merge all the pairs
glasso_pairs_merged <- do.call(rbind, glasso_pairs_list)

##5: save the pairs as a table
save_pairs(pairs_dt = glasso_pairs_merged,
           output_folder = output_folder,
           tissue_name = opt$tissue,
           chromosome_id = opt$chr_id)

cat("Session Summary")
cat("\n")

print(sessionInfo())
cat("\n")

cat("DONE!")
cat("\n")

sink()
