#! /usr/bin/env Rscript

##################
##load packages ##
##################

#suppressMessages(library(WGCNA)) ## faster cor()
suppressMessages(library(dplyr)) ## for the R pipes
suppressMessages(library(tidyr)) ## for tidying the dataframes
suppressMessages(library(data.table)) ## load files into R faster
#suppressMessages(library(parallel)) ## running code in parallel
suppressMessages(library(optparse)) ## adding arguments 
#suppressMessages(library(glasso)) ## running glasso on pairs
suppressMessages(library(cglasso)) ## running glasso on sparse data
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


process_cglasso <- function(gene_name, 
                           gene_bidir_tpms,
                           metadata_tissue, 
                           window = 1000000, 
                           distance_weight = TRUE){
    
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
    
    #filter lowly transcribed samples
    gene_bidir_sel_matrix_filt <- filter_low_transcribed(gene_bidir_sel_matrix,
                                                         nlimit = 2)

    gene_bidir_sel <- gene_bidir_sel[gene_bidir_sel$gene_transcript %in%
                                    colnames(gene_bidir_sel_matrix_filt),]
    
    ##-----------------------------------------------------------
    # log transform the matrix of tpms and run cov()
    ##-----------------------------------------------------------
        
    #log transform
    gene_bidir_sel_matrix_log10 <- log(gene_bidir_sel_matrix_filt+1, base=10)

    #assigned missing values as NAs
    gene_bidir_sel_matrix_log10_NAs <- gene_bidir_sel_matrix_log10
    gene_bidir_sel_matrix_log10_NAs[gene_bidir_sel_matrix_log10_NAs == 0] <- NA

    #create a cggm data class
    gene_bidir_cggm <- datacggm(Y=gene_bidir_sel_matrix_log10_NAs)
        
    #rename the columns and rows
    #rename rows and columns
    rowNames(gene_bidir_cggm)$Y <- paste0("u", seq_len(nobs(gene_bidir_cggm)))
    colNames(gene_bidir_cggm)$Y <- paste0("Y", seq_len(nresp(gene_bidir_cggm)))
    
    ##-----------------------------------------------------------
    # get distance matrix
    ##-----------------------------------------------------------
    sig_dist_matrix <- get_distance_matrix(gene_bidir_sel, 
                                           nsamples = nrow(gene_bidir_sel))
    row.names(sig_dist_matrix) <- colnames(sig_dist_matrix) <- colnames(gene_bidir_sel_matrix_filt)
    
    ##-----------------------------------------------------------
    # calculate the distance penalty
    ##-----------------------------------------------------------
    gene_bidir_rho_mat <- get_rho_mat(sig_dist_matrix)
    
    ##-----------------------------------------------------------
    # Run Graphical lasso with distance penalty
    ##-----------------------------------------------------------
    if (distance_weight == TRUE){
        
        #run cglasso with weight
        gene_bidir_cglasso <- cglasso(data = gene_bidir_cggm, weights.Tht=gene_bidir_rho_mat)
        
        } else {
        
        #run cglasso with weight
        gene_bidir_cglasso <- cglasso(data = gene_bidir_cggm)
        
    }
    
    # model selection with eBIC
    gene_bidir_cglasso_eBIC <- BIC(gene_bidir_cglasso, g = 0.5, mle = TRUE)
    
    # get optimal fitted model
    gene_bidir_optmdl <- select_cglasso(gene_bidir_cglasso,
                                           GoF = gene_bidir_cglasso_eBIC)
    
    # index for the optimal parameters
    rho_index <- match(gene_bidir_optmdl$rho,
                       gene_bidir_cglasso$rho)
    
    # extract the regression coefficient matrix (type = "B"), 
    # the covariance matrix (type = "Sigma") 
    # or the precision matrix (type = "Theta).
    gene_bidir_cglasso_coeff <- coef(gene_bidir_cglasso, 
                                     type = "all",
                                     drop = FALSE)

    # get the covariance matrix for the optimal parameters
    gene_bidir_cov <- gene_bidir_cglasso_coeff$Sigma[,,,rho_index]
    colnames(gene_bidir_cov) <- rownames(gene_bidir_cov) <- colnames(gene_bidir_sel_matrix_log10)

    # convert covariance to correlations coefficients
    gene_bidir_cors <- stats::cov2cor(gene_bidir_cov)
        
    ##-----------------------------------------------------------
    # restructure matrices
    ##-----------------------------------------------------------
    # correlation matrix
    cors_gene_bidir_dt <- restructure_matrix(gene_bidir_cors)
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
#genes_bidir_chr_df <- subset(tpms_datatable, chrom==opt$chr_id) #get the chromosomes genes and bidirectionals
gene_name <- tpms_datatable[grepl("_", tpms_datatable$gene_transcript),]$gene_transcript #get list of genes on chromosome

##2: subset the tissues metadata
metadata_celltype <- subset(metadata, tissue==opt$tissue) 

##3: process and run glasso on all genes in chromosome
glasso_pairs <- process_cglasso(gene_name=gene_name, 
                                    gene_bidir_tpms=tpms_datatable,
                                    metadata_tissue=metadata_celltype, 
                                    window = opt$window, 
                                    distance_weight = opt$exclude_missing_data)

##4: save the pairs as a table
save_pairs(pairs_dt = glasso_pairs,
           output_folder = output_folder,
           tissue_name = opt$tissue,
           chromosome_id = gene_name)

cat("Session Summary")
cat("\n")

print(sessionInfo())
cat("\n")

cat("DONE!")
cat("\n")

sink()
