#! /usr/bin/env Rscript

##################
##load packages ##
##################

suppressMessages(library(data.table)) ## load files into R faster
suppressMessages(library(optparse)) ##adding arguments

#############################################
## Initialize command options for script   ##
#############################################
# define input and output options
option_list = list(
        make_option(c("-t", "--tpm_file"), type="character", default=NULL,
              help="tab separated file with normalized counts.", metavar="character"),        
        make_option(c("-g", "--gene_counts"), type="character", default="NULL",
              help="table with gene counts", metavar="character"),
        make_option(c("-s", "--seed"), type="integer", default=42,
              help="Default [default = seed %default]", metavar="integer"),
        make_option(c("-o", "--out"), type="character", default="./",
              help="path to output directory [default= %default]", metavar="character"),
        make_option(c("-n", "--out_file"), type="character", default="gene_sampled_bidirs",
              help="path to output directory [default= %default]", metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$tpm_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input correlations).n", call.=FALSE)
}

date_time <- format(Sys.time(), "%Y_%B_%d_%H_%M_%S")

tpms_DT <- data.table::fread(opt$tpm_file)
genes <- data.table::fread(opt$gene_counts)
chroms <- unique(tpms_DT$chrom)

out_directory <- opt$out
out_file_name <- opt$out_file

sink_file_name <- base::basename(out_file_name)

sink(paste0(out_directory,"log_filter_",sink_file_name,"_",date_time,".txt"))

cat("Shuffling bidirectional counts across chromosomes")
cat("\n")

cat(paste0("Date                        : ", Sys.Date()))
cat("\n")

print("START")
cat("\n")

print(paste("Input normalized counts    :", as.character(opt$tpm_file)))
cat("\n")

print(paste("Gene counts                :", as.character(opt$gene_counts)))
cat("\n")

print(paste("Seed                       :", as.character(opt$seed)))
cat("\n")

swapped_chr_bidirs <- function(tpmDT, chrA, genenames, seed=opt$seed){
    
    set.seed(seed)
    
    #1. assign trascript type
    tpmDT$transcipt_type <- ifelse(tpmDT$gene_transcript %in% genenames,
                                                 "Gene","Bidirectional")
    
    #2. get chrA genes and bidirectionals
    tpmDT_chrA_genes <- subset(tpmDT, chrom == chrA & transcipt_type =='Gene')
    tpmDT_chrA_bidirs <- subset(tpmDT, chrom == chrA & transcipt_type !='Gene')
    
    #3. get non-chrA bidirectionals
    tpmDT_other_chr_bidirs <- subset(tpmDT, chrom != chrA & transcipt_type !='Gene')

    #4. sample non-chrA bidirectionals
    sample_other_chr_bidirs <- tpmDT_other_chr_bidirs[sample(1:nrow(tpmDT_other_chr_bidirs),
                                                             nrow(tpmDT_chrA_bidirs),
                                                            replace = FALSE),]
    #5. merge the annotations with sampled bidirectionals
    tpmDT_chrA_bidirs_new <- cbind(tpmDT_chrA_bidirs[,1:6], 
                                   sample_other_chr_bidirs[,7:ncol(sample_other_chr_bidirs)])
    
    #6. combine the gene and bidirectional counts
    tpmDT_chrA_genes_other_chr_bidirs <- rbind(tpmDT_chrA_genes,
                                               tpmDT_chrA_bidirs_new)
    
    #7. remove transcipt columns
    tpmDT_chrA_genes_sampled_bidirs <- tpmDT_chrA_genes_other_chr_bidirs[,transcipt_type:=NULL]
    
    print(paste0("Seed   : ", seed))
    cat("\n")

    print(paste("Genes in: ", 
          as.character(chrA)," = ", 
          as.character(nrow(tpmDT_chrA_genes))))
    cat("\n")

    print(paste("Bidirs in: ", 
          as.character(chrA)," = ", 
          as.character(nrow(tpmDT_chrA_bidirs_new))))
    cat("\n")

    print(paste("Total transcripts after shuffled bidirectionals (genes + bidirectionals): ",
                " = ", 
                as.character(nrow(tpmDT_chrA_genes_sampled_bidirs))))
    cat("\n")
    return(tpmDT_chrA_genes_sampled_bidirs)
}


##create a list for the data.tables
swap_list <- list()

##populate the list by chromosome
for (i in 1:length(chroms)){
    
    print(chroms[i])
    
    ##swap the bidirectionals by other chromosome bidirectionals
    swap_list[[i]] <- swapped_chr_bidirs(tpms_DT,
                                              chroms[i], 
                                              genes$gene_transcript)
    
}

##comnbine all transcript counts
swap_DT <- do.call(rbind, swap_list)

data.table::fwrite(swap_DT,
		 paste0(out_directory, out_file_name, ".tsv.gz"),
		 sep='\t')

cat("\n")
print("Session Summary")
cat("\n")

print(sessionInfo())
cat("\n")

print("DONE!")
cat("\n")

sink()
