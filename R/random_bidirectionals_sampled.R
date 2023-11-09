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
        make_option(c("-t", "--tpm_file"), type="character", default=NULL,
              help="tab separated file with normalized counts.", metavar="character"),        
        make_option(c("-g", "--gene_counts"), type="character", default="NULL",
              help="table with gene counts", metavar="character"),
        make_option(c("-s", "--seed"), type="integer", default=42,
              help="Default [default = seed %default]", metavar="integer"),
        make_option(c("-o", "--out"), type="character", default="./",
              help="path to output directory [default= %default]", metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$tpm_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input correlations).n", call.=FALSE)
}

tpms_DT <- data.table::fread(opt$tpm_file)
genes <- data.table::fread(opt$gene_counts)
chroms <- tpms_DT$chrom

out_directory <- opt$out

swapped_chr_bidirs <- function(tpmDT, chrA, genenames, seed=opt$seed){
    
    set.seed(seed)
    
    #assign trascript type
    tpmDT$transcipt_type <- ifelse(tpmDT$gene_transcript %in% genenames,
                                                 "Gene","Bidirectional")
    
    #get chrA genes and bidirectionals
    tpmDT_chrA_genes <- subset(tpmDT, chrom == chrA & transcipt_type =='Gene')
    tpmDT_chrA_bidirs <- subset(tpmDT, chrom == chrA & transcipt_type =='Gene')
    
    #get non-chrA bidirectionals
    tpmDT_other_chr_bidirs <- subset(tpmDT, chrom != chrA & transcipt_type !='Gene')

    #sample non-chrA bidirectionals
    sample_other_chr_bidirs <- tpmDT_other_chr_bidirs[sample(1:nrow(tpmDT_other_chr_bidirs),
                                                             nrow(tpmDT_chrA_bidirs),
                                                            replace = FALSE),]
    #merge the annotations with sampled bidirectionals
    tpmDT_chrA_bidirs_new <- cbind(tpmDT_chrA_bidirs[,1:6], 
                                   sample_other_chr_bidirs[,7:ncol(sample_other_chr_bidirs)])
    
    #combine the gene and bidirectional counts
    tpmDT_chrA_genes_other_chr_bidirs <- rbind(tpmDT_chrA_genes,
                                               tpmDT_chrA_bidirs_new)
    
    #remove transcipt columns
    tpmDT_chrA_genes_sampled_bidirs <- tpmDT_chrA_genes_other_chr_bidirs[,transcipt_type:=NULL]
    
    print(paste0("Seed   : ", seed))
    print(paste("Genes in: ", 
          as.character(chrA)," = ", 
          as.character(nrow(tpmDT_chrA_genes))))
    
    print(paste("Bidirs in: ", 
          as.character(chrA)," = ", 
          as.character(nrow(tpmDT_chrA_bidirs_new))))
    
    print(paste("Transcripts in Swapped Annotations: ",
                " = ", 
                as.character(nrow(tpmDT_chrA_genes_sampled_bidirs))))
    
    return(tpmDT_chrA_genes_sampled_bidirs)
}


##create a list for the data.tables
swap_list_inter <- list()

##populate the list by chromosome
for (i in 1:length(chroms)){
    
    print(chroms[i])
    
    ##swap the bidirectionals by other chromosome bidirectionals
    swap_list_inter[[i]] <- swapped_chr_bidirs(tpms_DT,
                                              chroms[i], 
                                              genes$gene_transcript)
    
}

##comnbine all transcript counts
swap_inter_DT <- do.call(rbind, swap_list_inter)

data.table::fwrite(swap_inter_DT,
		 paste0(out_directory, "genes_inter_sampled_bidir_by_chr.tsv.gz"),
		 sep='\t')