library(ggplot2)
library(biomaRt)
library(tidyverse)
library(openxlsx)
library(stringr)
library(ggpubr)
library(igraph)
library(Biostrings)
library(pheatmap)

#Write a function to create  reference files with mRNAs used for building target and control sets
#The argument of the function "dir" species the directory where input files are stored  in a folder "targetscan" 
build_reference <- function(dir){
  
  #Import table with predicted conserved targets of miRNA families obtained with TargetScan
  conserved_targets <- read.csv2(paste(dir, "targetscan/Predicted_Targets_Info.txt", sep=""), sep="\t", header=TRUE)
  
  #Keep only targets for Mus musculus
  conserved_targets <- conserved_targets[conserved_targets$Species.ID == 10090,]
  
  #Extract GC contents for the predicted targets using biomaRt
  ensembl <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  gene_ids <- paste(gsub("\\..*","",conserved_targets$Gene.ID),collapse=",")
  
  GC_contents <- getBM(attributes=c("ensembl_gene_id","percentage_gene_gc_content"),filters= 'ensembl_gene_id',values=gene_ids, mart=ensembl)
  conserved_targets$GC_content <- GC_contents$percentage_gene_gc_content[match(gsub("\\..*","",conserved_targets$Gene.ID),GC_contents$ensembl_gene_id )]
  
  #Add UTR length for the conserved targets
  UTR_sequences <- read.table(paste(dir, "targetscan/UTR_Sequences_Mm_with_UTR_length.txt", sep=""), header=T)
  conserved_targets$UTR_length <- UTR_sequences$UTR_length[match(conserved_targets$Transcript.ID, UTR_sequences$ENST)]
  
  #Add phyloP scores
  mm10_bed <- read.csv2(paste(dir,"targetscan/mm10_genes.bed", sep=""), header=FALSE, sep="\t")
  unique_transcript <- as.data.frame(unique(conserved_targets$Transcript.ID))
  colnames(unique_transcript)[1] <- "ENST"
  unique_mm10_bed <- mm10_bed[match(gsub("\\..*","",unique_transcript$ENST),gsub("\\..*","",mm10_bed$V4)),]
  write.table(unique_mm10_bed, file=paste(dir,"targetscan/unique_mm10_genes.bed", sep=""),row.names=F, sep="\t", quote=FALSE, col.names=FALSE)
  
  colnames(mm10_bed) <- c("chrom","chromStart","chromEnd","name","context++_score_percentile","strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
  predicted_targets_phylop <- read.csv2(paste(dir, "targetscan/unique_mm10_phylop_scores.tab", sep=""), sep="\t", header=T)
  colnames(predicted_targets_phylop) <- c("chrom","chromStart","chromEnd","phylop_score")
  
  predicted_targets_phylop$ENST <- mm10_bed$name[match(paste(predicted_targets_phylop$chrom,predicted_targets_phylop$chromStart,predicted_targets_phylop$chromEnd),paste(mm10_bed$chrom, mm10_bed$chromStart, mm10_bed$chromEnd))]
  write.table(predicted_targets_phylop, file=paste(dir,"targetscan/unique_mm10_phylop_scores_geneName.txt", sep=""), row.names=F, quote=F)
  conserved_targets$phylop <- predicted_targets_phylop$phylop_score[match(gsub("\\..*","",conserved_targets$Transcript.ID), gsub("\\..*","",predicted_targets_phylop$ENST))]
  
  #Remove entries with missing values
  conserved_targets <- conserved_targets[complete.cases(conserved_targets),]
  conserved_targets$phylop <- as.numeric(as.character(conserved_targets$phylop))
  
  #Create a data frame with unique targets used for building the control sets
  unique_targets <- conserved_targets[,c(2:4, 12:14)]
  #Remove duplicate entries 
  unique_targets <- unique_targets[-which(duplicated(unique_targets$Gene.Symbol)), ]
  
  
  #Obtain the miRNA families for the miRNAs in the targetscan database
  miR_fam_miRNAs <- read.delim(paste(dir,"targetscan/miR_Family_Info.txt", sep=""))
  miR_fam_miRNAs_mmu <- miR_fam_miRNAs[miR_fam_miRNAs$Species.ID == 10090,]
  write.csv(miR_fam_miRNAs_mmu, file=paste(dir,"targetscan/miRNA_Family_Info.csv", sep=""), row.names=F)
  
  
  return(list(conserved_targets=conserved_targets, unique_targets= unique_targets, miR_family = miR_fam_miRNAs_mmu))
}


#Function to preprocess data for cotargeting analysis
#The arguments of the function are a list of miRNAs and an output directory
miRNA_create_sets <- function(miRNA_table, outDir, seed, deMI_list=NULL){
  
  #Create output directory if it does not exist
  if (file.exists(outDir)==F){
    dir.create(outDir)
  } 
  #Get the miRNA family for the miRNAs provided in the input
  miRNA_table$miR_family <- miR_family$miR.family[match(miRNA_table$miRNA, miR_family$MiRBase.ID)]
  
  #Write a table with the miRNA families for each miRNA 
  write.xlsx(miRNA_table, file=paste0(outDir, "miRNA_family_list.xlsx"), rowNames=F)
  
  #Remove miRNAs without family assignment
  miRNA_table <- miRNA_table[!is.na(miRNA_table$miR_family),]
  
  #Remove "mmu"  pattern
  #miRNA_table$miRNA <- gsub("mmu-", "", miRNA_table$miRNA)
  
  #Revert to original miRNA names without the .1 suffix for some miRNAs
  miRNA_table$miRNA = miRNA_table$miRNA_renamed
  
  #Merge miRNAs which belong to the same conserved family
  families <- unique(miRNA_table$miR_family)
  targets<-vector()
  for (i in 1:length(families)){
    targets[i] <- paste(miRNA_table$miRNA[miRNA_table$miR_family==families[i]], collapse = "_")
  }
  
  miRNA_table <- data.frame(targets, families)
  colnames(miRNA_table) <- c("miRNA", "miR-family")
  
  
  #Create a matrix for the target count results
  target_counts <- data.frame(matrix(nrow = dim(miRNA_table)[1], ncol = 4))
  
  for (i in 1:nrow(miRNA_table)){
    targets <- unique(conserved_targets$Gene.Symbol[miRNA_table[i,2] == conserved_targets$miR.Family])
    target_counts[i, 1] <- miRNA_table[i,1]
    target_counts[i, 2] <- miRNA_table[i,2]
    target_counts[i, 3] <- length(targets)
    target_counts[i, 4] <- paste(targets,collapse=",") 
  }
  
  #Filter out all miRNA that have least amount of conserved target than cutoff	
  conserved_targets_save <- target_counts[target_counts[,3] >= cutoff,]
  
  #Obtain the seed sequence for each miR-family
  for(s in 1:dim(conserved_targets_save)[1]){
    conserved_targets_save[s,5] <- as.character(miR_family[match(conserved_targets_save[s,2], miR_family$miR.family),2])
  }
  
  colnames(conserved_targets_save) <- c("miRNA","miR_family","target_counts", "targets", "seed")
  write.xlsx(conserved_targets_save, file=paste0(outDir,"/conserved_targets_",cutoff,".xlsx"))
  
   

  conserved_targets_final <- conserved_targets_save

  result <- list()
  
  #Create a directory to store files with the target sets for each miRNA 
  dir_targets_sets <- paste0(outDir, "target_sets")
  
  if (file.exists(dir_targets_sets)==F){
    dir.create(dir_targets_sets)
  } 
  
  #Create output directory for control sets
  dir_controls <- paste0(outDir, "/control_sets")
  
  if (file.exists(dir_controls)==F){
    dir.create(dir_controls)
  }
  
  #Create an excel sheet with the target set for each microRNA and the corresponding control sets
  for (k in 1:nrow(conserved_targets_final)){
    target_names <- data.frame(strsplit(conserved_targets_final[k,4],","))
    colnames(target_names) <- "target_name"
    target_UTR_length <- as.data.frame(conserved_targets$UTR_length[match(target_names$target_name, conserved_targets$Gene.Symbol)])
    colnames(target_UTR_length) <- "target_UTR_length" 
    target_GC_content <- as.data.frame(conserved_targets$GC_content[match(target_names$target_name,conserved_targets$Gene.Symbol)])
    colnames(target_GC_content) <- "target_GC_content"
    target_phylop <- as.data.frame(conserved_targets$phylop[match(target_names$target_name,conserved_targets$Gene.Symbol)])
    colnames(target_phylop) <- "target_phylop"
    column.names <- c("target_name", "target_UTR_length", "target_GC_content", "target_phylop")
    matrix.names <- conserved_targets_final[k,1]
    targets_table <- data.frame(target_names,target_UTR_length,target_GC_content, target_phylop)

    write.xlsx(targets_table, paste0(dir_targets_sets,"/", conserved_targets_final$miRNA[k], ".xlsx"))
    
    set.seed(seed)
    
    control <- data.frame()
    
    #Filter out mRNAs from the unique conserved mRNAs which come in the target set of the current miRNA
    mRNAs <- unique_targets[unique_targets$Gene.Symbol %in% targets_table$target_name==F,]
    
    for (j in 1:nrow(targets_table)){
      hits <- mRNAs[mRNAs$UTR_length>0.85*targets_table$target_UTR_length[j] & mRNAs$UTR_length<1.15*targets_table$target_UTR_length[j],]
      hits <- hits[hits$GC_content>0.95*targets_table$target_GC_content[j] & hits$GC_content<1.05*targets_table$target_GC_content[j],]
      hits <- hits[hits$phylop>0.8*targets_table$target_phylop[j] & hits$phylop<1.2*targets_table$target_phylop[j],]
      
      r <- sample(1:nrow(hits),1)
      row <- hits[r,]
      control <- rbind(control, row)
      mRNAs <- mRNAs[mRNAs$Gene.Symbol %in% row$Gene.Symbol == F,]
      
    }
    control <- control[complete.cases(control),]
    write.xlsx(control, file=paste0(dir_controls, "/", conserved_targets_final$miRNA[k], "_control.xlsx"))
    
  }
  
  
}

miRNA_cotargeting <- function(outDir, dir_targets, dir_controls, miRNA_info){
  #Import target sets
  setwd(dir_targets)
  
  temp <- list.files(pattern = "*.xlsx")
  names <- str_remove(temp, pattern = ".xlsx")
  #names <- gsub("-", "_", names)
  
  for (i in 1:length(temp)) {
    assign(names[i], read.xlsx(temp[i], colNames = TRUE))
  }
  
  ####################################################################################
  #Import the control sets
  setwd(dir_controls)
  temp <- list.files(pattern = "*.xlsx")
  names_control <- str_remove(temp, pattern = ".xlsx")
  
  for (i in 1:length(temp)) {
    assign(names_control[i], read.xlsx(temp[i], colNames = TRUE))
  }
  
  ######################################################################################
  #Compare control and target sets in terms of UTR length, GC content and phylop score
  
  comparison <- matrix(nrow = length(names), ncol = 13, 0)
  colnames(comparison) <- c("control_set", "target_set", "number_controls", 
                            "number_targets", "avg_utr_ctr", "avg_utr_target", "p_utr",
                            "avg_gc_ctr", "avg_gc_target", "p_gc", "avg_phylop_ctr",
                            "avg_phylop_target", "p_phylop")
  
  
  for (i in 1:length(names)){
    ctr = get(names_control[i])
    target = get(names[i])
    
    comparison[i, 1] <- names_control[i]
    comparison[i,2] <- names[i]
    comparison[i,3] <- nrow(ctr)
    comparison[i,4] <- nrow(target)
    comparison[i,5] <- round(mean(ctr$UTR_length),2)
    comparison[i,6] <- round(mean(target$target_UTR_length),2)
    comparison[i,7] <- round(wilcox.test(ctr$UTR_length, target$target_UTR_length)$p.value,3)
    comparison[i,8] <- round(mean(ctr$GC_content),2)
    comparison[i,9] <- round(mean(target$target_GC_content),2)
    comparison[i,10] <- round(wilcox.test(ctr$GC_content, target$target_GC_content)$p.value,3)
    comparison[i,11] <- round(mean(ctr$phylop),2)
    comparison[i,12] <- round(mean(target$target_phylop),2)
    comparison[i,13] <- round(wilcox.test(ctr$phylop, target$target_phylop)$p.value,3)
    
  }
  
  #Export file
  write.xlsx(comparison, file = paste0(outDir, "/comparison_taget_control.xlsx"))
  
  
  ###############################################################################################
  #Sanity checks
  #Make sure that mRNA in the control set do not appear in the target set of a miRNA
  
  for (i in 1:length(names)){
    target <- get(names[i])
    ctr <- get(names_control[i])
    
    print(which(target$target_name %in% ctr$Gene.Symbol))
  }
  
  ###############################################################################################
  #Perform Fisher's exact tests to investigate if there are miRNA co-targeting pairs
  
  #Create a results matrix
  mirna_pairs <- matrix(nrow = length(names)*(length(names)-1)/2, ncol = 6)
  colnames(mirna_pairs) <- c("OR_A_B", "p_A_B", "adj.p_A_B", "OR_B_A", "p_B_A", "adj.p_B_A")
  
  mirna_pairs_names <- matrix(nrow=length(names)*(length(names)-1)/2, ncol=3)
  colnames(mirna_pairs_names) <- c("mirnaA", "mirnaB", "common_targets")
  
  chi.res <- matrix(nrow=length(names)*(length(names)-1)/2, ncol = 2)
  colnames(chi.res) <- c("p_A_B", "p_B_A") 
  
  ind=1
  for (i in 1:(length(names)-1)){
    targetA <- get(names[i])
    ctrA <- get(names_control[i])
    
    for (j in (i+1):length(names)){
        
    #Calculate co-targeting relationship if both miRNAs have distinct sequences
      seed1 = miRNA_info$seed[miRNA_info$miRNA == names[i]] 
      seed2 = miRNA_info$seed[miRNA_info$miRNA == names[j]] 
      
      x = Biostrings::matchPattern(seed1,seed2,max.mismatch=1)
      
      if(length(x)<1){
       
          targetB <- get(names[j])
          ctrB <- get(names_control[j])
          
          #miRNA A as reference (using control set of miRNA A)
          common_targets_ind <- which(targetA$target_name %in% targetB$target_name == TRUE)
          common_targets <- targetA$target_name[common_targets_ind]
          target_ab <- length(common_targets)
          ctr_ab <- length(which(ctrA$Gene.Symbol %in% targetB$target_name == TRUE))
          
          tab1 <- rbind(c(target_ab, nrow(targetA)-target_ab), c(ctr_ab, nrow(ctrA)-ctr_ab))
          f1 <- fisher.test(tab1)
          chi1 <- chisq.test(tab1)
          #miRNA B as reference (using control set of miRNA B)
          ctr_ba <- length(which(ctrB$Gene.Symbol %in% targetA$target_name == TRUE))
          
          tab2 <- rbind(c(target_ab, nrow(targetB)- target_ab), c(ctr_ba, nrow(ctrB) - ctr_ba))
          f2 <- fisher.test(tab2)
          chi2 <- chisq.test(tab2)
          #Write output to the results table
          mirna_pairs_names[ind,1] <- names[i]
          mirna_pairs_names[ind, 2] <- names[j]
          mirna_pairs_names[ind,3] <- paste(common_targets, collapse=",")
          mirna_pairs[ind,1] <- round(f1$estimate,3)
          mirna_pairs[ind, 2] <- f1$p.value
          mirna_pairs[ind, 4] <- round(f2$estimate, 3)
          mirna_pairs[ind, 5] <- f2$p.value
          chi.res[ind,1] <- chi1$p.value
          chi.res[ind,2] <- chi2$p.value  
          ind <- ind + 1
      }
      
    }
  }
  
  #Adjust p-values for multiple comparisons
  mirna_pairs_res <- data.frame(mirna_pairs_names, mirna_pairs)

  
  mirna_pairs_res$adj.p_A_B <- p.adjust(mirna_pairs_res$p_A_B, method = "BH")
  mirna_pairs_res$adj.p_B_A <- p.adjust(mirna_pairs_res$p_B_A, method = "BH")
  
  
  #Remove NA values 
  
  mirna_pairs_res = na.omit(mirna_pairs_res)
  
  #Calculate average log-odds ratio
  mirna_pairs_res$avg_log_OR <- rowMeans(data.frame(log(as.numeric(mirna_pairs_res$OR_A_B)), log(as.numeric(mirna_pairs_res$OR_B_A))))
  mirna_pairs_res$avg_log_pval <- rowMeans(data.frame(log(as.numeric(mirna_pairs_res$adj.p_A_B)), log(as.numeric(mirna_pairs_res$adj.p_B_A))))
  bidirectional_pairs <- mirna_pairs_res[mirna_pairs_res$adj.p_A_B<0.05 & mirna_pairs_res$adj.p_B_A<0.05,]
  write.xlsx(bidirectional_pairs[order(bidirectional_pairs$adj.p_A_B),], paste0(outDir, "/bidirectional_pairs.xlsx"))
  
  unidirectional_pairs <- mirna_pairs_res[mirna_pairs_res$adj.p_A_B<0.05 | mirna_pairs_res$adj.p_B_A<0.05,]
  write.xlsx(unidirectional_pairs[order(unidirectional_pairs$adj.p_A_B),], paste0(outDir, "/unidirectional_pairs.xlsx"))
  
  
  #Modify results matrices so that miRNAs are not combined to families
  #The tables are just collapsed meaning that miRNAs that belong to the same
  #family have identical co-targeting relationships and statistics, they are just listed separately
  
  bidirectional_pairs_collapsed = data.frame()
  
  for (i in 1:nrow(bidirectional_pairs)){
      miRNAs_A = unlist(strsplit(bidirectional_pairs$mirnaA[i], split = "_"))
      miRNAs_B = unlist(strsplit(bidirectional_pairs$mirnaB[i], split = "_"))
      
      tab1 = expand.grid(miRNAs_A, miRNAs_B)
      tab2 = do.call("rbind", replicate(nrow(tab1), bidirectional_pairs[i,3:11], simplify = FALSE))
      
      tab = cbind(tab1, tab2)
      
      bidirectional_pairs_collapsed = rbind(bidirectional_pairs_collapsed, tab)
      
      if(any(is.na(miRNAs_A))){
          print(i)
      }
  }
  
  colnames(bidirectional_pairs_collapsed) = colnames(bidirectional_pairs)
  bidirectional_pairs_collapsed[,1:2] = sapply(bidirectional_pairs_collapsed[,1:2],
                                               as.character)
  
  
  #---------------------------------------------
  #Unidirection pairs 
  unidirectional_pairs_collapsed = data.frame()
  
  for (i in 1:nrow(unidirectional_pairs)){
      miRNAs_A = unlist(strsplit(unidirectional_pairs$mirnaA[i], split = "_"))
      miRNAs_B = unlist(strsplit(unidirectional_pairs$mirnaB[i], split = "_"))
      
      tab1 = expand.grid(miRNAs_A, miRNAs_B)
      tab2 = do.call("rbind", replicate(nrow(tab1), unidirectional_pairs[i,3:11], simplify = FALSE))
      
      tab = cbind(tab1, tab2)
      
      unidirectional_pairs_collapsed = rbind(unidirectional_pairs_collapsed, tab)
      
      
  }
  
  colnames(unidirectional_pairs_collapsed) = colnames(unidirectional_pairs)
  unidirectional_pairs_collapsed[,1:2] = sapply(unidirectional_pairs_collapsed[,1:2],
                                               as.character)
  
  return(list(all_pairs = mirna_pairs_res, bidirectional = bidirectional_pairs, unidirectional = unidirectional_pairs,
              bidirectional_collapsed = bidirectional_pairs_collapsed,  unidirectional_collapsed = unidirectional_pairs_collapsed))
  
}

################################################################################
#FUNCTION TO CREATE adjacency MATRICES FOR PLOTTING OF RESULTS

create_adjacency_matrix = function(cotargeting_res){
    
    #Get names of miRNAs
    names = unique(c(cotargeting_res$mirnaA, cotargeting_res$mirnaB))
    
    ####################################
    #Create an adjacency matrix for cotargeting pairs
    adjacency_mat <- matrix(nrow=length(names), ncol=length(names), 0)
    rownames(adjacency_mat) <- names
    colnames(adjacency_mat) <- names
    
    for (i in 1:length(names)){
        for (j in 1:nrow(cotargeting_res)){
            if (names[i] == cotargeting_res$mirnaA[j]){
                partner <- cotargeting_res$mirnaB[j]
                ind <- which(colnames(adjacency_mat)==partner)
                adjacency_mat[i, ind] <- 1
                adjacency_mat[ind,i ] <- 1
            }
        }
    }
    
    ####################################
    #Create an adjacency matrix with odds ratios (all pairs irrespective of p-value)
    adjacency_mat_OR <- matrix(nrow=length(names), ncol=length(names), 0)
    rownames(adjacency_mat_OR) <- names
    colnames(adjacency_mat_OR) <- names
    
    for (i in 1:length(names)){
        for (j in 1:nrow(cotargeting_res)){
            if (names[i] == cotargeting_res$mirnaA[j]){
                partner <- cotargeting_res$mirnaB[j]
                ind <- which(colnames(adjacency_mat_OR)==partner)
                adjacency_mat_OR[i, ind] <- cotargeting_res$avg_log_OR[j]
                adjacency_mat_OR[ind,i ] <- cotargeting_res$avg_log_OR[j]
            }
        }
    }
    
    ####################################
    #Create an adjacency matrix with average log p-value 
    adjacency_mat_pval <- matrix(nrow=length(names), ncol=length(names), 0)
    rownames(adjacency_mat_pval) <- names
    colnames(adjacency_mat_pval) <- names
    
    for (i in 1:length(names)){
        for (j in 1:nrow(cotargeting_res)){
            if (names[i] == cotargeting_res$mirnaA[j]){
                partner <- cotargeting_res$mirnaB[j]
                ind <- which(colnames(adjacency_mat_pval)==partner)
                adjacency_mat_pval[i, ind] <- cotargeting_res$avg_log_pval[j]
                adjacency_mat_pval[ind,i ] <- cotargeting_res$avg_log_pval[j]
            }
        }
    }
    
    #If multiple miRNAs have been combined because they belong to the same family
    #Write each miRNA on a new line
    #colnames(adjacency_mat) <- gsub("_", "\n", colnames(adjacency_mat))
    #rownames(adjacency_mat) <- gsub("_", "\n", rownames(adjacency_mat))
    
    #Remove the "miR-" pattern from the miRNA names to shorten them
    colnames(adjacency_mat) <- gsub("miR-", "", colnames(adjacency_mat))
    rownames(adjacency_mat) <- gsub("miR-", "", rownames(adjacency_mat))
    
    colnames(adjacency_mat_OR) <- gsub("miR-", "", colnames(adjacency_mat_OR))
    rownames(adjacency_mat_OR) <- gsub("miR-", "", rownames(adjacency_mat_OR))
    
    return(list(adjacency_mat = adjacency_mat, adjacency_OR = adjacency_mat_OR, adjacency_pval = adjacency_mat_pval))
}

################################################################################
#ACTUAL ANALYSIS USING THE FUNCTONS DEFINED ABOVE

#Generate seeds for reproducibility when creating the control sets
set.seed(439530)
seeds <- sample(1:1000000000, size=30, replace =F)


#Run the function to obtain reference sets for building target and control sets
dir = "input_files"
reference <- build_reference(dir)
conserved_targets <- reference$conserved_targets
unique_targets <- reference$unique_targets
miR_family <- reference$miR_family
cutoff=300

#################################################################
##All modules processed together
#################################################################
#Run the mirna cotargeting function  using the miRNAs from modules
#Significantly correlated with developmental stage
#See module-trait relationship plot

#Perform the miRNA cotargeting analysis with all modules significantly correlated with developmental stage
sig_modules = c("black", "green", "lightgreen", "magenta", "tan", "yellow", "darkgreen", "darkred", "midnightblue", "lightcyan")

miRNAs_cotargeting = WCGNA_miRNAs[WCGNA_miRNAs$module %in% sig_modules,]
miRNAs_cotargeting$original_edited2 = gsub("mmu-", "", miRNAs_cotargeting$original)

#Change the following parameter to the desired directory
outDir="/output/directory"


#Edit input table for the co-targeting function
miRNAs_cotargeting = data.frame(miRNA = miRNAs_cotargeting$original_edited,
                                miRNA_original = miRNAs_cotargeting$original,
                                miRNA_renamed = miRNAs_cotargeting$original_edited2,
                                miRNA_renamed2 = miRNAs_cotargeting$miRNA,
                                module = miRNAs_cotargeting$module)


#Run cotargeting analysis
#miRNA_create_sets(miRNAs_cotargeting, outDir, seeds[7], deMI_list=miRNA_list_interest)
miRNA_create_sets(miRNAs_cotargeting, outDir, seeds[7])

miRNA_info = read.xlsx("input_files/conserved_targets_300.xlsx", colNames = TRUE)

all_modules_res <- miRNA_cotargeting(outDir, paste0(outDir, "target_sets"), paste0(outDir, "control_sets"), miRNA_info = miRNA_info)

#Create adjecency matrices
adjacency_matrices = create_adjacency_matrix(all_modules_res$bidirectional_collapsed)

###############################################################################
#Get number of co-targeting relationships for each miRNA

##############################################################
#Get number of bidirectional co-targeting relationships for each miRNA
num_cotarget_relationships <- data.frame(table(c(all_modules_res$bidirectional_collapsed$mirnaA, all_modules_res$bidirectional_collapsed$mirnaB)))


num_cotarget_relationships = merge(miRNAs_cotargeting, num_cotarget_relationships,
                                   by.x = "miRNA_renamed", by.y = "Var1", all.y =TRUE)

num_cotarget_relationships <- num_cotarget_relationships[order(num_cotarget_relationships$Freq, decreasing = T),]

#Create histogram with number of cotargeting relationships per miRNA


#Only green and black module
p = ggplot(num_cotarget_relationships[num_cotarget_relationships$module %in% c("green", "black"),], 
       aes(x = Freq, fill = module, col = module)) +
    geom_histogram(position = position_dodge(), alpha = 0.7, show.legend = F) +
    scale_fill_manual(values = c("black", "#2CA05A")) +
    scale_color_manual(values = c("black", "#2CA05A")) +
    theme_test(base_size = 19) + theme(axis.text = element_text(color="black")) +
    labs(x ="No. significant co-targeting relationships", y = "No. of miRNAs") +
    scale_y_continuous(expand = c(0,0), limits = c(0,8))
    



#Histogram with all modules
p = ggplot(num_cotarget_relationships, 
           aes(x = Freq)) +
    geom_histogram(fill = "gray70", col = "gray30") +  
    theme_test(base_size = 19) + theme(axis.text = element_text(color="black")) +
    labs(x ="No. significant co-targeting relationships", y = "No. of miRNAs") +  
    scale_y_continuous(expand = c(0,0), limits = c(0,10)) 
    

#-----------------------------------------------------------------------------------
#Create a plot with number of miRNAs per module with significant co-targeting relationships
miRNAs_module = data.frame(table(num_cotarget_relationships$module))
miRNAs_module = miRNAs_module[order(miRNAs_module$Freq, decreasing = T),]
miRNAs_module$Var1 = factor(miRNAs_module$Var1, levels = miRNAs_module$Var1)
miRNAs_module$color = as.character(miRNAs_module$Var1)
miRNAs_module$color[miRNAs_module$color =="green"] = "#2CA05A"
miRNAs_module$color[miRNAs_module$color =="yellow"] = "#FFDD55"

p = ggplot(miRNAs_module, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", fill = miRNAs_module$color, color ="black") + 
    theme_test(base_size=19) +
    theme(axis.text = element_text(color= "black"),
          axis.text.x = element_text(angle=35, hjust=1)) +
    labs(x ="",  y= "No. of miRNAs") +
    scale_y_continuous(expand = c(0,0), limits = c(0,40))
    



################################################################################
##################################################################
#Create a heatmap with co-targeting relationships

cols <- hcl.colors(n=5, palette = "Teal")

annotation_mat <- data.frame(miRNA = rownames(adjacency_matrices$adjacency_mat))
annotation_mat = merge(annotation_mat, miRNAs_cotargeting[,c("miRNA_renamed2", "module")],
                       by.x = "miRNA", by.y = "miRNA_renamed2", all.x = TRUE)
rownames(annotation_mat) = annotation_mat$miRNA

#Reorder to match the row names of the adjacency matrix
annotation_mat = annotation_mat[match(rownames(adjacency_matrices$adjacency_mat), rownames(annotation_mat)),]
annotation_mat = annotation_mat[2]

annot_cols = list(module = c(black = "black", green = "#2CA05A", darkgreen="darkgreen",
                  tan = "tan", yellow = "#FFDD55", lightgreen = "lightgreen",
                  darkred = "darkred", magenta = "magenta", 
                  lightcyan = "lightcyan", midnightblue = "midnightblue"))



ph <- pheatmap::pheatmap(adjacency_matrices$adjacency_mat,
               fontsize = 8, treeheight_col = 0, legend = F,
               color = c("white", "#3B809A"),
               annotation_row = annotation_mat,
               annotation_names_row = F, 
               annotation_legend_labels = "",
               annotation_colors = annot_cols)



#-------------------------------------------------------------------
#Create a heatmap with average odds ratios
cols <- hcl.colors(n=100, palette = "BluYl")

ph <- pheatmap::pheatmap(adjacency_matrices$adjacency_OR,
               fontsize = 8, treeheight_col = 0,
               color = c("white", rev(cols)),
               annotation_row = annotation_mat,
               annotation_names_row = F, 
               annotation_legend_labels = "",
               annotation_colors = annot_cols)

#-------------------------------------------------------------------
#Create a heatmap with miRNAs grouped by family

adjacency_matrices_fam = create_adjacency_matrix(all_modules_res$bidirectional)

rownames(adjacency_matrices_fam$adjacency_mat) = gsub("_", "/",
                rownames(adjacency_matrices_fam$adjacency_mat))

colnames(adjacency_matrices_fam$adjacency_mat) = gsub("_", "/",
                    colnames(adjacency_matrices_fam$adjacency_mat))


ph <- pheatmap::pheatmap(adjacency_matrices_fam$adjacency_mat,
               fontsize = 8, treeheight_col = 0, legend = F,
               color = c("white", cols[2]))



################################################################################
#       CREATE GRAPH WITH COTARGETING RELATIONSHIPS
g <- graph.adjacency(adjacency_matrices$adjacency_mat, mode="undirected")

hs <-hub_score(g, weights=NA)$vector


g_nodes <- data.frame(nodes = V(g)$name)
g_nodes = merge(g_nodes, miRNAs_cotargeting[,c("miRNA_renamed2", "module")], 
                by.x = "nodes", 
                by.y = "miRNA_renamed2", all.x = TRUE)
g_nodes$color = ifelse(g_nodes$module != "green", g_nodes$module, "#2CA05A")
g_nodes$color = ifelse(g_nodes$module== "yellow", "#FFDD55", g_nodes$color)

#Reorder data frame with nodes and colors to match order of nodes in the graph 
#as this is changes after merging the tables and then color assignment is not correct

g_nodes = g_nodes[match(V(g)$name, g_nodes$nodes),]
g_nodes$color_lab = ifelse(g_nodes$module =="black", "gray40", "black")


V(g)$color <- g_nodes$color

plot_dendrogram(cluster_fast_greedy(g), colbar=safe_palette[1:4], cex=0.7)


par(mar=c(0,0,0,0))
plot(g, vertex.size=hs*15, vertex.label.cex=1.5, vertex.label.font=2,
     vertex.label.color=g_nodes$color_lab, layout=layout_with_graphopt(g, charge=1.2))



################################################################################
#  CREATE A NETWORK ONLY WITH BLACK MODULE RELATIONSHIPS
#-------------------------------------------------------------------------------
cotargeting_black = all_modules_res$bidirectional_collapsed

#Subset to relationships containing miRNAs from the black module

mirnas_black = miRNAs_cotargeting$miRNA_renamed[miRNAs_cotargeting$module == "black"]

cotargeting_black = cotargeting_black[cotargeting_black$mirnaA %in% mirnas_black |
                                          cotargeting_black$mirnaB %in% mirnas_black,]


#Create an adjacency matrix for black module miRNAs
adjacency_matrices_black = create_adjacency_matrix(cotargeting_black)

#---------------------------
#Create the graph 
g <- graph.adjacency(adjacency_matrices_black$adjacency_mat, mode="undirected")

hs <-hub_score(g, weights=NA)$vector


g_nodes <- data.frame(nodes = V(g)$name)
g_nodes = merge(g_nodes, miRNAs_cotargeting[,c("miRNA_renamed2", "module")], by.x = "nodes", 
                by.y = "miRNA_renamed2", all.x = TRUE)
g_nodes$color = ifelse(g_nodes$module != "green", g_nodes$module, "#2CA05A")
g_nodes$color = ifelse(g_nodes$module== "yellow", "#FFDD55", g_nodes$color)

#Reorder data frame with nodes and colors to match order of nodes in the graph 
#as this is changes after merging the tables and then color assignment is not correct

g_nodes = g_nodes[match(V(g)$name, g_nodes$nodes),]


V(g)$color <- g_nodes$color


par(mar=c(0,0,0,0))

plot(g, vertex.size=hs*20, vertex.label.cex=1.5, vertex.label.font=2,
     vertex.label.color="gray20", layout=layout_with_graphopt(g, charge=1.2))


###########################################################################
#   CREATE A CIRCOS PLOT WITH CO-TARGETING BLACK MODULE ONLY
#--------------------------------------------------------------------------

miRNAs_circos = data.frame(miRNA = unique(c(cotargeting_black$mirnaA, 
                                            cotargeting_black$mirnaB)))


miRNAs_circos = merge(miRNAs_circos, miRNAs_cotargeting[,c("miRNA_renamed", "module")],
                      by.x ="miRNA", by.y ="miRNA_renamed", all.x = TRUE)

miRNAs_circos$color = miRNAs_circos$module
miRNAs_circos$color = ifelse(miRNAs_circos$color == "green", "#2CA05A",
                             miRNAs_circos$color)
miRNAs_circos$color = ifelse(miRNAs_circos$color == "yellow", "#FFDD55",
                             miRNAs_circos$color)


#-----------
miRNAs_circos = miRNAs_circos[order(miRNAs_circos$module),]
miRNAs_circos$miRNA = gsub("miR-", "", miRNAs_circos$miRNA)

#Format the cotargeting_black table to include the module and the color for the link
cotargeting_black$moduleA = NA
cotargeting_black$moduleB = NA
cotargeting_black$color_link = NA

for (i in 1:nrow(cotargeting_black)){
    
    cotargeting_black$moduleA[i] = miRNAs_cotargeting$module[miRNAs_cotargeting$miRNA_renamed == 
                                                  cotargeting_black$mirnaA[i]]
    cotargeting_black$moduleB[i] = miRNAs_cotargeting$module[miRNAs_cotargeting$miRNA_renamed == 
                                                  cotargeting_black$mirnaB[i]]
    
    if(cotargeting_black$moduleA[i] == "black" & 
       cotargeting_black$moduleB[i] == "black"){
        cotargeting_black$color_link[i] = "black"
            
    }
    
    if(cotargeting_black$moduleA[i] != "black"){
        cotargeting_black$color_link[i] = cotargeting_black$moduleA[i]
    }
    
    if(cotargeting_black$moduleB[i] != "black"){
        cotargeting_black$color_link[i] = cotargeting_black$moduleB[i]
    }
    
}


cotargeting_black$color_link = ifelse(cotargeting_black$color_link == "green",
                                      "#2CA05A", cotargeting_black$color_link)
                                      
cotargeting_black$color_link = ifelse(cotargeting_black$color_link == "yellow",
                                      "#FFDD55", cotargeting_black$color_link)


cotargeting_black = cotargeting_black[order(cotargeting_black$color_link),]
cotargeting_black$mirnaA = gsub("miR-", "", cotargeting_black$mirnaA)
cotargeting_black$mirnaB = gsub("miR-", "", cotargeting_black$mirnaB)

#Create a subtable with only black module cotargeting relationships
cotargeting_black_only = cotargeting_black[cotargeting_black$color_link == "black",]

#Create the circos plot
par(mar=c(0,0,0,0))
circos.clear()
circos.par(gap.after = 1)
circos.initialize(sectors = miRNAs_circos$miRNA, xlim = c(0,1))
circos.track(sectors = miRNAs_circos$miRNA, ylim = c(0,1),
             bg.col = miRNAs_circos$color, track.height=0.26)

for (i in 1:nrow(cotargeting_black)){
    circos.link(cotargeting_black$mirnaA[i], 0.5, cotargeting_black$mirnaB[i], 
                0.5, col = cotargeting_black$color_link[i])
    
}

#Replot the black module only cotargeting relationships
for (i in 1:nrow(cotargeting_black_only)){
    circos.link(cotargeting_black_only$mirnaA[i], 0.5, cotargeting_black_only$mirnaB[i], 
                0.5, col = cotargeting_black_only$color_link[i])
    
}

circos.track(track.index = 1, panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                col = "white",
                facing = "clockwise", cex=0.8, niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

dev.off()

#-------------------------------------------------------------------------------
# CREATE A CHORD PLOT WITH BLACK MODULE COTARGETING RELATIONSHIPS
#-------------------------------------------------------------------------------

#Create a data frame with frequencies of module combinations
freq_table = as.data.frame(table(cotargeting_black[,c("moduleA", "moduleB")]))

freq_table = freq_table[freq_table$Freq!=0,]

#Sum frequencies of relationships between both modules independent of order

freq_table_formatted = freq_table[1:10,]
freq_table_formatted$Freq2 = c(0, freq_table$Freq[11:19])
freq_table_formatted$Sum_Freq = freq_table_formatted$Freq + freq_table_formatted$Freq2

#Column with colors for chord plot
freq_table_formatted$color = as.character(freq_table_formatted$moduleA)
freq_table_formatted$color = ifelse(freq_table_formatted$color == "green",
                                    "#2CA05A", freq_table_formatted$color)
freq_table_formatted$color = ifelse(freq_table_formatted$color == "yellow",
                                  "#FFDD55", freq_table_formatted$color)
                                    

################################################################################
################################################################################
#  CREATE A NETWORK ONLY WITH GREE MODULE RELATIONSHIPS
#-------------------------------------------------------------------------------
cotargeting_green = all_modules_res$bidirectional_collapsed

#Subset to relationships containing miRNAs from the green module

mirnas_green = miRNAs_cotargeting$miRNA_renamed[miRNAs_cotargeting$module == "green"]

cotargeting_green = cotargeting_green[cotargeting_green$mirnaA %in% mirnas_green |
                                          cotargeting_green$mirnaB %in% mirnas_green,]


#Create an adjacency matrix for green module miRNAs
adjacency_matrices_green = create_adjacency_matrix(cotargeting_green)

#---------------------------
#Create the graph 
g <- graph.adjacency(adjacency_matrices_green$adjacency_mat, mode="undirected")

hs <-hub_score(g, weights=NA)$vector


g_nodes <- data.frame(nodes = V(g)$name)
g_nodes = merge(g_nodes, miRNAs_cotargeting[,c("miRNA_renamed2", "module")], by.x = "nodes", 
                by.y = "miRNA_renamed2", all.x = TRUE)
g_nodes$color = ifelse(g_nodes$module != "green", g_nodes$module, "#2CA05A")
g_nodes$color = ifelse(g_nodes$module== "yellow", "#FFDD55", g_nodes$color)

#Reorder data frame with nodes and colors to match order of nodes in the graph 
#as this is changes after merging the tables and then color assignment is not correct

g_nodes = g_nodes[match(V(g)$name, g_nodes$nodes),]


V(g)$color <- g_nodes$color


par(mar=c(0,0,0,0))

plot(g, vertex.size=hs*20, vertex.label.cex=1.5, vertex.label.font=2,
     vertex.label.color="gray10", layout=layout_with_graphopt(g, charge=1.2))


###########################################################################
#   CREATE A CIRCOS PLOT WITH CO-TARGETING BLACK GREEN ONLY
#--------------------------------------------------------------------------

miRNAs_circos = data.frame(miRNA = unique(c(cotargeting_green$mirnaA, 
                                            cotargeting_green$mirnaB)))

miRNAs_circos = merge(miRNAs_circos, miRNAs_cotargeting[,c("miRNA_renamed", "module")],
                      by.x ="miRNA", by.y ="miRNA_renamed", all.x = TRUE)

miRNAs_circos$color = miRNAs_circos$module
miRNAs_circos$color = ifelse(miRNAs_circos$color == "green", "#2CA05A",
                             miRNAs_circos$color)
miRNAs_circos$color = ifelse(miRNAs_circos$color == "yellow", "#FFDD55",
                             miRNAs_circos$color)


#-----------
miRNAs_circos = miRNAs_circos[order(miRNAs_circos$module),]
miRNAs_circos$miRNA = gsub("miR-", "", miRNAs_circos$miRNA)

#Format the cotargeting_green table to include the module and the color for the link
cotargeting_green$moduleA = NA
cotargeting_green$moduleB = NA
cotargeting_green$color_link = NA

for (i in 1:nrow(cotargeting_green)){
    
    cotargeting_green$moduleA[i] = miRNAs_cotargeting$module[miRNAs_cotargeting$miRNA_renamed == 
                                                          cotargeting_green$mirnaA[i]]
    cotargeting_green$moduleB[i] = miRNAs_cotargeting$module[miRNAs_cotargeting$miRNA_renamed == 
                                                          cotargeting_green$mirnaB[i]]
    
    if(cotargeting_green$moduleA[i] == "green" & 
       cotargeting_green$moduleB[i] == "green"){
        cotargeting_green$color_link[i] = "#2CA05A"
            
    }
    
    if(cotargeting_green$moduleA[i] != "green"){
        cotargeting_green$color_link[i] = cotargeting_green$moduleA[i]
    }
    
    if(cotargeting_green$moduleB[i] != "green"){
        cotargeting_green$color_link[i] = cotargeting_green$moduleB[i]
    }
    
}



                                      
cotargeting_green$color_link = ifelse(cotargeting_green$color_link == "yellow",
                                      "#FFDD55", cotargeting_green$color_link)
                                      

cotargeting_green = cotargeting_green[order(cotargeting_green$color_link),]
cotargeting_green$mirnaA = gsub("miR-", "", cotargeting_green$mirnaA)
cotargeting_green$mirnaB = gsub("miR-", "", cotargeting_green$mirnaB)

#Create a subtable with only black module cotargeting relationships
cotargeting_green_only = cotargeting_green[cotargeting_green$color_link == "#2CA05A",]


par(mar=c(0,0,0,0))
circos.clear()
circos.par(gap.after = 1)
circos.initialize(sectors = miRNAs_circos$miRNA, xlim = c(0,1))
circos.track(sectors = miRNAs_circos$miRNA, ylim = c(0,1),
             bg.col = miRNAs_circos$color, track.height=0.26)

for (i in 1:nrow(cotargeting_green)){
    circos.link(cotargeting_green$mirnaA[i], 0.5, cotargeting_green$mirnaB[i], 
                0.5, col = cotargeting_green$color_link[i])
    
}

#Replot the green module only cotargeting relationships
for (i in 1:nrow(cotargeting_green_only)){
    circos.link(cotargeting_green_only$mirnaA[i], 0.5, cotargeting_green_only$mirnaB[i], 
                0.5, col = cotargeting_green_only$color_link[i])
    
}

circos.track(track.index = 1, panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                col = "white",
                facing = "clockwise", cex=0.8, niceFacing = TRUE, adj = c(0, 0.5))

}, bg.border = NA)


#-------------------------------------------------------------------------------
# CREATE A CHORD PLOT WITH GREEN MODULE COTARGETING RELATIONSHIPS
#-------------------------------------------------------------------------------

#Create a data frame with frequencies of module combinations
freq_table = as.data.frame(table(cotargeting_green[,c("moduleA", "moduleB")]))

freq_table = freq_table[freq_table$Freq!=0,]

#Sum frequencies of relationships between both modules independent of order

freq_table_formatted = data.frame(moduleA = unique(c(cotargeting_green$moduleA, cotargeting_green$moduleB)),
                                  moduleB = "green",
                                  frequency = NA)

for (i in 1:nrow(freq_table_formatted)){
    modA = freq_table_formatted$moduleA[i]
    modB = freq_table_formatted$moduleB[i]
    
    tab = freq_table[(freq_table$moduleA == modA & freq_table$moduleB == modB) |
                         (freq_table$moduleA == modB & freq_table$moduleB == modA),]
    
    freq_table_formatted$frequency[i] = sum(tab$Freq)
}

#Column with colors for chord plot
freq_table_formatted$color = as.character(freq_table_formatted$moduleA)
freq_table_formatted$color = ifelse(freq_table_formatted$color == "green",
                                    "#2CA05A", freq_table_formatted$color)
freq_table_formatted$color = ifelse(freq_table_formatted$color == "yellow",
                                   "#FFDD55", freq_table_formatted$color)
                                                                       
                                    
#Create chord plot
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
chordDiagram(freq_table_formatted, grid.col = freq_table_formatted$color,
        col=freq_table_formatted$color)

                                    

###############################################################################
#   CHECK IF AVERAGE EXPRESSION LEVEL IS CORRELATED WITH NUMBER OF COTARGETING RELATIONSHIPS
#------------------------------------------------------------------------------

#Calculate mean normalized counts

normalized_counts_global_mean = data.frame(miRNA = rownames(normalized_counts),
                                           mean_expression = rowMeans(normalized_counts),
                                           median_expression = apply(normalized_counts,1,median))


#Merge with table containing number of cotargeting relationships

num_cotarget_relationships$miRNA_renamed2 = gsub("miR-", "",
                                         num_cotarget_relationships$miRNA_renamed)


num_cotarget_relationships = merge(num_cotarget_relationships, 
                                   normalized_counts_global_mean,
                                   by.x ="miRNA_renamed2", by.y ="miRNA", all.x =TRUE)


#Add hub miRNAs
num_cotarget_relationships$hub = "no"
num_cotarget_relationships$hub = ifelse(num_cotarget_relationships$miRNA_renamed2 %in%
                                            KDA_green_drivers, "hub_green",
                                        num_cotarget_relationships$hub)
num_cotarget_relationships$hub = ifelse(num_cotarget_relationships$miRNA_renamed2 %in%
                                            KDA_black_drivers, "hub_black",
                                        num_cotarget_relationships$hub)

pval = round(cor.test(num_cotarget_relationships$Freq, log2(num_cotarget_relationships$mean_expression),
         method = "pearson")$p.value,4)


rcor = round(cor.test(num_cotarget_relationships$Freq, log2(num_cotarget_relationships$mean_expression),
                      method = "pearson")$estimate,2)


p = ggplot(num_cotarget_relationships, 
           aes(x = Freq, y = log2(mean_expression), label = miRNA_renamed2)) +
    geom_point(size=2.5, col = "gray60") + theme_test(base_size=18) +
    theme(axis.text = element_text(color = "black")) +
    labs(x ="No. significant cotargeting relationships", y = "Log2 mean expression")+
    stat_smooth(formula="y~x", method = "lm", color = "gray40") +
    annotate(geom = "text", x = 42.5, y = 23, 
             label = paste0("r = ", rcor, ", p = ", pval), fontface=3, size=6) +
    ylim(0,25) + xlim(0,85) +
    geom_text_repel(data = num_cotarget_relationships[num_cotarget_relationships$hub ==
                                                          "hub_green",],
                    col = "#2CA05A", min.segment.length = 0, size = 5) +
    geom_text_repel(data = num_cotarget_relationships[num_cotarget_relationships$hub ==
                                                          "hub_black",],
                    col = "black", min.segment.length = 0, size = 5)



#------------------------------------------------------------------------------
#Creat a heatmap with the top 20 most expressed miRNAs
normalized_counts_global_mean = normalized_counts_global_mean[
    order(normalized_counts_global_mean$mean_expression, decreasing = TRUE),
]

top_20_expressed_miRNAs = normalized_counts_global_mean$miRNA[1:20]

#Create a heatmap with the top 20 most expressed miRNAs

top20_counts = as.data.frame(normalized_counts[rownames(normalized_counts) %in% top_20_expressed_miRNAs,])

top20_counts = top20_counts[match(top_20_expressed_miRNAs, rownames(top20_counts)),]

colnames(top20_counts) = gsub("O", 0, colnames(top20_counts))


pheatmap::pheatmap(log2(top20_counts), cluster_rows = F,
         color = rev(hcl.colors(n=100, palette = "BurgYl")),
         angle_col = 315)




################################################################################
# ENRICHMENT OF INTRA-MODULAR COTARGETING RELATIONSHIPS
#-------------------------------------------------------------------------------

#The following analysis investigates if there are more co-targeting relationships
#Within a module than between modules to test if miRNAs forming cotargeting networks
#Are also co-expressed
#For this purpose, a normalized ratio of intra-/inter-modular cotargeting relatinships is calculated
#The observed number of co-targeting relationships is normalized by the total numnber of possible cotargeting relationships
#If a module contains n miRNAs, then there are n*(n-1)/2 possible cotargeting relationships within the module
#If there are m miRNAs altogether, then the total possible number of inter-module
#Relationships for a given module is n*(m-m). These are the normalizing factors
#For the observed number of intra- and inter-modular cotargeting relationships
#Only test for green and black module
#---------------------------------
# #First, obtain the total number of miRNAs used for the cotargeting analysis after filtering
# 
filtered_miRNAs = unique(c(all_modules_res$all_pairs$mirnaA, all_modules_res$all_pairs$mirnaB))

#Unlist miRNAs belonging to the same family to individual miRNAs
filtered_miRNAs = unlist(strsplit(filtered_miRNAs, split="_"))

#Obtain module assignment for filtered miRNAs

filtered_miRNAs = miRNAs_cotargeting[miRNAs_cotargeting$miRNA_renamed %in% filtered_miRNAs,]
# 
# 
# #----------------------------------------------------------
# #   BLACK MODULE
# #----------------------------------------------------------
# 
#First calculate the observed number of cotargeting relationships

#Within the module
intramodular_black = nrow(cotargeting_black[cotargeting_black$moduleA == "black" &
                                                cotargeting_black$moduleB == "black",])
#Black module to other modules
intermodular_black = nrow(cotargeting_black[(cotargeting_black$moduleA == "black" &
                                                cotargeting_black$moduleB != "black")|
                                                (cotargeting_black$moduleA != "black" &
                                                     cotargeting_black$moduleB == "black"),])

#Calculate the total possible number of intramodular relationships
n_black = nrow(filtered_miRNAs[filtered_miRNAs$module == "black",])
total_intra_black = n_black*(n_black - 1)/2

#Calculate the total possible number of intermodular black/other relationships
n_other = nrow(filtered_miRNAs[filtered_miRNAs$module != "black",])

total_inter = n_black*n_other

#Calculate the normalized ratio of intra/inter-modular relationships

ratio_normalized_black = (intramodular_black/total_intra_black)/
                            (intermodular_black/total_inter)

#Calculate Odd's ratios and assess association between inter-and intra-modular relationships
tab = rbind(c(intramodular_black, total_intra_black - intramodular_black),
            c(intermodular_black, total_inter - intermodular_black))

fisher_black = fisher.test(tab)
OR_black = round(fisher_black$estimate,3)
CI_black_down = fisher_black$conf.int[1]
CI_black_up = fisher_black$conf.int[2]
p_black = round(fisher_black$p.value,3)
# # ################################################################################
# # #       GREEN MODULE
# # #-------------------------------------------------------------------------------
# # 
#First calculate the observed number of cotargeting relationships

#Within the module
intramodular_green = nrow(cotargeting_green[cotargeting_green$moduleA == "green" &
                                                cotargeting_green$moduleB == "green",])
#green module to other modules
intermodular_green = nrow(cotargeting_green[(cotargeting_green$moduleA == "green" &
                                                 cotargeting_green$moduleB != "green")|
                                                (cotargeting_green$moduleA != "green" &
                                                     cotargeting_green$moduleB == "green"),])

#Calculate the total possible number of intramodular relationships
n_green = nrow(filtered_miRNAs[filtered_miRNAs$module == "green",])
total_intra_green = n_green*(n_green - 1)/2

#Calculate the total possible number of intermodular green/other relationships
n_other = nrow(filtered_miRNAs[filtered_miRNAs$module != "green",])

total_inter = n_green*n_other

#Calculate the normalized ratio of intra/inter-modular relationships

ratio_normalized_green = (intramodular_green/total_intra_green)/
    (intermodular_green/total_inter)


#Calculate Odd's ratios and assess association between inter-and intra-modular relationships
tab = rbind(c(intramodular_green, total_intra_green - intramodular_green),
            c(intermodular_green, total_inter - intermodular_green))

fisher_green = fisher.test(tab)

OR_green = round(fisher_green$estimate,3)
CI_green_down = fisher_green$conf.int[1]
CI_green_up = fisher_green$conf.int[2]
p_green = round(fisher_green$p.value,3)


#Plot results for intra/intermodular analysis
modular_relationships = data.frame(OR = c(OR_black, OR_green),
                                   CI_up = c(CI_black_up, CI_green_up),
                                   CI_down = c(CI_black_down, CI_green_down),
                                   p.val = c(p_black, p_green),
                                   module = c("black module", "green module"))


p = ggplot(modular_relationships, aes(y = OR, x = module, color = module, fill = module)) +
    geom_bar(stat = "identity", alpha = 0.5, show.legend = FALSE, width = 0.6) + 
    geom_errorbar(aes(ymin = CI_down , ymax = CI_up), width = 0.1, show.legend = F) +
    theme_classic(base_size = 19) +
    theme(axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.text.x = element_text(angle=35, hjust = 1)) +
    labs(y = "Odds ratio", x = "") +
    scale_color_manual(values = c("black", "#2CA05A")) +
    scale_fill_manual(values = c("black", "#2CA05A")) +
    ggtitle("Intra/inter-modular \ncotargeting relationships") +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.5)) +
    annotate(geom = "text", x = c(1,2), y = c(modular_relationships$CI_up+0.05),
             label = paste0("p = ", modular_relationships$p.val),
             fontface = 3, size = 4)
    


    
#Version 2    
p + coord_flip() 


###################################################################################

#Get the top genes involved in the biggest number of co-targeting relationships
targets_cotargeting_net = data.frame(table(unlist(strsplit(all_modules_res$bidirectional$common_targets, 
                                                           split = ","))))

targets_cotargeting_net = targets_cotargeting_net[order(targets_cotargeting_net$Freq, decreasing = TRUE),]

#Extract all targets included in the analysis

all_targets = unique(unlist(strsplit(all_modules_res$all_pairs$common_targets, 
                              split = ",")))


#Correlation plot between number of co-targeting relationships and UTR length

targets_cotargeting_net = merge(targets_cotargeting_net,
                                unique(conserved_targets[,c("Gene.Symbol", "UTR_length")]), by.x = "Var1",
                                by.y = "Gene.Symbol", all.x = TRUE)


targets_cotargeting_net = targets_cotargeting_net[order(targets_cotargeting_net$Freq,
                                                        decreasing = TRUE),]

r = round(cor.test(targets_cotargeting_net$Freq, 
         targets_cotargeting_net$UTR_length)$estimate,3)

p=ggplot(targets_cotargeting_net, aes(x = UTR_length, y = Freq, label = Var1)) +
    geom_point(size=2.5, fill = "gray50", color = "gray15", shape=21) +
    stat_smooth(method = "lm", color = "#48C2B4") +
    theme_test(base_size = 19) +
    theme(axis.text = element_text(color ="black")) +
    labs(x = "3' UTR length", y = "No. co-targeting pairs") +
    annotate(geom = "text", x = 20000, y = 1000,
             label = paste0("r = ", r, ", p < 0.0001"), size = 6, fontface="italic") +
    geom_text_repel(data=targets_cotargeting_net[1:10,], min.segment.length = 0,
                    size=5)


################################################################################
#Create a plot with top 10 cotargeting relationships

top10_cotarget = all_modules_res$bidirectional[order(all_modules_res$bidirectional$avg_log_OR,
                                                     decreasing = T),]

top10_cotarget = top10_cotarget[1:10,]

top10_cotarget$pair = paste0(top10_cotarget$mirnaA, " and\n", top10_cotarget$mirnaB)
top10_cotarget$pair = gsub("_", "/", top10_cotarget$pair)
top10_cotarget$pair = factor(top10_cotarget$pair, levels = rev(top10_cotarget$pair))

p = ggplot(top10_cotarget, aes(x = exp(avg_log_OR), y = pair)) +
    geom_bar(stat = "identity", color = "gray15", fill = "#CBCDD9") + 
    theme_test(base_size = 19) +
    theme(axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(y = "", x ="Average odds ratio") +
    ggtitle("Top 10 bidirectional pairs") +
    scale_x_continuous(expand = c(0,0), limits = c(0,10))




#Version 2
p + coord_flip() + theme(axis.text.x = element_text(angle=45, hjust=1))


#As a heatmap?
heatmap_top10 = create_adjacency_matrix(top10_cotarget)

matrix_data = exp(heatmap_top10$adjacency_OR)

rownames(matrix_data) = gsub("_", "/", rownames(matrix_data))
colnames(matrix_data) = gsub("_", "/", colnames(matrix_data))
matrix_data[matrix_data==1] = NA

cols <- hcl.colors(n=100, palette = "BluYl")


pheatmap(matrix_data, cluster_rows = FALSE, cluster_cols = FALSE,
         angle_col = 315, color = rev(cols), na_col = "gray91")