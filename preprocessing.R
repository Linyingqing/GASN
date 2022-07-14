# preprocess the MAF file
library(data.table)
library(stringr)
library(plyr)
library(h2o)
library(fitdistrplus)
library(SuppDists)
preprocessing <- function(M){
  
  #open UVM.mutation.maf.txt
  M <- as.data.frame(fread(file = "data/LAML.csv"))
  
  # mutation type dictionary file from MutSigCV
  dict <- as.data.frame(fread("data/mutation_type_dictionary_file.txt"))
  
  #Gene
  if (("gene" %in% colnames(M))&("Hugo_Symbol" %in% colnames(M))){
    cat("NOTE: Both 'gene' and 'Hugo_Symbol' are present in mutation_file. Using 'gene'. \n")
  }else if("Hugo_Symbol" %in% colnames(M)){
    M$gene <- M$Hugo_Symbol
  }else if("gene" %in% colnames(M)){
    
  }else{
    stop("mutation_file lacks 'gene' or 'Hugo_Symbol' column.")
  }
  
  if("gene" %in% colnames(M)){
    if((str_sub(M$gene[1],1,4) == "ENSG")&("Hugo_Symbol" %in% colnames(M))){
      M$gene <- M$Hugo_Symbol
    }
  }
  
  M <- M[which(M$gene != "TTN"),]
  
  #Patient
  if (("patient" %in% colnames(M))&("Tumor_Sample_Barcode" %in% colnames(M))){
    cat("NOTE: Both 'patient' and 'Tumor_Sample_Barcode' are present in mutation_file. Using 'patient'. \n")
  }else if("patient" %in% colnames(M)){
    #ok
  }else if("Tumor_Sample_Barcode" %in% colnames(M)){
    M$patient <- M$Tumor_Sample_Barcode
  }else{
    stop("mutation_file lacks 'patient' or 'Tumor_Sample_Barcode' column.")
  }
  
  if (length(unique(M$patient)) < 2){
    stop("identification is not applicable to single patients. \n")
  }
  
  #Variant_Classification
  if(!("Variant_Classification" %in% colnames(M))&("type" %in% colnames(M))){
    M$Variant_Classification <- M$type
  }
  if(!("Variant_Classification" %in% colnames(M))){
    stop("mutation_file is missing Variant_Classification")
  }
  
  #Chromosome
  if(!("chr" %in% tolower(colnames(M))) & "chromosome" %in% tolower(colnames(M))){
    colnum <- which(tolower(colnames(M)) %in% "chromosome")
    M$chr <- M[[colnum]]
  }

  if(!("start" %in% tolower(colnames(M))) & "start_position" %in% tolower(colnames(M))){
    colnum <- which(tolower(colnames(M)) %in% "start_position")
    M$start <- M[[colnum]]
  }
  if(!("end" %in% tolower(colnames(M))) & "end_position" %in% tolower(colnames(M))){
    colnum <- which(tolower(colnames(M)) %in% "end_position")
    M$end <- M[[colnum]]
  }
  if(!("ref_allele" %in% colnames(M)) & "reference_allele" %in% tolower(colnames(M))){
    colnum <- which(tolower(colnames(M)) %in% "reference_allele")
    M$ref_allele <- M[[colnum]]
  }
  if(!("newbase" %in% tolower(colnames(M)))){
    if("tumor_seq_allele1" %in% tolower(colnames(M))){
      colnum <- which(tolower(colnames(M)) %in% "tumor_seq_allele1")
      M$newbase <- M[[colnum]]
      if("tumor_seq_allele2" %in% tolower(colnames(M))){
        colnum.allele1 <- which(tolower(colnames(M)) %in% "tumor_seq_allele1")
        idx <- which(M$ref_allele == M[[colnum.allele1]])
        colnum.allele2 <- which(tolower(colnames(M)) %in% "tumor_seq_allele2")
        M$newbase[idx] <- M[[colnum.allele2]][idx]
      }
    }
  }

  
  if(length(which(!(M$Variant_Classification %in% unique(dict$Variant_Classification))))>0){
    M <- M[-which(!(M$Variant_Classification %in% unique(dict$Variant_Classification))),]
  }
  
  #Effect
  if ("effect" %in% colnames(M)){
    # cat("Will use the pre-existing effect column.\n")
    M$effect <- gsub("^flank.*","noncoding",M$effect)
    M_effectnames <- c("noncoding","silent","nonsilent","null")
    if (!any(unique(M$effect) %in% M_effectnames)){
      stop("in mutation_file, 'effect' must be one of
         noncoding/silent/nonsilent/null")
    }
  }else{
    flag_num <- match(toupper(M$Variant_Classification),
                      toupper(dict$Variant_Classification),nomatch =nrow(dict)+1)
    dict <- rbind(dict,dict[nrow(dict),])
    dict[nrow(dict)+1,] <- "unknown"
    
    M$effect <- dict$effect[flag_num]
    bad <- which(M$effect == "unknown")
    if(length(bad)>0){
      cat(sprintf("WARNING: %d/%d mutations could not be mapped to
                  effect using mutation_type_dictionary_file:\n",
                  length(bad),length(M$effect)))
      cat("They will be removed from the analysis. \n")
      M <- M[-bad]
    }
    if(nrow(M) == 0){
      stop("No mutations left!")
    }
  }
  
  M <- subset(M,select = c("gene","chr","start","end","Variant_Classification","ref_allele","newbase","patient","effect"))
  return(M)
}


