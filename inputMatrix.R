# get the input matrix for FI-net
inputMatrix <- function(M){
  
  #save the number of deleterious mutation effect of each gene
  M$deleteriousMut <- 0
  M$deleteriousMut[which(M$effect %in% c("null","nonsilent"))] <- 1
  deleterious.mut <- as.matrix(tapply(M$deleteriousMut, M$gene, sum))
  
  #take the mutation number of each gene as an explaination variable
  M$allMut <- 1
  all.mut <- as.matrix(tapply(M$allMut,M$gene,sum))
  
  # take the ratio of null and nonsilent mutation as an explaination variable
  deleterious.rate <- deleterious.mut/all.mut
  
  #tabke the null and nonsilent mutation ratio and mutation number of genes
  gene.effect <- data.frame(gene=rownames(deleterious.rate), totalMutNum=all.mut,deleteriousMutNum=deleterious.mut)
  
  genes <- unique(M$gene)
  patients <- unique(M$patient)
  gene.pat.matrix <- as.data.frame(matrix(data = 0,nrow = length(genes),ncol = length(patients)))
  dimnames(gene.pat.matrix) <- list(genes,patients)
  for (j in 1:length(patients)) {
    pat.name <- colnames(gene.pat.matrix)[j]
    pat.FIS <- M[which(M$patient %in% pat.name),]
    gene.FIS <- as.data.frame(tapply(pat.FIS$FIscore, pat.FIS$gene, sum))
    flag.gene <- match(rownames(gene.pat.matrix),rownames(gene.FIS))
    flag.gene.na <- which(!is.na(flag.gene))
    gene.pat.matrix[flag.gene.na,j] <- gene.FIS[flag.gene[flag.gene.na],1]
  }
  
  sdFIS <- as.data.frame(apply(gene.pat.matrix, 1, sd))
  sdFIS$gene <- rownames(sdFIS)
  colnames(sdFIS) <- c("sd","gene")
  
  gene.score <- ddply(M,"gene",summarise,sum(FIscore))
  colnames(gene.score) <- c("gene","mutationScore")
  
  gene.feature <- as.data.frame(fread("data/geneFeature.txt"))
  
  colnames(gene.feature)[1] <-"gene"
  input.matrix <- merge(gene.score,gene.feature,by="gene")
  input.matrix <- merge(input.matrix,gene.effect,by="gene")
  input.matrix <- merge(input.matrix,sdFIS,by="gene")
  #write.table (input.matrix, file ="input.matrix_BCRA.txt",row.names = FALSE, col.names =TRUE, quote =FALSE)
  
  return(input.matrix)
}
