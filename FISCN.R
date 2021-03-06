# FI-net main code
FInet <- function(input.matrix,N){
  
  compartment=read.csv('data/Gen_compartment_together.csv')
  seperate = compartment
  com=unique(seperate[,2])
  
  X <- subset(input.matrix,select = -c(mutationScore,gene))
  X <- scale(X)
  Y <- input.matrix$mutationScore
  gene.feature.score <- as.data.frame(cbind(Y,X))
  colnames(gene.feature.score)[1] <- "mutationScore"
  rownames(gene.feature.score) <- input.matrix$gene
  
  M.cluster <- gene.feature.score
  M.cluster[,1:(ncol(gene.feature.score)-1)] <- input.matrix[,2:ncol(input.matrix)]
  
  group=list()

  for(i in 1:length(com))#####divide matrix into 11 compartment subgroups将矩阵划分为11个隔室子群
  {
    com_i=seperate[which(seperate[,2]==com[i]),1]
    group_i =  M.cluster[which(rownames( M.cluster)%in%com_i),]
    group[[i]]=cbind(group_i,compartment=as.character(com[i]))
    
  }
  
  M.total <- as.data.frame(matrix(data = NA,nrow = 0,ncol = ncol(M.cluster)))
  M.total_finall = list()
  for(i in 1:length(group))###calculate the entrophy for each sub-group计算每个子组的熵
  {
     
    M.cluster.i <- group[[i]]
    estimate.FIS <- M.cluster.i$estimateFI
    obs.FIS = M.cluster.i$mutationScore
    flag.trun <- which(estimate.FIS < quantile(estimate.FIS,0.05))
    estimate.FIS.quantile <- estimate.FIS[-flag.trun]
    
    if(min(estimate.FIS.quantile) <= 0){
      min.FIS <- min(estimate.FIS.quantile)
      estimate.FIS.quantile <- estimate.FIS.quantile - min.FIS + 0.05
      obs.FIS <- obs.FIS - min.FIS + 0.05
    }
    

   
    FIS.distribution1 <- fitdist(estimate.FIS.quantile,distr ="gamma",method = "mle")

    
    
    for (j in 1:nrow(M.cluster.i)) {
      if(obs.FIS[j] <= 0){
        M.cluster.i$pgamma[j] <- 1
      }else{
        
        ks.result1 <- ks.test(obs.FIS[j],"pgamma",FIS.distribution1$estimate[1],FIS.distribution1$estimate[2],alternative = "less")
        M.cluster.i$pgamma[j] <- ks.result1$p.value
      }
    }
    
    M.cluster.i$qgamma <- stats::p.adjust(M.cluster.i$pgamma,method = "fdr",nrow(M.cluster.i))

    
    M.cluster.i$qValue = M.cluster.i$qgamma
    
    M.cluster.i = cbind(M.cluster.i,gene = rownames(M.cluster.i))
    M.total_finall = rbind(M.cluster.i,M.total_finall)
    
  }
  
  M.total <- data.frame(gene = M.total_finall $gene,qValue = M.total_finall$qValue)

  M.total=M.total[order(as.numeric(M.total[,2]),decreasing=F),]
  
  #M.total = unique(M.total[,])
  M.total1 = M.total[duplicated(M.total$gene)==F,]
  
  M.total2 = M.total1$gene[which(M.total1$qValue<=0.05)]
  
  return(M.total)
} 

