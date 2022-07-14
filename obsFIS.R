# calculate the observed functional impact score for each mutation in MAF file based on MutationAccessor

obsFIS <- function(M){

  M$mutation <- paste("hg19",M$chr,M$start,M$ref_allele,M$newbase,sep = ",")
  
  MA.scores.direct <- "MA_scores_rel3_hg19_full"
  #????ȡ?ļ???MA_scores_rel3_hg19_full?¸????ļ?
  MA.scores.filename <- list.files(MA.scores.direct)
  # list all the files of mutation accessor
  
  M.score <- as.data.frame(matrix(data = NA,ncol = ncol(M)+2))
  colnames(M.score) <- c(colnames(M),"mutation","score")
  M.score <- M.score[-1,]
  # build a data.frame M.score to store the mutaiton score
  
  MA.score.effect <- as.data.frame(matrix(data = NA,ncol = ncol(M)+2))
  MA.score.effect <- MA.score.effect[-1,]
  # build a data.frame MA.score.effect to store the mutation score of MutationAccessor
  
  for(i in c(1:22,"X","Y")){
    #print(i)
    # take the chr i to deal with the mutation score for every mutation
    M.chr.i <- M[which(M$chr == i),]
    # take the mutation score of chromosome i
    
    if(nrow(M.chr.i) > 0){
      # read the mutation data of MutationAccessor
      MA.scores.filename.i <- paste(MA.scores.direct,"/","MA_scores_rel3_hg19_chr",i,"_full.csv",sep = "")
      MA.scores <- fread(file = MA.scores.filename.i)
      
      # assign mutation score for missense mutation using the MutationAccessor FI score
      mutation.intersect <- intersect(M.chr.i$mutation,MA.scores$Mutation)
      M.filter <- M.chr.i[M.chr.i$mutation %in% mutation.intersect,]
      MA.scores.filter <- MA.scores[MA.scores$Mutation %in% mutation.intersect,]
      flag.match <- match(M.filter$mutation,MA.scores.filter$Mutation)
      M.filter$score <- MA.scores.filter$`FI score`[flag.match]
      
      # take the effective  mutation score of MutationAccessor, and store in MA.score.effect
      MA.score.effect <- rbind(MA.score.effect,M.filter)
      
    }
  }
  
  if(nrow(MA.score.effect) < nrow(M)*0.5){
    bad <- nrow(M) - nrow(MA.score.effect)
    cat(sprintf("WARNING: %d/%d mutations could not be mapped to
                mutation functional impact score using mutation_accessor_file:\n",
                bad,nrow(M)))
    #stop("The maf file does not match with the mutation impact file of hg19")
  }
  
  # remove the mutation data of NA in MA.score.effect
  effect.score <- MA.score.effect[which(!is.na(MA.score.effect$score)),]
  
  # match the mutation of M with the mutation of effect.score
  flag_num <- match(toupper(M$mutation),
                    toupper(effect.score$mutation),nomatch = nrow(effect.score)+1)
  
  # assign mutation FI score to M
  MA.score <- c(effect.score$score,NA)
  M$FIscore <- MA.score[flag_num]
  
  # calculate the mean mutation FI score of every mutation effect
  M.ex.NA <- M[which(!is.na(M$FIscore)),]
  mut.effect.mean <- tapply(M.ex.NA[,"FIscore"], M.ex.NA[,"effect"], mean)
  
  mut.effect.score <- data.frame(effect=c("noncoding","nonsilent","null","silent"),
                               score=c(1,2,3,0))
  flag_effect <- mut.effect.score$effect %in% names(mut.effect.mean)
  
  # assign the mutation FI score of effect which is not in the MutationAccessor file
  if(sum(flag_effect) < 4){
    mut.effect.mean <- data.frame(effect=c(names(mut.effect.mean),
                                         as.vector(mut.effect.score$effect)[!flag_effect]),
                                score=c(mut.effect.mean,as.vector(mut.effect.score$score)[!flag_effect]))
  }else{
    mut.effect.mean <- data.frame(effect=names(mut.effect.mean),score=mut.effect.mean)
  }
  
  # assign the mutation FI score of NA in mutationdata by the mean of mutation effect
  M_NA <- M[which(is.na(M$FIscore)),]
  flag.match.effect <- match(M_NA$effect,mut.effect.mean$effect)
  M_NA$FIscore <- mut.effect.mean$score[flag.match.effect]
  
  # rbind the mutation score with NA and without NA
  M.new <- rbind(M.ex.NA,M_NA)
  M<-M.new
  #write.table (M.new, file ="Mnew.txt",row.names = FALSE, col.names =TRUE, quote =FALSE)
  
  return(M.new)
}


