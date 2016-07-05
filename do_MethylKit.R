  library(methylKit)
  

par_methylkit_benchmark <- function(methsig_cutoffs,methsig_qvals,locs,dfs,max_distance,min_CpG,methdiff,methdiff_cutoff){

  foreach(i=1:length(methsig_cutoffs),.combine=cbind) %dopar% {
    methsig_dfs<-extract_blocks_methsig(methsig_qvals,methsig_cutoffs[i],min_CpG,locs,max_distance,methdiff,methdiff_cutoff)#
    x<-score_overlap(dfs,methsig_dfs[[2]],methsig_cutoffs[i])
    print(x)
  }
}

methylkit_for_bench <- function(data,number_of_samples,number_of_replicas,output_path){
  reformatted_for_methylkit <- reformat_for_methylkit(data,number_of_samples,number_of_replicas,output_path)
  return(reformatted_for_methylkit)
}

reformat_for_methylkit <- function(d,number_of_samples,number_of_replicas,output_path){


  total <- number_of_samples*number_of_replicas
  dm <- dim(d[[1]])
  # I can't quite remember when the cast to integer is necessary, but it doesn't do any harm - TJH
  pos<- as.integer(d[[1]][,1])
  states <- d[[1]][,2]
  l<-dm[1];
  results=list() 
  sample_type <- matrix(data=NA,nrow=total,ncol=1)
  sample_ids <- matrix(data=NA,nrow=total,ncol=1)
  sample_chrs <- matrix(data='chr1',nrow=l,ncol=1)
  sample_strand <- matrix(data='+',nrow=l,ncol=1)
  filelist <- list()
  treatment <-matrix(data=0,nrow=total,ncol=1)
  sampleid <- list()
  data_store <- list() 
  for(i in 1:number_of_samples){
      data_set <- list()      
      for(j in 1:number_of_replicas){
        col <- (((i-1)*(number_of_replicas*4))+((j-1)*4)+1)+2
        index <- ((i-1)*number_of_replicas)+j
        data_for_file <- matrix(data=NA,nrow=l,ncol=7)
        colnames(data_for_file)<-c('chrBase','chr','base','strand','coverage','freqC','freqT') 
        data_for_file[,1]<-pos
        data_for_file[,2]<-sample_chrs
        data_for_file[,3]<-pos
        data_for_file[,4]<-sample_strand
        data_for_file[,5]<-d[[1]][,(col+1)]
        print(col)
        x <- d[[1]][,col]
        print(x[1:10])
        n <- d[[1]][,(col+1)]
        print(n[1:10])
        freqC <- (x/n)*100
        freqT <- 100-freqC
        data_for_file[,6]<-freqC
        data_for_file[,7]<-freqT
        mfile<-paste0(output_path,"_","methylkit_",index)
        filelist<-c(filelist,mfile)
        sampleid[[index]]<-paste0("bench",index)
        treatment[index,1]<-i-1
        write.table(data_for_file, file = mfile, append = FALSE, quote = FALSE, sep = "\t",
          eol = "\n", na = "NA", dec = ".", row.names = FALSE,
          col.names = FALSE)
      }

     
  }
  # needs header = FALSE or the first cytosine goes missing - TJH
       myobj=read(filelist,sample.id=sampleid,assembly="bench",treatment=as.numeric(treatment),context="CpG", header = FALSE)
       meth=unite(myobj, destrand=FALSE)
       myDiff=calculateDiffMeth(meth)

  return(myDiff)
}


extract_blocks_methylKit <-function (prob_diff,thresh,cutoff,locs,max_dist,methdiff,methdiff_cutoff){
  
  diff_meth <- matrix(data=0,nrow=1,ncol=length(prob_diff))
  test <- matrix(data=0,nrow=1,ncol=length(prob_diff))
  meth_sites <- matrix(data=0,nrow=1,ncol=length(prob_diff))
  c<-1
  for(i in 2:length(prob_diff)){
 
    if((prob_diff[i] < thresh)&&(abs(methdiff[i]) > methdiff_cutoff)){
      diff_meth[i] = c;
      test[i] = 1
      c<-c+1
    }else{
      c<-1
    }
  }
  
  grace<-0
  for(i in length(diff_meth):1){
    if(grace > 0){
      diff_meth[i] = 1;
      grace<-grace-1;
    }else if(diff_meth[i] > cutoff){
      grace<-diff_meth[i]-1
      diff_meth[i] = 1;
    }else{
      diff_meth[i] <- 0;
      grace<-0
    }
  }
  
  on<-0
  off<-0
  points=list() 
  c<-0
  for(i in 1:length(diff_meth)){

    if(diff_meth[i] == 1){
      if(on != 0){
        off <- i
        c<-c+1
      }else{
        on <- i
        c<-c+1
      }
    }else{
      if(on != 0){
        if(off != 0){
          points[[(length(points)+1)]]<-c(on,off,c)
          meth_sites[on:off]<-1
          on<-0
          off<-0
          c<-0
        }
      }else{
        on<-0
        off<-0
        c<-0
      }
    }
  }
  if(on != 0){
      if(off != 0){
        points[[(length(points)+1)]]<-c(on,off,c)
        meth_sites[on:off]<-1
      }
  }
  return(list(points,test,cutoff))
}




  fdr_biSeq <- function(locCor,rej_fdr,clus_fdr){
  clusters.rej <- testClusters(locCor,FDR.cluster = rej_fdr)
  
  clusters.trimmed <- trimClusters(clusters.rej,FDR.loc = clus_fdr)
  DMRs <- findDMRs(clusters.trimmed,
                   max.dist = 1000,
                   diff.dir = TRUE)
  return(DMRs)
}

score_biSeq <- function(DMRs,locs){
  l <- length(locs)
  diff_meth <- matrix(data=0,nrow=1,ncol=l)
  for(i in 1:dim(DMRs)[1]){
    start<-DMRs[i,1]
    end<-DMRs[i,2]
    start_i<-which(locs==start)
    end_i<-which(locs==end)
    diff_meth[start_i:end_i]<-1
  }
  return(diff_meth)
}
