library(methylSig)

par_methsig_benchmark <- function(methsig_cutoffs,methsig_qvals,locs,dfs,max_distance,min_CpG,methdiff,methdiff_cutoff){
  
  #methsig_cutoffs<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
  #methsig_ROC <- matrix(data=0,nrow=2,ncol=length(methsig_cutoffs))
  foreach(i=1:length(methsig_cutoffs),.combine=cbind) %dopar% {
    #methsig_points<-extract_predicted_methsig_cDMR_points(methsig_qvals,methsig_cutoffs[i],8,locs,10)
    methsig_dfs<-extract_blocks_methsig(methsig_qvals,methsig_cutoffs[i],min_CpG,locs,max_distance,methdiff,methdiff_cutoff)#
    x<-score_overlap(dfs,methsig_dfs[[2]],methsig_cutoffs[i])
    print(x)
  }
}

par_methsig_benchmark_DMRs <- function(methsig_cutoffs,methsig_qvals,locs,dfs_locs,max_distance,min_CpG,methdiff,methdiff_cutoff){
  

  foreach(i=1:length(methsig_cutoffs),.combine=cbind) %dopar% {
  #for(i in 1:length(methsig_cutoffs)) {
    #methsig_points<-extract_predicted_methsig_cDMR_points(methsig_qvals,methsig_cutoffs[i],8,locs,10)
    methsig_dfs<-extract_blocks_methsig(methsig_qvals,methsig_cutoffs[i],min_CpG,locs,max_distance,methdiff,methdiff_cutoff)#
    x<-score_overlap_DMR(dfs_locs, methsig_dfs[[4]],methsig_cutoffs[i])
      print(x)
  }
}


methsig_benchmark <- function(methsig_cutoffs,methsig_qvals,locs,dfs){
  
  methsig_cutoffs<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
  #methsig_ROC <- matrix(data=0,nrow=2,ncol=length(methsig_cutoffs))
  s<-list()
  for(i in 1:length(methsig_cutoffs)){
    #methsig_points<-extract_predicted_methsig_cDMR_points(methsig_qvals,methsig_cutoffs[i],8,locs,10)
    methsig_dfs<-extract_predicted_methsig_cDMR(methsig_qvals,methsig_cutoffs[i],1,locs,10)#
    x<-score_overlap(dfs,methsig_dfs)
    s[[length(s)+1]] <- x
  }
  return(s)
}

methsig_for_bench <- function(data_set,number_of_samples,number_of_replicas,locs){
	#load a methsig object from file
	load("meth.rda")
	reformatted_for_methsig <- reformat_for_methsig(data_set,number_of_samples,number_of_replicas,locs,meth)

	return(reformatted_for_methsig)
}


reformat_for_methsig <- function (data_set,number_of_samples,number_of_replicas,locs,meth){
  data_set<-t(data_set)
  l<-length(data_set[[1]][,1])
  size_of_set <- number_of_samples*number_of_replicas
  cov <- matrix(data=0,nrow=l,ncol=size_of_set)
  numT <- matrix(data=0,nrow=l,ncol=size_of_set)
  numC <- matrix(data=0,nrow=l,ncol=size_of_set)

  ids <-  data_set[[2]]
  start <- data_set[[2]]
  end <-data_set[[2]]+1
  chr <- matrix(data=1,nrow=l,ncol=1)
  strand <- matrix(data=F,nrow=l,ncol=1)
  sample.ids <-  matrix(data='A',nrow=size_of_set,ncol=1)
  sample.filenames <- matrix(data=0,nrow=6,ncol=1)
  treatment <- matrix(data=1,nrow=size_of_set,ncol=1)
  #destranded
  destranded <- TRUE
  #resolution
  resolution <- "base"
  #options
  options <- "maxCount=500 & minCount=1 & assembly=hg18 & context=CpG"
  
  for(i in 1:number_of_samples){

  	for(j in 1:number_of_replicas){
  		row <- (((i-1)*(number_of_replicas*4))+((j-1)*4)+1)+2
  		index <- ((i-1)*number_of_replicas)+j
  		cov[,index]<-data_set[[1]][,(row+2)]
  		numT[,index]<-data_set[[1]][,row+1]-data_set[[1]][,row]
  		numC[,index]<-data_set[[1]][,row]
  		treatment[index]<-i
  		sample.ids[index]<-paste0(i,j)
  	}
  }
  meth@data.numCs <- numC
  meth@data.coverage <- cov
  meth@data.numTs <- numT
  meth@data.ids <- ids
  meth@data.start <- start
  meth@data.end <- end
  meth@treatment <- as.numeric(treatment)
  meth@sample.ids <- as.character(sample.ids)
  meth@data.chr <- meth@data.chr[1:l]
  meth@data.strand <- meth@data.strand[1:l]
  gs<-seq(1,number_of_samples)
  myDiffSigboth = methylSigCalc(meth, groups=gs, min.per.group=number_of_samples,local.disp=TRUE, winsize.disp=2000,num.cores=12)

  return(myDiffSigboth)
}

extract_blocks_methsig <-function (prob_diff,thresh,cutoff,locs,max_dist,methdiff,methdiff_cutoff){
  preds<-list()
  preds_ind<-1
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
          preds[[preds_ind]]<-c(on,off)
          preds_ind<-preds_ind+1
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
  return(list(points,test,cutoff,preds))
}

extract_predicted_methsig_cDMR <-function (prob_diff,thresh,cutoff,locs,max){
  diff_meth <- matrix(data=0,nrow=1,ncol=length(prob_diff))
  c<-1
  
  for(i in 2:length(prob_diff)){
    
    if((prob_diff[i] < thresh)&((locs[i]-locs[i-1]) < max)){
      diff_meth[i] = c;
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
  for(i in 1:length(diff_meth)){
    if(diff_meth[i] == 1){
      if(on != 0){
        off <- i
      }else{
        on <- i
      }
    }else{
      if(on != 0){
        if(off != 0){
          diff_meth[on:off]<-1
          on<-0
          off<-0
        }
      }else{
        on<-0
        off<-0
      }
    }
  }
  return(diff_meth)
}