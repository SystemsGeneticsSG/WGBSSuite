#sink(file="/dev/null")
bsmooth_for_bench <- function(data,number_of_samples,number_of_replicas,biseq_min_sites,biseq_max_dist){
    
    library("bsseq")
    reformatted_for_bsmooth <- reformat_for_bsmooth(data,number_of_samples,number_of_replicas)
    return(reformatted_for_bsmooth)
}


reformat_for_bsmooth <- function(d,number_of_samples,number_of_replicas){
	total <- number_of_samples*number_of_replicas
	dm <- dim(d)
	pos<- d[,1]
	states <- d[,2]
	l<-dm[1];
	chr <- matrix(data=1,nrow=l,ncol=1)
	BS=list() 
	groups <- list()
	pos_mat <- matrix(data=NA,nrow=l,ncol=number_of_samples)
	for(i in 1:number_of_samples){

    	for(j in 1:number_of_replicas){
    		index <- ((i-1)*number_of_replicas)+j
    		col <- (((i-1)*(number_of_replicas*4))+((j-1)*4)+1)+2
    		sample_id<-paste0(i,j, collapse = NULL)	
    		rep<-paste0("replicate",j, collapse = NULL)	
   			BS[[index]] <- BSseq(pos = as.numeric(as.character(pos)), chr = chr, M = as.matrix(as.numeric(as.character(d[,col])), ncol = 1),Cov = as.matrix(as.numeric(as.character(d[,col+2])), ncol = 1), sampleNames = paste("for",sample_id,sep=""))
			sampleNames(BS[[index]]) <- sample_id
			if(index == 1){
    			combo <- BS[[index]]
    			pData(combo)$Rep <- matrix(data=rep,nrow=1,ncol=total)	
    		}else{
    			combo <- combine(combo,BS[[index]])
    			pData(combo)$Rep[index] <-rep
    		}
    		if(j==1){
    			groups[[i]]<-matrix(data=sample_id,nrow=1,ncol=number_of_replicas)	
    		}else{
    			groups[[i]][j] <- sample_id
    		}
			
			
    	}
	}

			validObject(combo)
			BS.combo.smooth <-BSmooth(combo)
            #reset the timer
            ptm <- proc.time()
			BS.fish <- calculate_fisher(combo,groups[[1]],groups[[2]])
            fish_time<-(proc.time() - ptm)[1]
            #reset the timer
            ptm <- proc.time()
        noLocal <- FALSE

        # sometimes local.correct fails and the whole simulation crashes. I've put this catch in to recover from this for now - TJH
        BS.combo.smooth.tstat <- try(BSmooth.tstat(BS.combo.smooth,group1 = groups[[1]],group2 = groups[[2]],estimate.var = "same",local.correct = TRUE,verbose = TRUE))
        if(class(BS.combo.smooth.tstat) == "try-error") {
            BS.combo.smooth.tstat <- BSmooth.tstat(BS.combo.smooth,group1 = groups[[1]],group2 = groups[[2]],estimate.var = "same",local.correct = FALSE,verbose = TRUE)
            warning("Local correction in BSmooth failed.")
            noLocal <- TRUE
        }
        
            b_time<-(proc.time() - ptm)[1]
			return(list(BS.combo.smooth.tstat,BS.fish,combo,BS.combo.smooth,fish_time,b_time))
}

calculate_fisher <- function(bsmooth_object,g1,g2){
	results <- fisherTests(bsmooth_object, group1 = g1, group2 = g2,returnLookup = TRUE)
	return(results)

}




bsmooth_benchmark <- function(BS.chr1.smooth.tstat,bsmooth_thresh,dfs,locs,min_number_of_CpGs){
	l<-length(dfs)
	for(j in 1:length(bsmooth_thresh)){
			dmrs.chr1 <- dmrFinder(BS.chr1.smooth.tstat, qcutoff = c(0+bsmooth_thresh[j], 1-bsmooth_thresh[j]),column=c("tstat"))
			dmrs.chr1.thresh <- subset(dmrs.chr1, n >= min_number_of_CpGs & abs(meanDiff) >= 0.01)
			diff_meth <- matrix(data=0,nrow=1,ncol=l)
			if(length(dmrs.chr1.thresh[,2])>0){
				for(i in 1:length(dmrs.chr1.thresh[,2])){
					if((dim(dmrs.chr1.thresh)[1])>0){
						diff_meth[match(dmrs.chr1.thresh[i,2],locs):match(dmrs.chr1.thresh[i,3],locs)]<-1
					}
				}
			}
			x<-score_overlap(dfs,diff_meth,bsmooth_thresh[j])
		}
}

par_bsmooth_benchmark <- function(BS.chr1.smooth.tstat,bsmooth_thresh,dfs,locs,bsmooth_second_cutoff,min_number_of_CpGs){
        l<-length(dfs)
        foreach(j=1:length(bsmooth_thresh),.combine=cbind) %dopar% {
                        dmrs.chr1 <- dmrFinder(BS.chr1.smooth.tstat, qcutoff = c(0+bsmooth_thresh[j], 1-bsmooth_thresh[j]),stat="tstat")
                        dmrs.chr1.thresh <- subset(dmrs.chr1, n >= min_number_of_CpGs & abs(meanDiff) >= bsmooth_second_cutoff)
                        diff_meth <- matrix(data=0,nrow=1,ncol=l)
                        if(length(dmrs.chr1.thresh[,2])>0){
                                for(i in 1:length(dmrs.chr1.thresh[,2])){
                                        if((dim(dmrs.chr1.thresh)[1])>0){
                                                diff_meth[match(dmrs.chr1.thresh[i,2],locs):match(dmrs.chr1.thresh[i,3],locs)]<-1
                                        }
                                }
                        }
                        x<-score_overlap(dfs,diff_meth,bsmooth_thresh[j])
                        print(x)
                }
}

par_bsmooth_benchmark_DMRs <- function(BS.chr1.smooth.tstat,bsmooth_thresh,dfs_locs,locs,bsmooth_second_cutoff,min_number_of_CpGs){
        l<-length(locs)
        dmr_locs<-list()
        dmr_locs_ind<-1
        foreach(j=1:length(bsmooth_thresh),.combine=cbind) %dopar% {
                        dmrs.chr1 <- dmrFinder(BS.chr1.smooth.tstat, qcutoff = c(0+bsmooth_thresh[j], 1-bsmooth_thresh[j])	,stat="tstat")
                        dmrs.chr1.thresh <- subset(dmrs.chr1, n >= min_number_of_CpGs & abs(meanDiff) >= bsmooth_second_cutoff)
                        diff_meth <- matrix(data=0,nrow=1,ncol=l)
                        if(length(dmrs.chr1.thresh[,2])>0){
                                for(i in 1:length(dmrs.chr1.thresh[,2])){
                                        if((dim(dmrs.chr1.thresh)[1])>0){
                                                dmr_locs[[dmr_locs_ind]]<-c(dmrs.chr1.thresh[i,2],dmrs.chr1.thresh[i,3])
                                                dmr_locs_ind<-dmr_locs_ind+1
                                        }
                                }
                        }

                        x<-score_overlap_DMR(dfs_locs,dmr_locs,bsmooth_thresh[j])
                        print(x)
                }
}

par_fisher_benchmark <- function(fish_results,fish_threshs,dfs,locs,min_number_of_CpGs){
	l<-length(dfs)
	foreach(j=1:length(fish_threshs),.combine=cbind) %dopar% {
			diff_meth <- matrix(data=0,nrow=1,ncol=l)
				for(i in 1:length(fish_results)){
					if(fish_results[i]<=fish_threshs[j]){
						diff_meth[i]<-1
					}
				}
			x<-score_overlap(dfs,diff_meth,fish_threshs[j])
			print(x)
		}
	}
