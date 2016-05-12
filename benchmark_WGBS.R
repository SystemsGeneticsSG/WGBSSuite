#load the required libraries and functions
sink(file="/dev/null")
source("do_MethSig.R")
source("do_BSmooth.R")
source("do_MethylKit.R")
library(doMC)
sink()

###############################################################################
#reads a file formatted for the benchmarking function, this is normally the   #
#file created by the simulation script but this function could be called from #
#elsewhere with a userdefined filename                                                                            #
###############################################################################
load_data_from_file <- function(filename){
        data <- read.table(filename,header = T)
        locs <- data[,1]
        states <- data[,2]
        return(list(data,locs,states))
}

###############################################################################
#read a variable from the command line                                        #
###############################################################################
get_variable <- function (message,default){
#print the prompt to the command line
cat(message)
#capture the user input from stdin
variable <- as.numeric(readLines(con="stdin", 1))
#if the user doesn't anything then use the default value
if(is.na(variable)){
  variable <- default
}
return(variable)
}

###############################################################################
#read a text from the command line                                        #
###############################################################################
get_text <- function (message,default){
#print the prompt to the command line
cat(message)
#capture the user input from stdin
variable <- readLines(con="stdin", 1)
#if the user doesn't anything then use the default value
if(variable == ''){
variable <- default
}
return(variable)
}

###############################################################################
#interactively initiate the function                                          #
###############################################################################
init_interactive_benchmark <- function (){

#get the variables from the command line
filename <-get_text("Please enter the filename (default: simulated_WGBS_data.txt)","/tmp/simulated_WGBS_data.txt")
analysis_folder <-get_text("Please enter the analysis folder you wish to store the results (default: /tmp/)","/tmp/")
analysis_name <-get_text("Please enter a name you wish to store the resultsunder (default: example_WGBS_run)","example_WGBS_run")
output_path <- paste0(analysis_folder,analysis_name)
number_of_samples<-2
number_of_replicas <-get_variable("Please enter number_of_replicas (default: 3)",3)
min_number_CpGs <-get_variable("Please enter the min_number_CpGs (default: 1)",1)
max_distance_between_DMRs <-get_variable("Please enter the max_distance_between_DMRs (default: 9999999)",9999999)
lowest_order <-get_variable("Please enter the lowest_order (default: 15)",15)
breaks_per_order <- get_variable("Please enter the breaks_per_order (default: 10)",10)
FDR_restarts <-100
FDR_iterations <-10000
biseq_fdrs2 <-0
FDR_cutoffs <- sort(c(seq(1,9,10/breaks_per_order) %o% 10^-(1:lowest_order)))
Methsig_cutoffs <- seq(from=0,to=1,length.out=length(FDR_cutoffs))
biseq_min_sites <-0
biseq_max_dist <-0
bsmooth_second_cutoff <-get_variable("Please enter the bsmooth_second_cutoff (default: 0.1)",0.1)
methdiff_cutoff <-get_variable("Please enter the methylsig_cutoff (default: 5)",5)
methykit_cutoff <-get_variable("Please enter the methdiff_cutoff (default: 5)",5)
type <-get_text("Do you want to use the binomial or truncated simulator (default: binomial)","binomial")

#run the benchmark
rd<-auto_init_benchmark(filename,number_of_samples,number_of_replicas,min_number_CpGs,max_distance_between_DMRs,0,1,lowest_order,breaks_per_order,FDR_restarts,FDR_iterations,biseq_fdrs2,biseq_min_sites,biseq_max_dist,bsmooth_second_cutoff,methdiff_cutoff,output_path,methykit_cutoff,type)
return(rd)

}


###############################################################################
#run main function                                          #
###############################################################################
auto_init_benchmark <- function(filename,number_of_samples,number_of_replicas,min_number_CpGs,max_distance_between_DMRs,INLA_prior_mean,INLA_prior_variance,lowest_order,breaks_per_order,FDR_restarts,FDR_iterations,biseq_fdrs2,biseq_min_sites,biseq_max_dist,bsmooth_second_cutoff,methdiff_cutoff,output_path,methykit_cutoff,type){
        sink(file="/dev/null")
        #create an array of FDR values that are equally spaced accoress each order of magnitiude (ie for a log scale)
        FDR_cutoffs <- sort(c(seq(1,9,10/breaks_per_order) %o% 10^-(1:lowest_order)))

        names<-list()
        ROC<-list()
        ts<-list()
        ind<-1
        #for the methods that use a q-value take the same number of increments and place them between 0 and 1 for calculating the ROC curves
        Methsig_cutoffs <- seq(from=0.001,to=0.999,length.out=length(FDR_cutoffs))

        #set up the multithreading, hardcoded to 16 threads at the moment
        registerDoMC(16)

        #load the data from file
        dataset <- load_data_from_file(filename)
        #extract the differentially methylated bases from the benchmarking set
        diff_meth_out <- extract_DMR_phase(dataset[[3]],min_number_CpGs,dataset[[2]],max_distance_between_DMRs)

        #assign the presence absence array to a variable (this has the form 00111000111 where 1 is differentially methylated)
        diff_meth_locs<-diff_meth_out[[1]]
        #assign the boundaries to a variabe (this has the form ((1,100),(300,600)))
        dfs_locs<-diff_meth_out[[2]]



        #record the time so that we can keep track of how long each method takes
        ptm <- proc.time()

        #run methylkit returning the predicted DMRs and per base likiehood of differenential methylation
        methylkit_results<-methylkit_for_bench(dataset,number_of_samples,number_of_replicas,output_path)
        #calculate the ROC curve for methylkit, this also include the AUC and DMR efficiency details
        methylkit_ROC <-par_methylkit_benchmark(FDR_cutoffs,methylkit_results$pvalue,dataset[[2]],diff_meth_locs,max_distance_between_DMRs,min_number_CpGs,methylkit_results$meth.diff,methykit_cutoff)
        

        #record the time that methylkit took
        methykit_time<-(proc.time() - ptm)[1]
        names[[ind]]<-'MethylKit'
        ROC[[ind]]<-methylkit_ROC
        ts[[ind]]<-methykit_time
        ind<-ind+1
        ptm <- proc.time()
        if(number_of_replicas>1){
        #run methylsig of the benchmark data
        methsig_results <- methsig_for_bench(dataset,number_of_samples,number_of_replicas,dataset[[2]])
        #create the ROC, AUC and DMR efficiencies for the methylsig

        methsig_ROC <-par_methsig_benchmark(Methsig_cutoffs,methsig_results@results[,1],dataset[[2]],diff_meth_locs,max_distance_between_DMRs,min_number_CpGs,methsig_results@results[,"meth.diff"],FDR_cutoffs)
        methsig_DMR_score <-par_methsig_benchmark_DMRs(Methsig_cutoffs,methsig_results@results[,1],dataset[[2]],dfs_locs,max_distance_between_DMRs,min_number_CpGs,methsig_results@results[,"meth.diff"],FDR_cutoffs)

        #record the runtime for methysig
        methsig_time<-(proc.time() - ptm)[1]
        names[[ind]]<-'MethylSig'
        ROC[[ind]]<-methsig_ROC
        ts[[ind]]<-methsig_time
        ind<-ind+1
        }
        #reset the timer
        ptm <- proc.time()

        if(number_of_replicas>1){
        #run bsmooth on the benchmark data
        bsmooth_obj  <- bsmooth_for_bench(dataset[[1]],number_of_samples,number_of_replicas)
        #save the bsmooth and fishers exact test results and store in seperate variables.
        bsmooth_results <- bsmooth_obj[[1]]
        fish_result <- bsmooth_obj[[2]]
        fish_pvals<-bsmooth_obj[[2]]$results
        #create the fishers exact ROC curve
        fish_ROC <- par_fisher_benchmark(fish_pvals[,1],FDR_cutoffs,diff_meth_locs,dataset[[2]],min_number_CpGs)
        
        #create the bsmooth ROC
        bsmooth_ROC <-par_bsmooth_benchmark(bsmooth_results,Methsig_cutoffs,diff_meth_locs,dataset[[2]],bsmooth_second_cutoff,min_number_CpGs)
        
        #create the bsmooth DMR efficiency
        bsmooth_DMR_score <-par_bsmooth_benchmark_DMRs(bsmooth_results,Methsig_cutoffs,dfs_locs,dataset[[2]],bsmooth_second_cutoff,min_number_CpGs)
        #save the bsmooth time
        bsmooth_time<-((proc.time() - ptm)[1])-bsmooth_obj[[5]]
        fish_time<-((proc.time() - ptm)[1])-bsmooth_obj[[6]]
        names[[ind]]<-'Fishers exact'
        ROC[[ind]]<-fish_ROC
        ts[[ind]]<-fish_time
        ind<-ind+1
        names[[ind]]<-'BSmooth'
        ROC[[ind]]<-bsmooth_ROC
        ts[[ind]]<-bsmooth_time
        #reset the timer
        ptm <- proc.time()
        }

        names_for_plot<-array(unlist(names), dim = c(nrow(names[[1]]), ncol(names[[1]]), length(names)))
        times_for_plot<-array(unlist(ts), dim = c(nrow(ts[[1]]), ncol(ts[[1]]), length(ts)))
        
        #plot the AUC/ROC/efficieny results
        plot_AUC(ROC,names_for_plot,output_path)
        plot_ROC(ROC,names_for_plot,output_path)
        plot_Time(times_for_plot,names_for_plot,output_path)
        save(ROC,ts,dfs_locs,file=paste0(output_path,".RData"))
        sink()
        return(list(ROC,ts,output_path))
}


###############################################################################
#takes the results from each of the benchmarkings and plots the results on an #
#ROC curve so that the performance can be compared.                           #
###############################################################################
plot_ROC <- function (ROCs,names,filename){
        #get the number of techniques that are going to be tested
        s <- length(ROCs)

        #open a pdf to save the results
        pdf(paste(filename,".ROC.pdf",sep=""))

        #loop through all of the techniques
        for(i in 1:s){
                #check if its the first and plot the ROC or else add to the ROC plot.
                if(i==1){
                plot(1-as.double(ROCs[[i]][6,]),as.double(ROCs[[i]][5,]),,ylim=c(0,1),xlim=c(0,1),col=i, main="ROC performance of methylation software", xlab="1-specificity", ylab="sensitivity",type='s')
                }else{
                lines(1-as.double(ROCs[[i]][6,]),as.double(ROCs[[i]][5,]),col=i,ylim=c(0,1),xlim=c(0,1),type='s')
                }
                #insert the legend to the bottom right.
                legend("bottomright", names, col = seq(1,s), lwd = 1, title="key")
        }
        #close the pdf
        dev.off()
}

###############################################################################
#pliots the efficiecy with which the DMRs are detected                        #
###############################################################################
plot_efficiency <- function (DMRs,names,filename){
        s <- length(DMRs)
        pdf(paste(filename,".efficiency.pdf",sep=""))
        for(i in 1:s){
                if(i==1){
                plot(as.double(DMRs[[i]][1,]),col=i,ylim=c(0,1), main="efficiency of methylation software", ylab="efficiency",type='s')
                }else{
                lines(as.double(DMRs[[i]][1,]),col=i,type='s')
                }
                legend("topright", names, col = seq(1,s), lwd = 1, title="key")
        }
        dev.off()
}

###############################################################################
#plots the area under the curve analysis for the benchmark                    #
###############################################################################
plot_AUC <- function (ROCs,names,filename){
        s <- length(ROCs)
        pdf(paste(filename,".AUC.pdf",sep=""))
        AUC<-matrix(data=0,nrow=1,ncol=s)
        for(i in 1:s){
                AUC[i]<-calculate_AUC(1-as.double(ROCs[[i]][6,]),as.double(ROCs[[i]][5,]))
        }
        barplot(AUC, main="Area under the curve analysis",xlab="DMR detection method",names.arg=names) 
        dev.off()
}

###############################################################################
#plots the time taken by each technique                                       #
###############################################################################
plot_Time <- function (times,names,filename){
        pdf(paste(filename,".times.pdf",sep=""))
        barplot(times, main="Runtime of each of the techniques",xlab="DMR detection method",ylab="Time for execution (in secs)",names.arg=names) 
        #legend("topright", names, col = seq(1,s), lwd = 1, title="key")
        dev.off()
}

###############################################################################
#calculates the AUC for a given ROC                                           #
###############################################################################
calculate_AUC <- function(x,y){
  l<-length(x)
  total = x[1]*(y[1]/2)
 for(i in 2:l){        
    x_dist = x[i]-x[i-1]
    y_dist = y[i]
    area = x_dist*y_dist
    total = total + area
  }
  total = total + ((1-x[l])*y[l])
  return(total)
}


###############################################################################
#extracts the DMRs from the state vector                                      #
###############################################################################
extract_DMR_phase <-function (states_a,thresh,locs,max_distance){
diff_meth <- matrix(data=0,nrow=1,ncol=length(states_a))
dmr_locs<-list()
dmr_locs_ind<-1
on<-0
        off<-0
        count<-0
        points=list()
        for(i in 2:length(states_a)){
                if((states_a[i] == 1)){
                        if(on == 0){
                                on<-i
                                count<-count+1
                        }else{
                                off<-i
                                count<-count+1  
                        }
                }else if(on != 0){
                        if((off != 0)){
                                diff_meth[on:off]<-1
                                dmr_locs[[dmr_locs_ind]]<-c(locs[on],locs[off])
                                dmr_locs_ind<-dmr_locs_ind+1
                                on<-0
                                off<-0
                        count<-0
                        }else{
                                on<-0
                                off<-0
                                count<-0
                        }
                }
     }
        if((off != 0)&(on != 0)){
                diff_meth[on:off]<-1
                dmr_locs[[dmr_locs_ind]]<-c(locs[on],locs[off])
                dmr_locs_ind<-dmr_locs_ind+1
        }
        return(list(diff_meth,dmr_locs))
}

get_DMRs <- function(states){
        
}

score_overlap <-function (real,pred,cutoff){
        l <- length(real)
        l <- length(pred)
        TP<-0
        FN<-0
        FP<-0
        TN<-0
        for(i in 1:l){
                if(real[i]==1 & pred[i]==1){
                        TP<-TP+1
                }else if(real[i]==1 & pred[i]==0){
                        FN<-FN+1
                }else if(real[i]==0 & pred[i]==1){
                        FP<-FP+1
                }else{
                        TN<-TN+1
                }
        }
        sensitivity = TP/(TP+FN)
        specificity = TN/(TN+FP)
        result<-list(TP,FP,FN,TN,sensitivity,specificity,cutoff)
        return(result)
}

isoverlap <- function (loc1,loc2){
        a1<-loc1[1]
        a2<-loc1[2]
        b1<-loc2[1]
        b2<-loc2[2]
        if(a1>b2){
                return(0)
        }else if(a2<b1){
                return(0)
        }else{
                return(1)
        }
}

score_overlap_DMR <-function (real,pred,cutoff){
        l_real <- length(real)
        l_pred <- length(pred)
        R_t<-0
        R_f<-0
        #loop through each real DMR
        #loop through each predicted DMR
        #check for the overlap
        #if its above threshold TP+1
        if(l_pred==0){
                R_f<-l_real
        }else{
        for(i in 1:l_real){
                id<-0

                for(j in 1:l_pred){
                        is_overlap<-isoverlap(real[[i]],pred[[j]])
                        if(is_overlap==1){
                                id<-1
                        }
                }
                if(id==0){
                        R_f <- R_f + 1
                }else{
                        R_t <- R_t + 1
                }
        }
        }

        result <- R_t/(R_f+l_pred)
        return(list(result,R_t,R_f,l_pred,l_real,pred,cutoff))
}

args <- commandArgs(trailingOnly = TRUE)
if(length(args)==18){
        if(as.numeric(args[3]) == 1){
                cat("#######################################################################################################\n")
                cat("# WARNING!                     Only Methylkit works with 1 replica.                           WARNING!#\n")
                cat("#######################################################################################################\n")
        }
        
  rd<-auto_init_benchmark(args[1],as.numeric(args[2]),as.numeric(args[3]),as.numeric(args[4]),as.numeric(args[5]),as.numeric(args[6]),as.numeric(args[7]),as.numeric(args[8]),as.numeric(args[9]),as.numeric(args[10]),as.numeric(args[11]),as.numeric(args[12]),as.numeric(args[13]),as.numeric(args[14]),as.numeric(args[15]),as.numeric(args[16]),args[17],as.numeric(args[18]),'binomial')
          cat("\n##########################BENCHMARK COMPLETE############################\n")
      cat(paste0("\nThe ROC analysis is here: \n",rd[[3]],".ROC.pdf\n"))
      cat(paste0("\nThe AUC analysis is here: \n",rd[[3]],".AUC.pdf\n"))
      cat(paste0("\nThe Runtime analysis is here: \n",rd[[3]],".times.pdf\n"))
}else if(length(args)==1){
    if(args[1]=='interactive'){
       rd<-init_interactive_benchmark()
       cat("\n##########################BENCHMARK COMPLETE############################\n")
      cat(paste0("\nThe ROC analysis is here: \n",rd[[3]],".ROC.pdf\n"))
      cat(paste0("\nThe AUC analysis is here: \n",rd[[3]],".AUC.pdf\n"))
      cat(paste0("\nThe Runtime analysis is here: \n",rd[[3]],".times.pdf\n"))
    }
}




###################################################################################################################################################################################################################


