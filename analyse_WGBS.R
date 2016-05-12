###############################################################################
#ANALYSE WGBS SEQUENCE DATA TO DEFINE PARAMS FOR HMM SIM-Owen Rackham, MRC CSC#
###############################################################################
#This scripts loads a WGBS dataset and performs some simple analysis in order #
#to paramerterise a simulation to produce data with similar characteristics.  #
#for more information contact owen.rackham@imperial.ac.uk                     #
###############################################################################



###############################################################################
#interactively initiate the function                                          #
###############################################################################
init_interactive_analysis <- function (){
  cat("#######################################################################################################\n")
  cat("#                              WELCOME TO THE WGBS INTERACTIVE MODE                                   #\n")
  cat("#######################################################################################################\n\n")
  cat("#      For further information please visit www.wgbssuite.org.uk for the docs or contact details      #\n")
  cat("#######################################################################################################\n")
#get the variables from the command line/tmp/humanWGBS/1.in
filename <-get_text("Please enter the filename you wish to analyse (default: ./test.in)","test.in")
analysis_folder <-get_text("Please enter the analysis folder you wish to store the results (default: /tmp/)","/tmp/")
analysis_name <-get_text("Please enter a name you wish to store the results under (default: example_WGBS_run)","example_WGBS_run")
output_path <- paste0(analysis_folder,analysis_name)
number_of_samples <-get_variable("Please enter the number_of_samples in your dataset (default: 2)",2)
number_of_replicas <-get_variable("Please enter the number_of_replicas in your dataset (default: 2)",2)
chr_column <-get_variable("Please enter the column number that contain chromosome information (default: 1)",1)
position_column <-get_variable("Please enter the column number that contains the position information (default: 2)",2)
strand_column <-get_variable("Please enter the column number that contians the strand information (default: 4)",4)
starting_column <-get_variable("Please enter the column number on which the first replica starts (default: 5)",5)
ct_column <-get_variable("Please enter the column number relative to the sample start that contains the c2t count info (default: 0)",0)
total_column <-get_variable("Please enter the column number relative to the sample start that contains the total count info (default: 1)",1)
total_columns_per_sample <-get_variable("Please enter the total number of columns per replica (default: 2)",2)

#run the simulation script

d<-extract_data(filename,number_of_samples,number_of_replicas,chr_column,position_column,strand_column,starting_column,ct_column,total_column,total_columns_per_sample)
summarise_real_data(d,number_of_replicas,number_of_samples,filename,output_path)
#run the summary script
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
#extract the data from the file                                               #
###############################################################################
extract_data <- function (filename,number_of_samples,number_of_replicas,chr_column,position_column,strand_column,starting_column,ct_column,total_column,total_columns_per_sample){
  cat("#######################################################################################################\n")
  write("Extracting data from the files.", stderr())
	pb <- txtProgressBar(style=3)
	data <- as.matrix(read.table(filename,header = T))
	number_of_sites <- dim(data)[1]
	total_samples <- number_of_samples*number_of_replicas
	d_size <- 2+(4*total_samples)
	d <- matrix(data=0,nrow=number_of_sites,ncol=d_size)
	d[,1] <-as.numeric(data[,position_column])

	current_index<-3

	for(i in 1:number_of_samples){
		for(j in 1:number_of_replicas){
			index <- ((i-1)*number_of_replicas)+j
			setTxtProgressBar(pb, (index/total_samples))
			x_ind <- ((index-1)*total_columns_per_sample)+starting_column
			d[,current_index] <- as.numeric(data[,(x_ind+ct_column)])
			current_index <- current_index+1
			d[,current_index] <- as.numeric(data[,(x_ind+total_column)])-as.numeric(data[,(x_ind+ct_column)])
			current_index <- current_index+1
			d[,current_index] <- as.numeric(data[,(x_ind+total_column)])
			current_index <- current_index+1
			d[,current_index] <- as.numeric(data[,(x_ind+ct_column)])/as.numeric(data[,(x_ind+total_column)])
			current_index <- current_index+1
		}
	}
	
  return(d)
}

###############################################################################
#Function to summarise the  data from the simulation. This should show: 1)    #
#distribution of methlated/unmethylated sites 2) distribution of distances    #
#between CpGs 2) Distribution of the size of differentially methylated regions#
# 3) Distribution of read counts 4) 
###############################################################################
summarise_real_data <- function(generated_set,number_of_replicas,number_of_samples,filename,output_path){
  #define the layout of the figures (actually has no effect when saving to pdf)
  mat <- matrix(1:6, 2, 3)
  layout(mat)
  #take the name from the simulated set and append the pdf extension
  
  pdf(paste(output_path,".analysis_of_real_data.pdf",sep=""))
  #create the graphs

  read_count_distribution(generated_set,number_of_replicas,number_of_samples)
  methylation_settings <- methylation_distribution(generated_set,number_of_replicas,number_of_samples)
  meths <- methylation_settings[[5]]
   cat("\n#######################################################################################################\n")
  cat("##########################THE ANALYSIS IS COMPLETED, SEE BELOW.########################################\n")
  cat("#######################################################################################################\n\n")
  cat("The following information can be used to parametrise the simulation of data to match your dataset","\n")
  cat(paste("The average difference in methylation proportion: ",methylation_settings[[1]]),"\n")
  cat(paste("The standard deviation in difference in methylation proportion: ",methylation_settings[[2]]),"\n")
  cat(paste("The max difference in methylation proportion: ",methylation_settings[[4]]),"\n")
  cat(paste("The probability of success in a methylated region: ",methylation_settings[[3]][[1]]),"\n")
  cat(paste("The probability of success in a unmethylated region: ",methylation_settings[[3]][[2]]),"\n")
  cat(paste("The average coverage in a methylated region: ",round(methylation_settings[[3]][[3]])),"\n")
  cat(paste("The average coverage in an unmethylated region: ",round(methylation_settings[[3]][[4]])),"\n")
  distances<-distance_distribution(generated_set)
  distance_coefs<-distances[[1]]
  l<-length(distances[[2]])

  cat(paste("The fast decay rate for CpGs: ",distance_coefs[[1]]),"\n")
  cat(paste("The slow decay rate for CpGs: ",distance_coefs[[2]]),"\n\n")
  dev.off()
  cat(paste0("Further details can be found in: ",output_path,".analysis_of_real_data.pdf"),"\n\n")
  cat(paste0("nb. in Unix: evince ",output_path,".analysis_of_real_data.pdf"),"\n")
  cat(paste0("or in OsX: open ",output_path,".analysis_of_real_data.pdf"),"\n\n")
  cat(paste("To create a binomially simulated dataset with the default settings execute the following:\n"))
  command<- paste0("Rscript simulate_WGBS.R ",5000," ",methylation_settings[[3]][[1]]," ",methylation_settings[[3]][[2]]," ",0.2," ",0.2," ",methylation_settings[[3]][[4]]," ",methylation_settings[[3]][[4]]," ",3," ",2," ",0.05," ",0.5," ",distance_coefs[[1]],",",distance_coefs[[2]]," ",output_path," ","binomial","\n")
	cat(command)
  cat(paste("To create a negative binomially simulated dataset with the default settings execute the following:\n"))
  command<- paste0("Rscript simulate_WGBS.R ",5000," ",methylation_settings[[3]][[1]]," ",methylation_settings[[3]][[2]]," ",0.2," ",0.2," ",methylation_settings[[3]][[4]]," ",methylation_settings[[3]][[4]]," ",3," ",2," ",0.05," ",0.5," ",distance_coefs[[1]],",",distance_coefs[[2]]," ",output_path," ","truncated","\n")
  cat(command)
}

###############################################################################
#fit spacing distribution to real data                                        #
###############################################################################
fit_distance_distribution <- function(length_set,nearby){
  #seperate the nearby CpG and the total set and then fit two seperate distributions to these
  #two sets.
  nearby_cpgs <- length_set[length_set<nearby]
  hist_nearby <- hist(nearby_cpgs,col=rgb(1,0,0,0.2),breaks=(nearby/2),ylab="frequency",xlab="coverage",main="histogram of distances between CpGs locally")
  x<-nearby_cpgs
  sink(file="/dev/null")
  library(MASS)
  sink()
  fitted_exponential<-fitdistr(x, 'exponential')
  fast_decay <- coef(fitted_exponential)
  nearby_cpgs <- length_set[length_set>nearby]
  hist_nearby <- hist(nearby_cpgs,col=rgb(1,0,0,0.2),ylab="frequency",xlab="coverage",main="histogram of distances between CpGs locally")
  x<-nearby_cpgs
  fitted_exponential<-fitdistr(x, 'exponential')
  slow_decay<-coef(fitted_exponential)
  return(list(fast_decay,slow_decay))
}




###############################################################################
#Distribution of proportion of methylation                                    #
###############################################################################
methylation_distribution <- function(generated_set,number_of_replicas,number_of_samples){
  #calcualte the proportion of methylated vs demethylated reads in each replicate
  #and plot the valus as a histogram

  all_meth<-vector()
  all_fail<-vector()
  all_success<-vector()
  all_coverage<-vector()
  group_meth<-list()
  comp_meth<-list()
  write("\nCalculating methylation distributions.\n",stderr())
  pb <- txtProgressBar(style=3)
  for(i in 1:number_of_samples){
    for(j in 1:number_of_replicas){
      index <- ((i-1)*number_of_replicas)+j
      setTxtProgressBar(pb, (index/(number_of_replicas*number_of_samples)))
      col <- (((i-1)*(number_of_replicas*4))+((j-1)*4)+1)+2
      meth_prop <- (generated_set[,col]/generated_set[,col+2])
      diff_as_prop <- (generated_set[,col]-generated_set[,col]+1)/generated_set[,col+2]
      meth_prop[is.na(meth_prop)] <- 0.0000000001
      meth_prop[meth_prop>=1] <- 0.99999999
      hist(meth_prop,col=rgb(0,0,0,0.2),ylab="frequency",xlab="proportion of methylated CpGs",main=paste("histogram of methylation proportion in sample",i,"replicate",j))
      all_meth<-c(all_meth,as.numeric(meth_prop))
      all_success<-c(all_fail,as.numeric(generated_set[,col]))
      all_fail<-c(all_fail,as.numeric(generated_set[,col+1]))
      all_coverage<-c(all_coverage,as.numeric(generated_set[,col+2]))
      if(j==1){
      		group_meth[[i]]<-meth_prop
      		comp_meth[[i]]<-meth_prop
      	}else{
      		group_meth[[i]]<-c(meth_prop,group_meth[[i]])
      		comp_meth[[i]]<-meth_prop+comp_meth[[i]]
      	}
    }
  }


  comp_meth[[1]] <- comp_meth[[1]]/number_of_replicas
  comp_meth[[2]] <- comp_meth[[2]]/number_of_replicas


  av_diff <- abs(comp_meth[[1]]-comp_meth[[2]])
  av_diff_mean <- mean(av_diff,na.rm=T)
  av_diff_sd <- sd(av_diff,na.rm = T)
  max_diff <- max(av_diff,na.rm = T)
  hist(av_diff,main="distribution of meth proportion diffs")

  hist((all_fail/all_coverage),main="failures over coverage")
  hist((all_success/all_coverage),main="success over coverage")
  meth_fitted_results <- fit_methylation_distribution(all_meth,0.5,all_coverage)
  #negbinomial_for_success <- fit_negativebinomial_to_methylation_distribution(all_meth,0.2,0.8,all_success, "methylated reads")
  return(list(av_diff_mean,av_diff_sd,meth_fitted_results,max_diff,av_diff))
}

###############################################################################
#fit spacing distribution to real data                                        #
###############################################################################
fit_methylation_distribution <- function(meth_set,split,coverage){
  #seperate the nearby CpG and the total set and then fit two seperate distributions to these
  #two sets.
  methylated <- meth_set[meth_set>split]
  un_methylated <- meth_set[meth_set<split]
  methylated_average <- mean(methylated,na.rm = T)
  un_methylated_average <- mean(un_methylated,na.rm = T)
  coverage[coverage>100]<-100
  methylated_coverage <- coverage[meth_set>split]
  unmethylated_coverage <- coverage[meth_set<split]
  methylated_coverage_average <- mean(methylated_coverage,na.rm = T)
  unmethylated_coverage_average <- mean(unmethylated_coverage,na.rm = T)
  hist(methylated_coverage,main="methylated coverage",breaks=50)
  hist(unmethylated_coverage,main="unmethylated_coverage",breaks=50)

  return(list(methylated_average,un_methylated_average,methylated_coverage_average,unmethylated_coverage_average))
}

###############################################################################
#fit negative binomial to real data                                        #
###############################################################################
fit_negativebinomial_to_methylation_distribution <- function(meth_set,split_low,split_high,fails,name){
  #seperate the nearby CpG and the total set and then fit two seperate distributions to these
  #two sets.
  sink(file="/dev/null")
  library(MASS)
  sink()

  methylated_fails <- fails[meth_set>split_high]
  middle_fails <- fails[meth_set>split_low & meth_set<split_high]
  un_methylated_fails <- fails[meth_set<split_low]

  hist(fails,main=paste("all distribution",name),breaks=200)
  hist(methylated_fails,main=paste("methylated distribution",name),breaks=200)
  hist(middle_fails,main=paste("middle methylated distribution",name), breaks=200)
  hist(un_methylated_fails,main=paste("un methylated distribution",name), breaks=200)

  meth_fitted_nbin<-fitdistr(methylated_fails,  "negative binomial")
  middle_fitted_nbin<-fitdistr(middle_fails,  "negative binomial",)
  unmeth_fitted_nbin<-fitdistr(un_methylated_fails,  "negative binomial")

  meth_fitted_nbin_gz<-fitdistr(methylated_fails[methylated_fails>0],  "negative binomial")
  middle_fitted_nbin_gz<-fitdistr(middle_fails[middle_fails>0],  "negative binomial",)
  unmeth_fitted_nbin_gz<-fitdistr(un_methylated_fails[un_methylated_fails>0],  "negative binomial")

  methylated_coefs<-coef(meth_fitted_nbin)
  middle_coefs<-coef(middle_fitted_nbin)
  unmethylated_coefs<-coef(unmeth_fitted_nbin)

  methylated_coefs_gz<-coef(meth_fitted_nbin_gz)
  middle_coefs_gz<-coef(middle_fitted_nbin_gz)
  unmethylated_coefs_gz<-coef(unmeth_fitted_nbin_gz)

  hist(rnbinom(5000, size = methylated_coefs[1], mu = methylated_coefs[2]),main=paste("fitted methylated distribution",name,methylated_coefs))
  hist(rnbinom(5000, size = middle_coefs[1], mu = middle_coefs[2]),main=paste("fitted middle distribution",name,middle_coefs))
  hist(rnbinom(5000, size = unmethylated_coefs[1], mu = unmethylated_coefs[2]),main=paste("fitted demethylated distribution",name,unmethylated_coefs))

  hist(rnbinom(5000, size = methylated_coefs_gz[1], mu = methylated_coefs_gz[2]),main=paste("fitted methylated distribution greater than 0",name,methylated_coefs_gz))
  hist(rnbinom(5000, size = middle_coefs_gz[1], mu = middle_coefs_gz[2]),main=paste("fitted middle distribution greater than 0",name,middle_coefs_gz))
  hist(rnbinom(5000, size = unmethylated_coefs_gz[1], mu = unmethylated_coefs_gz[2]),main=paste("fitted unmethylated distribution greater than 0",name,unmethylated_coefs_gz))

  return(list(methylated_fails,un_methylated_fails,methylated_coefs,middle_coefs,unmethylated_coefs,methylated_coefs_gz,middle_coefs_gz,unmethylated_coefs_gz))
}

###############################################################################
#fit beta to real data                                        #
###############################################################################
fit_beta_to_distribution <- function(meth_set,name){
  sink(file="/dev/null")
  library(MASS)
  sink()
  hist(meth_set,main=paste("methlaytion distribution",name),breaks=200)

}


###############################################################################
#CpG locations on a line                                                      #
###############################################################################
methylation_postions <- function(generated_set){
  #plot the location of each CpG as a transparent bar on a line
  plot(generated_set[,1],matrix(data=1,ncol=1,nrow=length(generated_set[,1])),pch='|',col=rgb(0,0,0,0.01),ylab="",xlab="location",yaxt='n', ann=FALSE)

}
###############################################################################
#Distribution of distance between sites                                       #
###############################################################################
distance_distribution <- function(generated_set){
  #get the number of CpGs
  test<-sort(unique(generated_set[,1]))
  CpG_length <- length(test)
  #initialise a matrix to store the lengths in

  length_set <- matrix(data=0,nrow=CpG_length,ncol=1)
  
  #loop through and save the lengths between each site
  for (i in 2:CpG_length){
    length_set[i]<-test[i]-test[i-1]
  }
  #plot the distribution
  hist(length_set,col=rgb(1,0,0,0.2),breaks=100,ylab="frequency",xlab="coverage",main="histogram of distances between CpGs")
  distance_coefs<-fit_distance_distribution(length_set,200)
  mid_distance_coefs<-fit_distance_distribution(length_set,1000)

  return(list(distance_coefs,length_set,mid_distance_coefs))
}

###############################################################################
#Distribution of read counts                                                  #
###############################################################################
read_count_distribution <- function(generated_set,number_of_replicas,number_of_samples){
  #plot a histogram for each read count set.
  write("\nCalculating read count distribution.\n",stderr())
  pb <- txtProgressBar(style=3)
  for(i in 1:number_of_samples){
    for(j in 1:number_of_replicas){
      index <- ((i-1)*number_of_replicas)+j
      setTxtProgressBar(pb, (index/(number_of_replicas*number_of_samples)))
      col <- ((i-1)*(number_of_replicas*4))+((j-1)*4)+3
      hist(as.numeric(generated_set[,col+2]),col=rgb(1,0,0,0.2),ylab="frequency",xlab="coverage",main=paste("histogram of read counts from sample",i,"replicate",j))
    }
  }
}


args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=8){
  init_interactive_analysis()
}else{
  init(args)
}