###############################################################################
#SIMULATED WGBS SEQUENCE DATA USING A HIERACHICAL HMM - Owen Rackham, MRC CSC #
###############################################################################
#This script is designed to be run using RScript from the terminal, it will   #
#produce a simulated dataset for WGBS containing 2 sample types WITH 3 reps.  #
#At present both this and the number of states in the HMM are fixed but v2    #
#will contain a multistate HMM for methylation status allowing smooth         #
#transition from methylated to de-methylated states as is seen in real data.  #
#particularly in the promoter region of genes. V2 will also allow you to set  #
#the number of sample types and the number of replicates in each case.        #
#for more information contact owen.rackham@imperial.ac.uk
###############################################################################

###############################################################################
#main body function which takes lots of params that control the simulated data#
###############################################################################
generate_sim_set <- function(n,m,Pi_m,Pi_d,mean_m,mean_d,prob_m,prob_d,error_m,error_d,phase_diff,outfile,type_of_locations,number_of_replicas,number_of_samples,rates_for_HMM_for_CpG_locations,transition_size,balance,output_path,type,tranMPar, tranNPar, nonconvPar) {
#init_simulation <- function (theta,n_mean,n_standard_deviation,number_of_CpGs,probability_of_success_in_differentially_methylated_region,probability_of_success_in_non_differentially_methylated_region,error_rate_in_differentially_methylated_region,mean_number_of_reads_in_differentially_methylated_region,mean_number_of_reads_in_non_differentially_methylated_region,number_of_replicas,number_of_samples,phase_diff,balance,cpg_matrix,){

  pdf(paste(output_path,"data.pdf",sep=""))
  #load required libraries
  sink(file="/dev/null")
  library(zoo)
  library(HiddenMarkov)
  sink()
  #create a set of CpG locations
  locs<-create_simulated_locations(n,Pi_m,rates_for_HMM_for_CpG_locations)

#    tranMPar; The transition probability in a methylated region:  
#    tranNPar; The transition probability in an unmethylated region:  0.0577772070946612
#    message(tranMPar)
#    message(tranNPar)
                                        #simulate the state transition based on the location of the CpGs
    a<-simulate_state_transition((n*m),c(tranMPar[1],0.99,tranNPar[1],0.99),locs,0.5,transition_size, rates = c(tranMPar[2], 1.895e-02, tranNPar[2], 1.895e-02))
    
    write.table(a, file = "a.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE) 
   

  #extract the postions of the blocks
  state_blocks<-find_the_blocks(a)


  #this needs to be changed to reflect the changes in the find_the_blocks function - TJH
  
  #set the percentage of DM in each block type
  percs_for_diff<-c(0.5,0,0,0)

  #update the blocks to be differentially methylated
  diff_methed <- make_differential(state_blocks,percs_for_diff,a)

                                        # this bit is optional I guess, but it seems to me to make more sense that if a region is DM, then the transition states on either side are also DM. Otherwise, depending on your parameters, you could have the transition spikes be more methylated than the differentially methylated region for one of the experimental conditions. If you allow differential behaviour in non-methylated regions, or at random transition locations (i.e., if any value except the first in percs_for_diff to be non-zero) then this code needs some more tweaking. - TJH
  spla <- split(as.integer(a[1,]), rep(1:nrow(state_blocks), state_blocks[,4]))
  spldm <- split(as.integer(diff_methed[1,]), rep(1:nrow(state_blocks), state_blocks[,4]))
  dmBlocks <- which(sapply(spldm, function(x) all(x == 1)))
   lT <- (sapply(spla[dmBlocks - 1], function(x) all(!is.null(x)) & all(x == 4))); leftTransition <- dmBlocks[which(lT)] - 1
    rT <- (sapply(spla[dmBlocks + 1], function(x) all(!is.null(x)) & all(x == 2))); rightTransition <- dmBlocks[which(rT)] + 1
    spldm[leftTransition] <- lapply(spldm[leftTransition], function(x) rep(1, length(x)))
    spldm[rightTransition] <- lapply(spldm[rightTransition], function(x) rep(1, length(x)))

# according to the balance, makes change in differential methylation all in one direction (balance = 0), all in the other direction(balance = 1) or something in between.
    mult <- sample(c(-1, 1), size = length(dmBlocks), replace = TRUE, prob = c(balance, 1-balance))
    mult_n <- mult_m <- spldm;

    zMult <- c(dmBlocks[intersect(which(lT), which(mult == -1))] - 1,
               dmBlocks[mult == -1],
               dmBlocks[intersect(which(lT), which(mult == -1))] + 1)
    mult_m[zMult] <- lapply(mult_m[zMult], function(x) rep(0, length(x)))
    mult_m <- unlist(mult_m)
    
    zMult <- c(dmBlocks[intersect(which(lT), which(mult == 1))] - 1,
                   dmBlocks[mult == 1],
                   dmBlocks[intersect(which(lT), which(mult == 1))] + 1)
    mult_n[zMult] <- lapply(mult_n[zMult], function(x) rep(0, length(x)))
    mult_n <- unlist(mult_n)

    
    diff_methed[1,] <- unlist(spldm)
    
    write.table(diff_methed, file = "diff_methed.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)



  #create the simulated reads methylated/unmethylated at each CpG, at the moment this is hard coded to be 3 replicates of each type
  #The phase diff param control how different the methylation is in the differentially methylated regions.
  ds <- (number_of_replicas*number_of_samples)*4
  d <- matrix(data=0,nrow=ds,ncol=n)
  idx <- -(number_of_replicas*4) + 1
  
      assignProbs <- function(prob, diff, probMod) {
        pm <- rep(prob, length(probMod))
        pm[probMod == 1] <- prob + diff
        pm <- pmin(1, pm); pm <- pmax(0, pm)
        pm
    }

    prob_d_mod <- assignProbs(prob_d, phase_diff[1], mult_m)
    prob_m_mod <- assignProbs(prob_m, -phase_diff[2], mult_m)

    prob_d_nonmod <- assignProbs(prob_d, phase_diff[1], mult_n)
    prob_m_nonmod <- assignProbs(prob_m, -phase_diff[2], mult_n)
  
                                        #prob_d_mod<-hypo_hyper(prob_d,phase_diff[2],a,balance)
#    prob_d_mod<-hypo_hyper_diffs(prob_d,phase_diff[2],diff_methed,balance)
#    prob_m_mod<-hypo_hyper_diffs(prob_m,phase_diff[2],diff_methed,balance)
#    hist(prob_m_mod,main="prob_d_mod")
                                        #prob_d_nonmod<-hypo_hyper(prob_d,phase_diff[1],a,balance)
#    prob_d_nonmod<-hypo_hyper_diffs(prob_d,phase_diff[1],diff_methed,balance)
#    prob_m_nonmod<-hypo_hyper_diffs(prob_m,phase_diff[1],diff_methed,balance)
    
#    hist(prob_d_nonmod,main="prob_d_nonmod")

    save(prob_d, prob_m, prob_d_mod, prob_m_mod, prob_d_nonmod, prob_m_nonmod, phase_diff, diff_methed, balance, mult_n, mult_m, a, file = "probmod.RData")

  probs<-list()
    # probs[[1]][[2]] should be (prob_m_mod+prob_d_mod) / 2, not ((prob_m_nonmod+prob_d_nonmod)/2) - TJH
  probs[[1]]<-list(prob_m_mod,((prob_m_mod+prob_d_mod)/2),prob_d_mod,((prob_m_mod+prob_d_mod)/2))
  probs[[2]]<-list(prob_m_nonmod,((prob_m_nonmod+prob_d_nonmod)/2),prob_d_nonmod,((prob_m_nonmod+prob_d_nonmod)/2))
  errors<-list(error_m,((error_d+error_m)/2),error_d,((error_d+error_m)/2))
  means<-list(mean_m,((mean_d+mean_m)/2),mean_d,((mean_d+mean_m)/2))

  write("\nCreating the simulated dataset.\n",stderr())
  pb <- txtProgressBar(style=3)
  pdf("locd.pdf")
  for(i in 1:number_of_samples){
    idx <- idx + (4*number_of_replicas)
    for(j in 1:number_of_replicas){
      index <- ((i-1)*number_of_replicas)+j
      setTxtProgressBar(pb, (index/(number_of_replicas*number_of_samples)))
      idx_rep <- idx + (j*4) - 4
      #d[idx_rep:(idx_rep+3),] <- generate_replicat_methyl_nbin_data(a,probs,means,0,errors,locs,i)
      if(type=='binomial'){
          d[idx_rep:(idx_rep+3),] <- generate_replicat_methyl_bin_data(a,probs,means,0,errors,locs,i,output_path)
        }else if(type=='truncated'){
          d[idx_rep:(idx_rep+3),] <- generate_replicat_methyl_truncated_nbin_data(a,probs,means,0,errors,locs,i,20,output_path)
        }else{
          d[idx_rep:(idx_rep+3),] <- generate_replicat_methyl_nbin_data_model_3(a,probs,means,0,errors,locs,i,30)
        }
      if((i == 1)&&(j==1)){
          plot(locs,d[(idx_rep+3),],col=i,type='b',ylim=c(0,1),pch=diff_methed+1)
        }else{
          lines(locs,d[(idx_rep+3),],col=i,type='b',pch=diff_methed+1)
        }
    }
  }
  dev.off()
  
  # non-conversion
      if(!any(is.na(nonconvPar))) {
        Cs <- d[1:(nrow(d) / 4) * 4 - 3,]
        Ns <- d[1:(nrow(d) / 4) * 4 - 1,]
        nonconv <- rbeta(nrow(Cs), nonconvPar[1], nonconvPar[2])
        write.table(nonconv, file = "nonconv.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
        NC <- matrix(rbinom(prod(dim(Cs)), size = Ns - Cs, prob = matrix(nonconv, nrow = nrow(Cs), ncol = ncol(Cs))), nrow = nrow(Cs), ncol = ncol(Cs))
        d[1:(nrow(d) / 4) * 4 - 3,] <- Cs + NC
        d[1:(nrow(d) / 4) * 4,] <- d[1:(nrow(d) / 4) * 4 - 3,] / d[1:(nrow(d) / 4) * 4 - 1,]
    }
  
  d<-rbind(locs,diff_methed,d)
  #wrap up the simulated data into a single variable
  #transpose
  d<-t(d)
  outname <- paste(c(outfile, n,m,Pi_m,Pi_d,mean_m,mean_d,prob_m,prob_d,error_m,error_d,phase_diff[1],phase_diff[2],balance,".txt"), collapse = "_")
  #write to file
  write.table(d,outname, row.names = F,quote=F,sep = "\t")
  dev.off()
  return(list(locs,d,outname))
  
}


###############################################################################
#find the  start and end of the blocks in a vector of states                  #
###############################################################################
find_the_blocks<-function(a){
  l<-length(a)
  coords<-list()
  index<-1
  on<-1
  off<-0
 # this needs to go to l + 1 or the last cytosine is not included - TJH
  for(i in 2:(l + 1)){
      # check introduced to capture last cytosine - TJH
      if(a[i-1] != a[i] | i > 1){
          # I think the check for on = 0 is redundant because it never can be - TJH

          # some tweaks needed here - it seems to make more sense if the coordinate structure matches the 'a' vector - TJH
          # also adding ' + 1' to the width matches GRanges conventions - TJH
          off<-i - 1          
          coords[[index]]<-c(on,off,a[i-1],off-on + 1)
          index<-index+1
          on<-i
          off<-0
      }
  }
  coords_mat <- do.call(rbind,coords)
  return(coords_mat)
}

###############################################################################
#make differentially expressed by picking a preset amount of each state type  #
###############################################################################
make_differential<-function(state_blocks,percentage,a){
  l<-length(a)
  states<-unique(state_blocks[,3])
  n <-length(states)
  diff_meth <- matrix(data=0,nrow=1,ncol=l)
  for(i in 1:n){
    blocks <- state_blocks[state_blocks[,3]==i,]
    block_length <- sum(blocks[,4])
    cutoff_length <- percentage[i]*block_length
    so_far <- 0
    while((so_far <= cutoff_length)&&!is.null(dim(blocks)[1])){
      picked<- sample(1:dim(blocks)[1], 1)
      if(so_far+blocks[picked,4] < cutoff_length){
        so_far <- so_far + blocks[picked,4]
        diff_meth[blocks[picked,1]:blocks[picked,2]]<-1
      }
      blocks<- blocks[-(picked),]
    }
  }
  return(diff_meth)
}






###############################################################################
#Alter demethylated probs to be beta distributed                              #
###############################################################################
create_beta_probs<-function(a,prob_m,prob_d){
l<-length(a)
probs <- matrix(data=0,nrow=1,ncol=l)
probs[1]<-trun_beta_probs(prob_m,prob_d)

for(i in 2:l){
  if(a[i-1]==a[i]){
      probs[i] <- probs[i-1]
    }else{
      probs[i] <- trun_beta_probs(prob_m,prob_d)
    }
}
return(probs)
}

###############################################################################
#rbeta but not 0 or 1                              #
###############################################################################
trun_beta_probs<-function(prob_m,prod_d){


  a<-rbeta(1,0.4,0.4)
  if(a<prod_d){
    a<-prod_d
  }
  if(a>prob_m){
    a<-prob_m
  }

return(a)
}



###############################################################################
#Create the CpG locations                                                     #
###############################################################################
create_simulated_locations<-function(n,Pi,rates_for_HMM_for_CpG_locations){
#set up the variables
delta <- c(0, 1)

#create a HMM model that will simulate the gaps between CpGs. it is a 2 state HMM
#in order to simulate CpG islands and CpG deserts. TODO: multistates for CpG shores
#,cliffs etc etc etc to be added later.
x <- dthmm(NULL, Pi, delta, "exp", list(rate=rates_for_HMM_for_CpG_locations))

#simulate this HMM for n steps
x <- simulate(x, nsim=n)

#extract the length
l<-length(x$x)

#initialise a matrix to store the results
locs <- matrix(data=1,nrow=1,ncol=l)

#loop over the simulated output and create the locations
for(i in 2:l){
  locs[i]<-round(locs[i-1]+x$x[i])+1
}

return(locs)
}

###############################################################################
#Randomly create hyper/hypo methylation                                       #
###############################################################################
hypo_hyper<-function(prob_m,phase_diff,states,split){
l<-length(states)
current_state<-states[1]
  pos<-0
  neg<-0
prob_m_mod <- matrix(data=prob_m[1],nrow=1,ncol=l)

  if(length(prob_m) > 1){
      multiply<-plus_minus(split,phase_diff,prob_m[1])
    }else{
      multiply <- plus_minus(split,phase_diff,prob_m)
    }

for(i in 2:l){
  if((states[i]==2)&(current_state==1)){
   if(length(prob_m) > 1){
      multiply<-plus_minus(split,phase_diff,prob_m[i])
    }else{
      multiply <- plus_minus(split,phase_diff,prob_m)
    }
  }
  prob_m_mod[i] <-  multiply
  current_state<-states[i]
}

return(prob_m_mod)

}


###############################################################################
#Randomly create hyper/hypo methylation                                       #
###############################################################################
hypo_hyper_diffs<-function(prob_m,phase_diff,states,split){
l<-length(states)
current_state<-states[1]
  pos<-0
  neg<-0
prob_m_mod <- matrix(data=prob_m[1],nrow=1,ncol=l)


  if(length(prob_m) > 1){
      multiply<-plus_minus(split,phase_diff,prob_m[1])
    }else{
      multiply <- plus_minus(split,phase_diff,prob_m)
    }

for(i in 2:l){
  if((states[i]==1)&(current_state==0)){
      multiply <- plus_minus(split,phase_diff,prob_m)
    }
  if(states[i]==1){
    prob_m_mod[i] <-  multiply
  }
  current_state<-states[i]
}

return(prob_m_mod)

}

###############################################################################
#Randomly choose plus or minus                                                #
###############################################################################
plus_minus<-function(split,phase_diff,prob_m){

  if(runif(1, 0, 1) >= split){
    mod_phase_diff <- phase_diff

    multiply <- prob_m+mod_phase_diff

    if(multiply > 1){
      multiply <- prob_m-mod_phase_diff
    }
  }else{
            

    mod_phase_diff <- phase_diff
    multiply <- prob_m-mod_phase_diff
    if(multiply < 0){
      multiply <- prob_m+mod_phase_diff
    }
  }

return(multiply)
}

###############################################################################
#Create file from real data set                                               #
###############################################################################
create_real_data<-function(filename,number_of_samples,number_of_replicas){
#set up the variables
delta <- c(0, 1)
#create a HMM model that will simulate the gaps between CpGs. it is a 2 state HMM
#in order to simulate CpG islands and CpG deserts. TODO: multistates for CpG shores
#,cliffs etc etc etc to be added later.
x <- dthmm(NULL, Pi, delta, "exp", list(rate=rates_for_HMM_for_CpG_locations))
#simulate this HMM for n steps
x <- simulate(x, nsim=n)
#extract the length
l<-length(x$x)
#initialise a matrix to store the results
locs <- matrix(data=1,nrow=1,ncol=l)
#loop over the simulated output and create the locations
for(i in 2:l){
  locs[i]<-round(locs[i-1]+x$x[i])+1
}
return(locs)
}

###############################################################################
#add a random fluctuation to probability of success                           #
###############################################################################
add_random_noise <- function (theta,n_mean,n_standard_deviation){

#create the delta value with defined mean and standard deviation
delta <- rnorm(1, mean = n_mean, sd = n_standard_deviation)
theta_bar <-0


if(theta<=0){
  theta<-0.0000000000001
}
if(theta>=1){
  theta<-0.9999999999999
}
#convert theta to theta_bar based on the delta value calculated above.
theta_bar <-exp(log(theta/(1-theta))+delta)/(1+exp(log(theta/(1-theta))+delta))

#convert 0 and 1s
if(theta_bar>=1){theta_bar<-0.9999999999999}
if(theta_bar<=0){thata_bar<-0.0000000000001}
return(theta_bar)
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
init_interactive <- function (){
#get the variables from the command line
cat("Leave blank to use the default value:\n")
number_of_CpGs <-get_variable("Please enter the number of CpGs to simulate (default: 1000)",1000)

probability_of_success_in_differentially_methylated_region<-get_variable("Please enter probability of success in  methylated region (default: 0.9)",0.9)

probability_of_success_in_non_methylated_region<-get_variable("Please enter probability of success in non methylated region (default: 0.1)",0.1)

error_rate_in_differentially_methylated_region <-get_variable("Please enter error rate in differentially methylated region (default: 0.1)",0.1)

error_rate_in_non_differentially_methylated_region <-get_variable("Please enter error rate in non differentially methylated region (default: 0.1)",0.1)

mean_number_of_reads_in_differentially_methylated_region <-get_variable("Please enter the mean number of reads in non differentially methylated region (default: 30)",30)

mean_number_of_reads_in_non_differentially_methylated_region <-get_variable("Please enter the mean number of reads in differentially methylated region (default: 30)",30)

number_of_replicas <-get_variable("Please enter the number of replicas (default: 3)",3)

number_of_samples <-2

phase_diff <-get_variable("Please enter the phase diff (default: 0.05)",0.05)

balance <-get_variable("Please enter the balance (default: 0.5)",0.5)

cpg_matrix <-convert_to_array(get_text("Please enter the cpg matrix for the transitions (default: 0.05,3.480543e-04,0.05,0.025)","0.05,3.480543e-04,0.05,0.025"))

analysis_folder <-get_text("Please enter the analysis folder you wish to store the results (default: /tmp/)","/tmp/")

analysis_name <-get_text("Please enter a name you wish to store the resultsunder (default: example_WGBS_run)","example_WGBS_run")

output_path <- paste0(analysis_folder,analysis_name)

type <-get_text("Do you want to use the binomial or truncated simulator (default: binomial)","binomial")

nonconvPar <-convert_to_array(get_text("Please enter shape parameters for beta-distribution on nonconversion rates (default: NA)","NA"))

tranMPar <-convert_to_array(get_text("Please enter transition probabilities and length rate moderation for methylation regions (default: 0.01,0.001895)","0.01,0.001895"))

tranNPar <-convert_to_array(get_text("Please enter transition probabilities and length rate moderation for non-methylated regions (default: 0.08,0.001895)","0.08,0.001895"))

phase_difference_between_two_samples <- c(0,phase_diff)

file_to_write_results_to <- '/tmp/simulated_WGBS_data'
probs <- c(1,1,0.9,0.8,0.7,0.6)

transition_size <- 0
#####setting that should not be set by the user
number_of_repeats <- 1
transition_matrix_for_CpG_locations <-matrix(c(0.65, 0.35,0.2, 0.8),byrow=TRUE, nrow=2)
transition_matrix_for_probability_distributions <- matrix(c(0.9,0.1,0.1,0.9),byrow=TRUE, nrow=2)
type_of_locations <- 2

cat(paste(c(number_of_CpGs,probability_of_success_in_differentially_methylated_region,probability_of_success_in_non_methylated_region,error_rate_in_differentially_methylated_region,error_rate_in_non_differentially_methylated_region,mean_number_of_reads_in_differentially_methylated_region,mean_number_of_reads_in_non_differentially_methylated_region,number_of_replicas,number_of_samples,phase_diff,balance,paste(as.numeric(cpg_matrix), collapse = ","),output_path,type,nonconvPar), collapse = "\n"))

#do the simulation
run <- init_simulation(number_of_CpGs,probability_of_success_in_differentially_methylated_region,probability_of_success_in_non_methylated_region,error_rate_in_differentially_methylated_region,error_rate_in_non_differentially_methylated_region,mean_number_of_reads_in_differentially_methylated_region,mean_number_of_reads_in_non_differentially_methylated_region,number_of_replicas,number_of_samples,phase_diff,balance,cpg_matrix,output_path,type,tranMPar, tranNPar, nonconvPar)
return(list(run,number_of_replicas,output_path))

}

###############################################################################
#initiate the function based on args                                          #
###############################################################################
init_simulation <- function (number_of_CpGs,probability_of_success_in_differentially_methylated_region,probability_of_success_in_non_methylated_region,error_rate_in_differentially_methylated_region,error_rate_in_non_differentially_methylated_region,mean_number_of_reads_in_differentially_methylated_region,mean_number_of_reads_in_non_differentially_methylated_region,number_of_replicas,number_of_samples,phase_diff,balance,cpg_matrix,output_path,type,tranMPar, tranNPar, nonconvPar){
#get the variables from the command line

phase_difference_between_two_samples <- c(0,phase_diff)


#error_rate_in_differentially_methylated_region <- 0.1

file_to_write_results_to <- output_path
probs <- c(1,1,0.9,0.8,0.7,0.6)


transition_size <- 0
#####probably not set by the user
number_of_repeats <- 1
transition_matrix_for_CpG_locations <-matrix(c(0.65, 0.35,0.2, 0.8),byrow=TRUE, nrow=2)
transition_matrix_for_probability_distributions <- matrix(c(0.9,0.1,0.1,0.9),byrow=TRUE, nrow=2)
type_of_locations <- 2

#run the simulation script
simulated_data <-generate_sim_set(number_of_CpGs,
                number_of_repeats,
                transition_matrix_for_CpG_locations,
                transition_matrix_for_probability_distributions,
                mean_number_of_reads_in_differentially_methylated_region,
                mean_number_of_reads_in_non_methylated_region,
                probability_of_success_in_differentially_methylated_region,
                probability_of_success_in_non_differentially_methylated_region,
                error_rate_in_differentially_methylated_region,
                error_rate_in_non_differentially_methylated_region,
                phase_difference_between_two_samples,
                file_to_write_results_to,
                type_of_locations,
                number_of_replicas,
                number_of_samples,
                cpg_matrix,
                transition_size,
                balance,
                output_path,
                type, tranMPar, tranNPar, nonconvPar
                  )
#run the summary script
summarise_simulated_data(simulated_data,number_of_replicas,number_of_samples)

filename<-simulated_data[[3]]

return(filename)
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
#Function to summarise the  data from the simulation. This should show: 1)    #
#distribution of methlated/unmethylated sites 2) distribution of distances    #
#between CpGs 2) Distribution of the size of differentially methylated regions#
# 3) Distribution of read counts 4) 
###############################################################################
summarise_simulated_data <- function(generated_set,number_of_replicas,number_of_samples){
  #define the layout of the figures (actually has no effect when saving to pdf)
  mat <- matrix(1:6, 2, 3)
  layout(mat)
  #take the name from the simulated set and append the pdf extension
  outname <- generated_set[[3]]
  pdf(paste(outname,".pdf",sep=""))
  #create the graphs
  read_count_distribution(generated_set,number_of_replicas,number_of_samples)
  methylation_distribution(generated_set,number_of_replicas,number_of_samples)
  distance_distribution(generated_set)
  dmr_size(generated_set)
  methylation_postions(generated_set)
  dev.off()
}
###############################################################################
#Distribution of proportion of methylation                                    #
###############################################################################
methylation_distribution <- function(generated_set,number_of_replicas,number_of_samples){
  #calcualte the proportion of methylated vs demethylated reads in each replicate
  #and plot the valus as a histogram
  
  for(i in 1:number_of_samples){
    for(j in 1:number_of_replicas){
      col <- (((i-1)*(number_of_replicas*4))+((j-1)*4)+1)+2
      hist(generated_set[[2]][,col]/generated_set[[2]][,col+1],col=rgb(0,0,0,0.2),ylab="frequency",xlab="proportion of methylated CpGs",main="histogram of methylation proportion")
    }
  }
}

###############################################################################
#CpG locations on a line                                                      #
###############################################################################
methylation_postions <- function(generated_set){
  #plot the location of each CpG as a transparent bar on a line
  plot(generated_set[[2]][,1],matrix(data=1,ncol=1,nrow=length(generated_set[[2]][,1])),pch='|',col=rgb(0,0,0,0.01),ylab="",xlab="location",yaxt='n', ann=FALSE)

}
###############################################################################
#Distribution of distance between sites                                       #
###############################################################################
distance_distribution <- function(generated_set){
  #get the number of CpGs
  CpG_length <- length(generated_set[[2]][,1])
  #initialise a matrix to store the lengths in
  length_set <- matrix(data=0,nrow=6,ncol=CpG_length)
  #loop through and save the lengths between each site
  for (i in 2:CpG_length){
    length_set[1,i]<-generated_set[[2]][i,1]-generated_set[[2]][i-1,1]
  }
  #plot the distribution
  hist(length_set[1,],col=rgb(1,0,0,0.2),breaks=100,ylab="frequency",xlab="coverage",main="histogram of distances between CpGs")
}
###############################################################################
#Distribution of size of differentially methylated regions                    #
###############################################################################
dmr_size <- function(generated_set){
  #get the number of CpGs
  CpG_states <- length(generated_set[[2]][,2])
  #initialise a matrix to store the 
  length_set <- matrix(data=NA,nrow=6,ncol=CpG_states)
  #initialise a variable to store the number of CpGs in each state
  counter <- 1
  #initialise a variable with the first state
  previous <- generated_set[[2]][1,2]
  #loop through all states and store the counter variable each time the state
  #changes.
  for (i in 2:CpG_states){
    if(previous == generated_set[[2]][i,2]){
      counter <- counter + 1
    }else{
      counter <- counter + 1
      length_set[i] <- counter
      counter <- 1
      previous <- generated_set[[2]][i,2]
    }
  }
  #plot the distribution
  hist(length_set,breaks=100,ylab="frequency",xlab="size of DMR",main="histogram of DMR sizes")
}
###############################################################################
#Distribution of read counts                                                  #
###############################################################################
read_count_distribution <- function(generated_set,number_of_replicas,number_of_samples){
  #plot a histogram for each read count set.
  for(i in 1:number_of_samples){
    for(j in 1:number_of_replicas){
      col <- ((i-1)*(number_of_replicas*4))+((j-1)*4)+2
      
      hist(as.numeric(generated_set[[2]][,col+2]),col=rgb(1,0,0,0.2),ylab="frequency",xlab="coverage",main=paste("histogram of simulation read counts ",col))

    }
  }
}
###############################################################################
#Simple function which reduces the probability of changing state based on the #
#distance between two CpGs.                                                   #
###############################################################################
reduce_by_distance <- function (x,b){
  #given a value and decay value decrease a value by the the decay function.
  y <- exp(-b * log(x))
  return(y)
}

###############################################################################
#Simple function which takes a comma seperated list and converts it into a    #
#square matrix for use as a transition matirx in the HMM package              #
###############################################################################
convert_to_matrix <- function (x){
  #split the string by commas
  list <- unlist(strsplit(x, split=","))
  size <- sqrt(length(list))
  mat <- matrix(list,byrow=TRUE, nrow=size)
  mode(mat) <- "numeric"
   return(mat)
}

###############################################################################
#Simple function which takes a comma seperated list and converts it into a    #
#array for use as a transition matirx in the HMM package              #
###############################################################################
convert_to_array <- function (x){
  #split the string by commas
  list <- unlist(strsplit(x, split=","))
  size <- sqrt(length(list))
  mode(list) <- "numeric"
   return(list)
}

###############################################################################
#This function simulates the transitions from non and differentially methyd   #
#states. It requires n = the number of CpG, Pi = probability of changing state#
#, locs = the postion of the CpGs and start which is the initial starting prob#
# TO DO:                                                                      #
#At the moment the reduction by space is using the hard coded parameter       #
#-3.895e-02 but this should be updated in the future and helper function to   #
#calculate it from the data should be written                                 #
###############################################################################
simulate_state_transition <- function (n,Pi,locs,start,transition_size,rates){
  #initially set the current state to 2
  trans_up <- seq(from = 3, to = transition_size+2)
  trans_down <- seq(from = transition_size+2, to = 3)
  current_state<-2
  #randomly reassign based the on the start probability
  if(runif(1,0,1)<start){
    current_state<-1
  }
  #declare the variables
  state_transition <- matrix(data=1,nrow=1,ncol=n)
  prob_transition <- matrix(data=0,nrow=1,ncol=n)
  dist_transition <- matrix(data=0,nrow=1,ncol=n)
  
  #set the first state as that calculated above
  state_transition[1]<-current_state
  
  #loop through all of the CpGs
  i <- 2
  issue <- 0;
  while(i < length(state_transition)){
    #calculate the distance from the previous CpG
    dist<-locs[i]-locs[i-1];
    #store the distance
    dist_transition[i]<-dist
    #calcualte the probability of changing state based on the distance (and given rates - TJH)
    prob<-Pi[current_state]*reduce_by_distance(dist,rates[current_state])
    #store the value
    prob_transition[i]<-prob
    #randomly sample to see if the state should change
    if(runif(1,0,1)<prob){
      #if no: then save the current 
      if(current_state < 4){
      current_state<-current_state+1
      state_transition[i]<-current_state
      }else{
      current_state<-1
      state_transition[i]<-current_state
      }
    }
    state_transition[i]<-current_state
    i <- i +1
  }
  #return this value
  return(state_transition)
}


###############################################################################
#This function takes the CpG locations and state transitions and produces the #
#simulated read counts at each location. The function requires: states = the #
#progression of states for each CpG, prob_m the probability of success in the #
#not differentially methylated state, prob_d = probability of success in the #
#differentially methylated state, mean_m = the mean in the non differentially #
#state, mean_d = the mean in the differentially methylated state, s = ?, error#
#_m/d is the error sd, locs is the location of the CpGs. #
###############################################################################
generate_replicat_methyl_oldbin_data <-function (states,prob_m_a,prob_d,mean_m,mean_d,s,error_m,error_d,locs){
#calculate the number of sites
  l<-length(states);
  #create a vector of coverage for the non_diffs
  coverages<-list()
  for (s in 1:4){
  coverages[[s]]<-rpois(l, means[[s]])+1
  }
  #create a vecotr of coverage for the diffs
  #set the coverage to the meth vector
  coverage <- coverages[[1]]
  #initialse the vectors
  all_data <- matrix(data=0,nrow=1,ncol=l)
  #prob_m_s = prob_m_a
  #prob_d = matrix(data=prob_d,nrow=l,ncol=1)
  #loop through and add random noise the probs of success
  for (s in 1:4){
    for (i in 1:l){
      probs[[j]][[s]][i]<- add_random_noise(probs[[j]][[s]][i],0,errors[[s]])
    }
  }
#loop through and add random noise the probs of success
for (i in 1:l){
prob_m_s[i]<- add_random_noise(prob_m_s[i],0,error_m)
}
for (i in 1:l){prob_d[i]<- add_random_noise(prob_d[i],0,error_d)}
#loop through and simulate the number of methylated reads at each site
for(i in 1:l){
if(states[i]==1){
c<-rbinom(1, coverage_m[i], prob_m_s[i])
#c<-rbinom(1, coverage_d[i], probs[states[i]])
all_data[i]<-c
coverage[i]<-coverage_m[i]
}else{
c<-rbinom(1, coverage_d[i], prob_d[i])
all_data[i]<-c
coverage[i]<-coverage_d[i]
}
}
pos <- locs
#create one variable to return
proportion <- all_data/coverage;
final_store<-rbind(all_data,coverage,coverage,proportion)
rownames(final_store)<-c("x","n","coverage","proportion")
return(final_store)
}

generate_replicat_methyl_bin_data <-function (states,probs,means,s,errors,locs,j,output_path){
  #calculate the number of sites
  l<-length(states);
  #create a vector of coverage for the non_diffs
  coverages<-list()
  for (s in 1:4){
  coverages[[s]]<-rpois(l, means[[s]])+1
  #coverages[[s]]<-rep(100,l)
  }
  #create a vecotr of coverage for the diffs
  #set the coverage to the meth vector
  coverage <- coverages[[1]]
  #initialse the vectors
  all_data <- matrix(data=0,nrow=1,ncol=l)
  prs <- matrix(data=0,nrow=1,ncol=l)
  #prob_m_s = prob_m_a
  #prob_d = matrix(data=prob_d,nrow=l,ncol=1)
  #loop through and add random noise the probs of success
  for (s in 1:4){
    for (i in 1:l){
      probs[[j]][[s]][i]<- add_random_noise(probs[[j]][[s]][i],0,errors[[s]])
    }
    hist(probs[[j]][[s]],main=paste("histogram of probs in state",s))
  }
  mt <- list()
  for(s in 1:4){
    sj<-rbinom(l, coverages[[s]], probs[[j]][[s]])
    hist(sj,main=paste("histogram of sj in state",s))
    mt[[s]]<-sj
  }

  #loop through and simulate the number of methylated reads at each site
  for(i in 1:l){
      c<-mt[[states[i]]][i]
      all_data[i]<-c
      coverage[i]<-coverages[[states[i]]][i]
      prs[i]<-probs[[j]][[states[i]]][i]
  }
  pos <- locs
  #create one variable to return
  proportion <- all_data/coverage;
   hist(prs,main=paste("histogram of propotion gdgdfgf"))
  final_store<-rbind(all_data,coverage,coverage,proportion)
  plot(proportion[1:1000],main="methylation profile")
  rownames(final_store)<-c("x","n","coverage","proportion")
  return(final_store)
  
}



###############################################################################
#This function takes the CpG locations and state transitions and produces the #
#simulated read counts at each location. The function requires: states = the  #
#progression of states for each CpG, prob_m the probability of success in the #
#not differentially methylated state, prob_d = probability of success in the  #
#differentially methylated state, mean_m = the mean in the non differentially #
#state, mean_d = the mean in the differentially methylated state, s = ?, error#
#_m/d is the error sd, locs is the location of the CpGs.                      #
###############################################################################
generate_replicat_methyl_nbin_data_model_3 <-function (states,probs,means,s,errors,locs,j,size){
   #calculate the number of sites
  l<-length(states);
  #create a vector of coverage for the non_diffs
  coverages<-list()
  for (s in 1:4){
  coverages[[s]]<-rpois(l, means[[s]])+1
  }
  #create a vecotr of coverage for the diffs
  #set the coverage to the meth vector
  #initialse the vectors
  all_data <- matrix(data=0,nrow=1,ncol=l)
  coverage <- matrix(data=0,nrow=1,ncol=l)
  #prob_m_s = prob_m_a
  #prob_d = matrix(data=prob_d,nrow=l,ncol=1)
  #loop through and add random noise the probs of success
  for (s in 1:4){
    for (i in 1:l){
      probs[[j]][[s]][i]<- add_random_noise(probs[[j]][[s]][i],0,errors[[s]])
    }
  }
  mt <- list()
  for(s in 1:4){
    sj<-rbinom(l, coverages[[s]], probs[[j]][[s]])
    coverages[[s]]<-ceiling((sj+1) / probs[[j]][[s]])
    hist(sj,main=paste("histogram of sj in state",s))
    hist(coverages[[s]],main=paste("histogram of coverage in state",s))
    mt[[s]]<-sj
  }

  #loop through and simulate the number of methylated reads at each site
  for(i in 1:l){
      c<-mt[[states[i]]][i]
      all_data[i]<-c
      coverage[i]<-coverages[[states[i]]][i]
  }
  pos <- locs
  #create one variable to return
  proportion <- all_data/coverage;
  plot(proportion[1:1000],main="methylation profile")
  final_store<-rbind(all_data,coverage,coverage,proportion)
  rownames(final_store)<-c("x","n","coverage","proportion")

  return(final_store)
  
}

###############################################################################
#This function takes the CpG locations and state transitions and produces the #
#simulated read counts at each location. The function requires: states = the  #
#progression of states for each CpG, prob_m the probability of success in the #
#not differentially methylated state, prob_d = probability of success in the  #
#differentially methylated state, mean_m = the mean in the non differentially #
#state, mean_d = the mean in the differentially methylated state, s = ?, error#
#_m/d is the error sd, locs is the location of the CpGs.                      #
###############################################################################
generate_replicat_methyl_truncated_nbin_data <-function (states,probs,means,s,errors,locs,j,size,output_path){
   #calculate the number of sites
  l<-length(states);
  #create a vector of coverage for the non_diffs

  #create a vecotr of coverage for the diffs
  #set the coverage to the meth vector
  #initialse the vectors
  all_data <- matrix(data=0,nrow=1,ncol=l)
  coverage <- matrix(data=0,nrow=1,ncol=l)
  coverages<-list()
  prs <- matrix(data=0,nrow=1,ncol=l)
  #prob_m_s = prob_m_a
  #prob_d = matrix(data=prob_d,nrow=l,ncol=1)
  #loop through and add random noise the probs of success
  for (s in 1:4){
    for (i in 1:l){
      probs[[j]][[s]][i]<- add_random_noise(probs[[j]][[s]][i],0,errors[[s]])
    }
  }


  mt <- list()
  for(s in 1:4){
    sj <- matrix(data=0,nrow=1,ncol=l)
    coverages[[s]]<-matrix(data=0,nrow=1,ncol=l)
    for(i in 1:l){
    a<-means[[s]]*probs[[j]][[s]][i]*(probs[[j]][[s]][i]/(1-probs[[j]][[s]][i]))
    b<-(1-probs[[j]][[s]][i])/probs[[j]][[s]][i]
    N <- rpois(1,means[[s]])
    coverages[[s]][i]<-N
    repeat{
      lambda_j <- rgamma(1,shape = a,scale=b)
      x<-rpois(1,lambda_j)
      
      if(x<=N){
        sj[i]<-x
        
        break
      }

    }
    }
    mt[[s]]<-sj
  }

  #loop through and simulate the number of methylated reads at each site
  for(i in 1:l){
      c<-mt[[states[i]]][i]
      all_data[i]<-c
      coverage[i]<-coverages[[states[i]]][i]
      prs[i]<-probs[[j]][[states[i]]][i]
  }
  pos <- locs
  #create one variable to return
  proportion <- all_data/coverage;
  final_store<-rbind(all_data,coverage,coverage,proportion)
  rownames(final_store)<-c("x","n","coverage","proportion")
    plot(proportion[1:1000],main="methylation profile")
  save(prs,file=paste0(output_path,"_probabilities_",j,".RData"))
  return(final_store)
  
}



options(warn=1)
#init the values 

#run the interactive launch script if the correct number of command line arguments has not been supplied
args <- commandArgs(trailingOnly = TRUE)
if(length(args)>=14){
      filename<-init_simulation(as.numeric(args[1]),as.numeric(args[2]),as.numeric(args[3]),as.numeric(args[4]),as.numeric(args[5]),as.numeric(args[6]),as.numeric(args[7]),as.numeric(args[8]),as.numeric(args[9]),as.numeric(args[10]),as.numeric(args[11]),convert_to_array(args[12]),args[13],args[14],convert_to_array(args[15]),convert_to_array(args[16]),convert_to_array(args[17]))      cat("\n##########################SIMULATION COMPLETE############################\n")
      cat(paste("\nThe simulated data can be found here:\n",filename,"\n"))
      cat(paste0("\nA pdf summarising the simulated data can be found here:\n",filename,".pdf\n"))
      cmd<-paste("\nRscript benchmark_WGBS.R",filename,"2",as.numeric(args[8]),"1 999999 1 1 15 30 100 1000000 0 0 200 0 0",args[13],"0")
      cat(paste("\nTo run the benchmarking script on this dataset with the default settings use:\n",cmd,"\n"))
}else if(length(args)==1){
    if(args[1]=='interactive'){
       rd<-init_interactive()
      cat("\n##########################SIMULATION COMPLETE############################\n")
      cat(paste("\nThe simulated data can be found here:\n",rd[[1]],"\n"))
      cat(paste0("\nA pdf summarising the simulated data can be found here:\n",rd[[1]],".pdf\n"))
      cmd<-paste("\nRscript benchmark_WGBS.R",rd[[1]],"2",rd[[2]],"1 999999 1 1 15 30 100 1000000 0 0 200 0 0",rd[[3]],"0")
      cat(paste("\nTo run the benchmarking script on this dataset with the default settings use:\n",cmd,"\n"))
    }
}else{
  print("else")
  if(args[1]=='multi'){
      source("benchmark_WGBS.R")
      reps <- as.numeric(args[16])
      files<-list()
      rocs<-list()
      times<-list()
      for(rep in 1:reps){
        files[[rep]]<-init_simulation(as.numeric(args[2]),as.numeric(args[3]),as.numeric(args[4]),as.numeric(args[5]),as.numeric(args[6]),as.numeric(args[7]),as.numeric(args[8]),as.numeric(args[9]),as.numeric(args[10]),as.numeric(args[11]),as.numeric(args[12]),convert_to_array(args[13]),paste0(args[14],"rep_",rep),args[15],convert_to_array(args[16]),convert_to_array(args[17]),convert_to_array(args[18]))        index<-length(rocs)+1
        dat<-auto_init_benchmark(files[[rep]],2,3,1,999999,1,0.1,15,100,100,1000000,0.5,1,200,0.0,0,files[[rep]],0,'binomial')
        rocs[[rep]]<-dat[[1]]
        times[[rep]]<-dat[[2]]
        save(rocs,times,file=paste0(args[14],".ROCs.RData"))
      }
      lreps<-length(rocs)
      lcuts<-length(rocs[[1]])
      ltimes<-length(times[[1]])
      average_times<-list()
      for(t in 1:ltimes){
         average_times[[t]] = 0
        for(rep in 1:lreps){
        average_times[[t]]<-average_times[[t]]+times[[rep]][[t]]
        }
        average_times[[t]]<-average_times[[t]]/lreps
      }
      dim_rocs<-dim(rocs[[1]][[1]])
      average<-list()
      for (cut in 1:lcuts){
      average[[cut]]<-matrix(data=0,nrow=6,ncol=dim_rocs[2])
        for(rep in 1:lreps){
          for (j in 1:6){
            for (k in 1:dim_rocs[2]){
                average[[cut]][j,k]<-average[[cut]][j,k]+rocs[[rep]][[cut]][j,k][[1]]
            }
          }        
        }
        average[[cut]]<-average[[cut]]/lreps}
      }
      at<-array(unlist(average_times), dim = c(nrow(average_times[[1]]), ncol(average_times[[1]]), length(average_times)))
      save(average,at,file=paste0(args[14],".ROCs.Average.RData"))
      plot_ROC(average,c('methsig','bsmooth','fishers exact','MethyKit'),paste0(args[14],"average",".ROC.pdf"))
      plot_AUC(average,c('methsig','bsmooth','fishers exact','MethyKit'),paste0(args[14],"average",".AUC.pdf"))
      plot_Time(at,c('methylKit','bsmooth','fishers exact','MethyKit'),paste0(args[14],"average",".times.pdf"))

    }

###################################################################################################
#IF YOU WANT TO RUN THIS FROM WITHIN R THEN SOURCE THIS FILE AND THEN USE THE FOLLOWING CALL      #
#simulated_data <-generate_sim_set(...)
#                  
#
#summarise_simulated_data(simulated_data)
####################################################################################################
