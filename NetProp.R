####################################################################
#
#	Propagate the mutation through the network
#
#	Author:	Ma Chengchegn
#
#	Date:	2017-3-1
#
####################################################################


library(getopt)

spec <- matrix(c(
  'verbose', 'v', 2, "integer", "Verbose option",  
  'input', 'i', 1, "character", "input data with 0 indicating normal, i indicating mutation",
  'network', 'n', 1,"character", "interaction network matrix gene by gene",
  'step', 's', 2, "numeric", "Step length for each iteration of propagation[deault 0.5]",
  'contract', 'c', 2, "numeric", "contraction limit number[default 1e-6]",
  'output', 'o', 1, "character", "Output matrix file"
  ), byrow=T, ncol=5)

opts <- getopt(spec)

if(is.null(opts$input) | is.null(opts$network) |is.null(opts$output)){
	cat(getopt(spec, usage=T))
	q(status =1 )
}
#print(spec)
#print(opts)
if ( is.null(opts$step) ) { opts$step= 0.5 }
if ( is.null(opts$contract) ) { opts$contract= 1e-6}
input <- opts$input
network <- opts$network
step <- opts$step
contract <- opts$contract
output <- opts$output

cat("Reading Input data...\n")
if(!file.exists(input)){
  print("Input file does not exist. Please check!")
  exit()
}else{
  people <- as.matrix(read.table(input,as.is=F, header=T))
  if(!is.numeric(people)){
    print("Non numeric input data. Please Chech!")
    q(status=1)
  }
}
cat("Reading network matrix...\n")
if(!file.exists(network)){
  print("Network file does not exist. Pliease check!")
  exit()
}else{
  network <- read.table(network,as.is=F)
  network <- as.matrix(network)
}
########################################## 
#   Functions
##########################################

###### Make network matrix

#makematrix <- function(network){
#  id <- as.character(unique(c(network[,1], network[,2])))
#  int <- matrix(,length(id),length(id))
#  colnames(int) <- rownames(int) <- id
#  for(i in 1:nrow(network)){
#    #if(i % 10000==0){print(i)}
#    print(i)
#    int[as.character(network[i,1]), as.character(network[i,2])] <- network[i,3]
#
#    int[as.character(network[i,2]), as.character(network[i,1])] <- network[i,3]
#  }
#  int[which(is.na(int))] <- 0
#  return(int)
#}

##### Propagation function
propagate <- function(F0, A, s){
  Fs <- Ft<- F0 <- t(F0)
  count = 1
  while(TRUE){
    Fs <- Ft
    Ft <- s*Ft%*%A + (1-s)*F0
    con <- norm(Ft) - norm(Fs)
    print(paste(count, con))
    if(abs(con) <= contract){
      break()
    }
    count = count+1
  }
  return(t(Ft))
}

##### Excecution
#int <- makematrix(network)
cat("Normalizing the interaction matrix...\n")
int_norm <- network/norm(network)   ### interaction matrix normalization by the matrix norm
cat("Network propagation...\n")
PropNet <- propagate(people, int_norm,0.5)
cat("Writing the propagated data...\n")
write.table(PropNet, file=output,row.names=F)
cat("Network propagation done!\n")




