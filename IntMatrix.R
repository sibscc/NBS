########################################################################
#
#	Making gene interaction matrix from data like:
#	gene1	gene2	interaction_score
#
#	Author: Ma Chengcheng
#
#	Date: 2017-3-1
#######################################################################
library(getopt)
spec <- matrix(c(
  'verbose', 'v', 2, "integer", "Verbose option",
  'input_network', 'i', 1,"character", "Network file with genes in the first two columns and score in the third column",
  'output', 'o', 1, "character", "Output matrix file"
  ), byrow=T, ncol=5)

opts <- getopt(spec)

if(is.null(opts$input_network) |is.null(opts$output)){
        cat(getopt(spec, usage=T))
        q(status =1 )
}

net <- opts$input_network
output <- opts$output
network <- read.table(net)

makematrix <- function(network){
  id <- as.character(unique(c(network[,1], network[,2])))
  int <- matrix(,length(id),length(id))
  colnames(int) <- rownames(int) <- id
  for(i in 1:nrow(network)){
    if(i %% 10000==0){print(i)}
    #print(i)
    int[as.character(network[i,1]), as.character(network[i,2])] <- round(network[i,3],3)

    int[as.character(network[i,2]), as.character(network[i,1])] <- round(network[i,3],3)
  }
  int[which(is.na(int))] <- 0
  return(int)
}

cat("Transforming interactions to matrix...\n")
int_net <- makematrix(network)
cat("Transforming done! Writing interaction matrix...\n")
write.table(int_net, file=output)
cat("Interaction to matrix Done!!\n")
