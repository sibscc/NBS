#######################################################################
#
#	Transforming VCF files or gVCF file to gene by sample matrix
#
#	Author:	Ma Chengcheng
#
#	Date:	2017-3-2
#
#######################################################################
library(getopt)
spec <- matrix(c(
  'verbose', 'v', 2, "integer", "Verbose option",
  'input', 'i', 1, "character", "input vcf path or gvcf file[must]",
  'network', 'n', 1,"character", "Network file with genes or Entrez IDs to be use[must]",
  'output', 'o', 1, "character", "Output matrix file[must]",
  'gvcf',   'g', 2, "logical", "whether gvcf file to use",
  'germline', 'm', 2, "logical", "need germline mutation (0 1 2) or not"
), byrow=T, ncol=5)


opts <- getopt(spec)

if(is.null(opts$input) | is.null(opts$network) |is.null(opts$output)){
  cat(getopt(spec, usage=T))
  q(status =1 )
}

library(org.Hs.eg.db)
library(data.table)
input <- opts$input
network <- opts$network
output <- opts$output
if(is.null(opts$gvcf)){gvcf = F}else{gvcf = opts$gvcf}
if(is.null(opts$germline)){germline = F}else{germline = opts$germline}
print(germline)

# input <- "~/workdir/Data/NBS/VCFs/"
# network <- "~/workdir/Data/NBS/network/HumanIntMatrix"
# out <- "~/workdir/Data/NBS/Results/Vcf2Matrix"

if(!file.exists(input)){
  pirnt("Input path/file not exist!!")
  q(status=1)
}

#### Load EntrezID to gene
x <- org.Hs.egSYMBOL2EG
xx <- mappedkeys(x)
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
####### get the nonsynonymous mutated genes in a vcf file
cat("Data loaded. Start ...!!!\n")
if(gvcf == F){
	samples <- dir(input, pattern="*.vcf")
	IDs <- fread(network, select=1)[[1]]    ## read the network matrix to get the rownames of VCF matrix
	out <- matrix(0, length(IDs), length(files))
	colnames(out) <- files
	rownames(out) <- IDs
	if(germline == F){
	    for(i in 1:length(samples)){
		  f = read.table(file.path(input, files[i]))
		  non <- grep("nonsynonymous",f[,8],value=T)
		  nn <- unlist(strsplit(non, ";"))
		  refgene <- grep("Gene.refGene=", nn, value=T)
		  refgene <- unique(refgene)
		  gene <- unlist(strsplit(refgene, "="))[seq(2, 2*length(refgene), 2)]
		  ############## Search in org.Hs.eg.fb
		  #gene <- unlist(strsplit(gene, ","))
		  #gene <- unlist(strsplit(gene, "\\\\x3b"))
		  gene <- unique(gene)
		  yy <- unlist(xx[gene])
		  zz <- intersect(rownames(out), yy)
		  out[zz, files[i]] <- 1
		}
	}else{
	    for(i in 1:length(samples)){
	          f = read.table(file.path(input, files[i]))
	          non <- grep("nonsynonymous",f[,8],value=T)
		  genelist <- numeric()
	          for(j in 1:length(non)){
			item <- strsplit(non[j], ";")[[1]]
			AC <- strsplit(grep("^AC=",item, value=T), "=")[[1]][2]
			gene <-  strsplit(grep("^Gene.refGene=",item, value=T), "=")[[1]][2]
			if(AC ==1){
				if(!(gene %in% names(genelist))){
					genelist[gene] <- 1
				}
			}else if(AC == 2){
				genelist[gene] <- 2
			}
		  }
		  Symbo2Id <- xx[names(genelist)]
		  Symbo2Id[sapply(Symbo2Id, is.null)] <- NA
		  EZs <- unlist(Symbo2Id)
		  names(genelist) <- EZs
	          ############## Search in org.Hs.eg.fb
	          zz <- intersect(rownames(out), EZs)
	          out[zz, files[i]] <- genelist[zz]
	        }
	}
	write.table(out, file=output,sep="\t", col.names=NA)
}else{
    files <- read.table(input, header=T,fill=T,sep="\t")
    start = which(colnames(files) == "FORMAT")
    samples <- colnames(files)[(start + 1) : ncol(files)]
    IDs <- fread(network, select=1)[[1]]    ## read the network matrix to get the rownames of VCF matrix
    out <- matrix(0, length(IDs), length(samples))
    colnames(out) <- samples
    rownames(out) <- IDs
    splicing_i <-grep("splicing", files[,"Func.refGene"])     # use genes with splicing change or nonsynonymous change
        nonsy_i <- grep("nonsyn", files[,"ExonicFunc.refGene"])
        stop_i <- grep("stopgain", files[,"ExonicFunc.refGene"])
        frame_i <- grep("frameshift", files[,"ExonicFunc.refGene"])
        change_i <- unique(c(splicing_i, nonsy_i, stop_i, frame_i))
        change_m <- files[change_i, ]
        genes <- unique(change_m[,"Gene.refGene"])
        geneIndex <- list()
	for(i in as.character(genes)){
                geneIndex[[i]] <- which(i == change_m[,"Gene.refGene"])
        }
    if(germline == F){						  ## gVCF somatic mutation or not consider homo or heter mutation
	for(i in samples){
		gene <- character()
		for(j in names(geneIndex)){
			gt <- change_m[geneIndex[[j]], i]
			homo <- grep("0/0", gt)
			unk <- grep("\\./\\.", gt)
			if(length(homo) + length(unk) < length(gt)){
				gene <- c(gene, j)
			}
		}
		yy <- unlist(xx[gene])
                zz <- intersect(rownames(out), yy)
                out[zz, i] <- 1
	}
    }else{
	genelist <- numeric()
		

    }    
}


cat("VCF to Matrix Done!!!\n")
















