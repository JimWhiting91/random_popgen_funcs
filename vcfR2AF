# Function calculates allele frequencies from an input VCF for all pops in the popmap
vcfR2AF <- function(vcf,popmap,n_cores=1,allele="REF"){
  library(parallel)
  
  # vcf should be a vcfR object
  # popmap should be a data.frame with 2 columns, column 1 for ind IDs and column 2 for pop ID
  # returns a matrix of allele frequencies for each pop in popmap
  
  # First calculate all population allele frequencies
  # Extract gt manually from VCF
  gt <- extract.gt(vcf)
  
  # Convert to numbers
  gt[gt == "0/0"] <- 0
  gt[gt == "0/1"] <- 0.5
  gt[gt == "1/1"] <- 1
  class(gt) <- "numeric"
  
  # Filter the genotypes based on popmap
  gt <- gt[,colnames(gt) %in% popmap[,1]]
  
  # For all populations, calculate the allele frequency of the ALT allele (it doesn't matter which we use)
  pops <- unique(unlist(popmap[,2]))
  
  pop_AF <- mclapply(pops,function(pop){
    
    # Get only these pops...
    gt_tmp <- gt[,popmap[popmap[,2]==pop,1]]
    
    # Get the AFs
    allele_counts <- rowSums(gt_tmp,na.rm = T)
    
    # Get counts, and factor in missingness...
    total_counts <- ncol(gt_tmp) - apply(gt_tmp, 1, function(x) sum(is.na(x)))
    
    return(allele_counts/total_counts)
    
  },mc.cores = n_cores)
  
  # Make a new matrix
  AF_mat <- matrix(ncol=length(pops),nrow=nrow(gt))
  for(i in 1:ncol(AF_mat)){
    AF_mat[,i] <- pop_AF[[i]]
  }
  colnames(AF_mat) <- pops
  
  # Get REF if needed
  if(allele=="REF"){
    AF_mat <- 1 - AF_mat
  }
  return(AF_mat)
}
