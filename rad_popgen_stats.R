# This function takes as input:
# rad_vcf = vcfR object containing radseq data
# unseq_gff = This is a gff file that includes information for unsequenced regions, important for averaging across sequenced sites
# popmap = 2 column dataframe with vcf inds in column 1 and pop assignment in column 2
# min_snp = The minimum number of SNPs per rad locus in order to calculate statistics for it

# Function returns:
# fst, dxy and pi (or whichever subset specified), per rad-locus

# # Usage example, can choose any combo of stats
# popgen_res <- rad_popgen_stats(test_vcf,test_gff,stats_to_calc=c("dxy","pi","fst"),popmap,n_cores = 4,min_snp = 0)
# popgen_res$pi
# popgen_res$dxy
# popgen_res$fst

# Use these
library(vcfR)
library(parallel)
library(data.table)
library(hierfstat)
library(adegenet)

# Outline the function
rad_popgen_stats <- function(rad_vcf,unseq_gff,stats_to_calc=c("fst","dxy","pi"),popmap,n_cores=1,min_snp=0){
  
  # Fetch contig identifiers...
  contigs <- sapply(strsplit(test_vcf@fix[,3],":"),"[[",1)
  contig_labs <- unique(contigs)
  
  # Fetch contig specific gffs
  contig_gffs <- lapply(contig_labs,function(contig){
    unseq_gff[grep(paste0(contig,"_pos"),unseq_gff$V9),]
  })
  names(contig_gffs) <- contig_labs
  
  # Set our population labs
  pops <- unique(unlist(popmap[,2]))
  
  # Work through contigs
  contigs_out <- mclapply(contig_labs,function(contig){
    message(paste0("Processing contig ",contig))
    # Get contig size and unsequenced size
    tmp_gff <- data.frame(contig_gffs[as.character(contig)])
    seq_size <- length(unlist(sapply(1:nrow(tmp_gff),function(x){
      return(seq(tmp_gff[x,4],tmp_gff[x,5]))
    })))
    
    # Subset VCF and for individuals in popmap
    vcf_sub <- rad_vcf[which(contigs == contig),]
    vcf_sub <- vcf_sub[,c(1,which(colnames(vcf_sub@gt) %in% popmap[,1]))]
    
    # Skip over these if you want
    if(nrow(vcf_sub@fix) < min_snp){
      return(NULL)
    } else {
      
      # Get allele counts
      gt <- extract.gt(vcf_sub)
      
      # Convert to numbers
      gt[gt == "0/0"] <- 0
      gt[gt == "0/1"] <- 0.5
      gt[gt == "1/1"] <- 1
      suppressWarnings(class(gt) <- "numeric")
      
      # Filter the genotypes based on popmap
      gt <- gt[,colnames(gt) %in% popmap[,1]]
      
      # Coerce if needed
      if(is.null(nrow(gt))){
        gt <- t(as.matrix(gt))
      }
      
      # First we'll get our count_dd
      if("pi" %in% stats_to_calc | "dxy" %in% stats_to_calc){
        pop_counts <- lapply(1:length(pops),function(pop){
          
          # Take column
          tmp_gt <- gt[,popmap[popmap[,2] == pops[pop],1]]
          
          # Coerce if needed
          if(is.null(nrow(tmp_gt))){
            tmp_gt <- t(as.matrix(tmp_gt))
          }
          
          # Merge to data.frame
          count_dd <- data.frame(count1 = rowSums(tmp_gt,na.rm = T),
                                 count2 = NA)
          count_dd$count2 <- ncol(tmp_gt) - count_dd$count1 - rowSums(is.na(tmp_gt))
          
          # Sum
          count_dd$sum <- rowSums(count_dd)
          
          return(count_dd)
        })
        names(pop_counts) <- pops
      }
      
      # Use these to get pi
      if("pi" %in% stats_to_calc){
        pop_pis <- sapply(pop_counts,function(counts){
          pi_tmp <-  (2*(counts$count1 * counts$count2))/(counts$sum*(counts$sum-1))
          return(sum(pi_tmp)/seq_size)
        })
        pop_pis <- data.frame(pop_pis)
      } else {
        pop_pis = NULL
      }
      
      # And use to get dxy
      if("dxy" %in% stats_to_calc){
        dxy_mat <- matrix(ncol=length(pops),nrow=length(pops))
        rownames(dxy_mat) <- pops
        colnames(dxy_mat) <- pops
        comps_to_make <- combn(1:length(pops),2)
        for(i in 1:ncol(comps_to_make)){
          pop1 = pops[comps_to_make[1,i]]
          pop2 = pops[comps_to_make[2,i]]
          
          pop1_count <- data.frame(pop_counts[pop1])
          pop2_count <- data.frame(pop_counts[pop2])
          
          dxy_tmp <- 1-(
            ((pop1_count[,1] * pop2_count[,1])/(pop1_count[,3] * pop2_count[,3]))+
              ((pop1_count[,2] * pop2_count[,2])/(pop1_count[,3] * pop2_count[,3])))
          
          dxy_mat[pop1,pop2] <- sum(dxy_tmp)/seq_size
        }
        
        # Convert for output
        dxy_tmp <- reshape2::melt(dxy_mat)
        dxy_tmp$contig <- contig
      } else {
        dxy_tmp <- NULL
      }
      
      # hierfstat for Fst
      if("fst" %in% stats_to_calc){
        fst_mat <- matrix(ncol=length(pops),nrow=length(pops))
        rownames(fst_mat) <- pops
        colnames(fst_mat) <- pops
        comps_to_make <- combn(1:length(pops),2)
        for(i in 1:ncol(comps_to_make)){
          pop1 = pops[comps_to_make[1,i]]
          pop2 = pops[comps_to_make[2,i]]
          
          # Get individuals to keep
          inds_to_keep <- popmap[popmap[,2] %in% c(pop1,pop2),1]
          
          # Calculate
          dat <- suppressWarnings(vcfR2genind(vcf_sub[,c(1,which(colnames(vcf_sub@gt) %in% inds_to_keep))]))
          pop(dat) <- popmap[popmap[,1] %in% rownames(dat$tab),2]
          dat2 <- genind2hierfstat(dat)
          stats <- basic.stats(dat2,diploid = 2,digits = 2)
          fst_mat[pop1,pop2] <- stats$overall["Fst"]
        }
        
        # Convert for output
        fst_tmp <- reshape2::melt(fst_mat)
        fst_tmp$contig <- contig
      } else {
        fst_tmp <- NULL
      }
      
      # Return it all
      return(list(pi=pop_pis,
                  dxy=dxy_tmp,
                  fst=fst_tmp))
    } 
  },mc.cores=n_cores)
  
  # Return single lists for each type
  pi_out <- data.frame(rbindlist(lapply(contigs_out,"[[",1)))
  if("pi" %in% stats_to_calc){
    pi_out$contig <- rep(contig_labs,each=length(pops))
    pi_out$pop <- rep(pops)
  }
  dxy_out <- na.omit(data.frame(rbindlist(lapply(contigs_out,"[[",2))))
  if("dxy" %in% stats_to_calc){
    colnames(dxy_out) <- c("pop1","pop2","dxy","contig")
  }
  
  fst_out <- na.omit(data.frame(rbindlist(lapply(contigs_out,"[[",3))))
  if("fst" %in% stats_to_calc){
    colnames(fst_out) <- c("pop1","pop2","fst","contig")
  }
  
  return(list(pi=pi_out,
              dxy=dxy_out,
              fst=fst_out))
}
