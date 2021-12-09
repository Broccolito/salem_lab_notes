## Converting dosage files to vcf files

### Need:

- Dosage file

```
5 5:34973992:C:G 5:34973992:C:G 34973992 C G 0 1 1.97501 0.99699 0.97099 0.99899 0.978 0 0 1.99699 1 0.987 1 0.99899 1.99699 0.998 0.998 0.998 0.994 1.963 0 0.99899 0.998 0.998 0.953 0.98202 0.978 0 0 0 0 1.955 1.45398 0.974 0.979 0.966 0.822 1.88501 0.98199 1.46799 0.99699 1.86698 0.989 0.97099 0.99899 0 0 0.953 0.97301 0 0 0.998 0.957 0.998 1 1.957 0.96201 0 0.98 0 0.965 0 0 0 0.989 1 0 0.998 0.989 1 0 0 0 1 0 1 1 0 1 0.99699 0.934 0.92 0 
...
```



- Sample file

```
ID_1 ID_2 missing
0 0 0
1_3226 1_3226 0
1_8228 1_8228 0
2_2294 2_2294 0
2_3416 2_3416 0
2_6790 2_6790 0
2_6952 2_6952 0
...
```



- VCF file name

```
xxx_unphased.vcf
```



### Use the R function

```R
library(dplyr)
library(R.utils)
library(vcfR)

dosage_to_vcf = function(dosage_dir = "FHS_EA_TOPHIT_SNPs_subset.dosage",
                         sample_dir = "FHS_EA_TOPHIT_SNPs_subset.sample",
                         vcf_filename = "FHS_EA_TOPHIT_SNPs_subset_unphased.vcf",
                         gzip_compression = TRUE){
  
  # Read dosage files and sample file
  dosage = read.delim(dosage_dir, sep = " ", header = FALSE)
  sample_info = read.delim(sample_dir, sep = " ", header = FALSE)
  samples = paste0("SUBJECT", sample_info[,1][-(1:2)])
  
  # Get rid of redundant columns for the dosage file
  dosage = dosage[,-c(2,3)]
  names(dosage) = c("chr", "pos", "ref", "alt", samples)
  
  # Solidify genotype calls
  dosage_called = dosage %>%
    mutate_at(vars(matches("SUBJECT")), round) %>%
    mutate_at(vars(matches("SUBJECT")), function(x){
      ifelse(is.na(x),"./.",
             ifelse(x==0, "0/0", ifelse(x==1, "0/1", ifelse(x==2, "1/1","./.")))
      )
    })
  
  
  line1 = "##fileformat=VCFv4.2"
  line2 = "##FORMAT=<ID=GT,Number=1,Type=Integer,Description=Genotype>"
  line3 = paste("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
                paste(samples, collapse = "\t"), sep = "\t")
  
  line_gt = vector()
  for(i in 1:dim(dosage_called)[1]){
    line_meta = paste(
      dosage_called[i,1], # Chromosome
      dosage_called[i,2], # Position
      paste0("chr", dosage_called[i,1],
             "_", dosage_called[i,2]), # ID
      dosage_called[i,3], #Ref
      dosage_called[i,4], # Alt
      ".", "PASS", ".", "GT", # Quality, filter, info, and format
      sep = "\t")
    line_x = paste(line_meta, 
                   paste(dosage_called[i,5:dim(dosage_called)[2]],
                         collapse = "\t"),
                   sep = "\t")
    line_gt = c(line_gt, line_x)
  }
  
  writeLines(c(line1, line2, line3, line_gt),
             con = vcf_filename,
             sep = "\n")
  if(gzip_compression){
    gzip(vcf_filename)
  }
  
}

### NOT RUN
# dosage_to_vcf(dosage_dir = "FHS_EA_TOPHIT_SNPs_subset.dosage",
#               sample_dir = "FHS_EA_TOPHIT_SNPs_subset.sample",
#               vcf_filename = "FHS_EA_TOPHIT_SNPs_subset_unphased.vcf",
#               gzip_compression = TRUE)
# vcf = read.vcfR("FHS_EA_TOPHIT_SNPs_subset_unphased.vcf.gz")
```



### Generate haplotype frequency from phased VCF files

```R
library(dplyr)
library(vcfR)
library(gtools)

vcf_to_haplotype = function(phased_vcf = "FHS_EA_TOPHIT_SNPs_subset_phased.vcf.gz",
                            filename = "FHS_EA_TOPHIT_SNPs_subset_haplotype.csv"){
  
  vcf = read.vcfR(phased_vcf)
  
  n_snps = dim(vcf@fix)[1]
  possible_haplotypes = apply(permutations(n=2,r=n_snps,repeats.allowed=TRUE)-1,1,
                              FUN = function(x){paste(x, collapse = "")})
  
  all_haplotypes = vector()
  for(hap in possible_haplotypes){
    gt_df = as.data.frame(vcf@gt[,-1])
    haplotypes = vector()
    for(gt in gt_df){
      paternal = gt %>%
        strsplit(split = "|") %>%
        lapply(function(x){x[1]}) %>%
        unlist() %>%
        paste(collapse = "")
      maternal = gt %>%
        strsplit(split = "|") %>%
        lapply(function(x){x[3]}) %>%
        unlist() %>%
        paste(collapse = "")
      haplotype = sum(c(paternal==hap, maternal==hap))
      haplotypes = c(haplotypes, haplotype)
    }
    all_haplotypes = rbind(all_haplotypes, haplotypes)
  }
  all_haplotypes = as.data.frame(t(all_haplotypes))
  
  to_nucleotide = function(allele = "000"){
    allele = unlist(strsplit(allele, split = ""))
    for(i in 1:length(allele)){
      allele[i] = ifelse(allele[i]=="0",getREF(vcf)[i],
                         ifelse(allele[i]=="1",
                                getALT(vcf)[i],"."))
    }
    allele = paste(allele,collapse = "")
    return(allele)
  }
  
  names(all_haplotypes) = paste0("HAP_", sapply(possible_haplotypes, to_nucleotide))
  all_haplotypes = cbind.data.frame(tibble(SUBJECT = colnames(vcf@gt)[-1]),
                                    all_haplotypes)
  write.csv(all_haplotypes, file = filename, quote = FALSE, row.names = FALSE)
  
}

###NOT RUN
# vcf_to_haplotype(phased_vcf = "FHS_EA_TOPHIT_SNPs_subset_phased.vcf.gz",
#                  filename = "FHS_EA_TOPHIT_SNPs_subset_haplotype.csv")

```

