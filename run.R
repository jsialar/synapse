library(snpStats)
setwd('d:/synapse')

pathM <- paste("ROSMAP_arrayGenotype", c(".bed", ".bim", ".fam"), sep = "")
SNP_M <- read.plink(pathM[1], pathM[2], pathM[3])
samplenames=read.csv('ROSMAP_arrayGenotype.csv')
covariates_info=read.csv('DATA - neuroticism & demographics.csv')

#Sample filter
Sampstats <- row.summary(SNP_M$genotype)
callRate <- 0.95
Sampkeep=Sampstats$Call.rate>0.95& 
  (samplenames[[1]] %in% covariates_info$gwas_id&
     is.na())

#SNP filter
SNPstats <- col.summary(SNP_M$genotypes)
maf <- 0.01
geno.callrate=0.95
hardy=0.001
SNPkeep=with(SNPstats, Call.rate>geno.callrate& MAF>=maf&!is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ))
sum(SNPkeep)

SNP_mat_filtered=SNP_M$genotype[Sampkeep, SNPkeep]
rn_SNP_mat_filtered=samplenames[[1]][Sampkeep]
idx=match(rn_SNP_mat_filtered, covariates_info$gwas_id)

covariates_info_subset=covariates_info[idx,]

source('GWAA.R')
rownames(SNP_mat_filtered)=rn_SNP_mat_filtered
phenodata=covariates_info_subset[,c(2,10:16)]
phenodata$age_death=NULL
colnames(phenodata)[colnames(phenodata)%in% 'gwas_id']='id'
colnames(phenodata)[colnames(phenodata)%in% 'anxiety_20items']='phenotype'
phenodata$msex=as.factor(phenodata$msex)
phenodata$ros=as.factor(phenodata$ros)

SNP_mat_filtered=SNP_mat_filtered[!is.na(phenodata$phenotype),]
phenodata=phenodata[!is.na(phenodata$phenotype),]
phenodata$phenotype=rankNorm(phenodata$phenotype)
phenodata$id=as.character(phenodata$id)
library(RNOmni)
#hist(phenodata$phenotype)
#hist(rankNorm(phenodata$phenotype))
start <- Sys.time()
GWAA(genodata=SNP_mat_filtered,  phenodata=phenodata, family = gaussian, filename='GWAAresults',
                 append=FALSE, workers=4, flip=TRUE,
                 select.snps=NULL, hosts=NULL, nSplits=10)
elapsed=Sys.time() - start
library(installr)
os.sleep()
GWAAresults=read.table('GWAAresults', header=T)
GWAAresults$padjust=p.adjust(GWAAresults$p.value, 'fdr')
sort(GWAAresults$padjust)
