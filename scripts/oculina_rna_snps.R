######################################### 
# Analyses of RNA-seq derived snp for Oculina arbuscula using R
# This script is associated with the manuscript 
# Rivera, H.E. and Davies, S.W. 
#
# Data are available from https://github.com/hrivera28/Oculina_sym-apo_gene_expression


######################################### 
# Load libraries
library(plyr)
library(dplyr)
library(adegenet)
library(vcfR)
library(hierfstat)
library(ape)
library(StAMPP)
library(pegas)
library(ggplot2)
library(poppr)
library(ggpubr)

######################################### 
# Import data 
# This snp matrix contains bi-allelic SNPs produced from running the Stacks pipeline using a P.lutea 
# draft genome, see SI for full details on filtering steps that produced the final dataset used here 

# This contains 1,289 bi-allelic SNPs from 75 individuals 
snp_vcf<-read.vcfR("Oculina_snps.recode.vcf")
# Check metadata
head(snp_vcf)

# Convert vcf to genind and genlight object for using adegenet package
snp_gen<-vcfR2genlight(snp_vcf, n.cores=5)
snp_genid<-vcfR2genind(snp_vcf)

genind2df(snp_genid)

# Import population strata dataframe
Strata<-c("A", "B", "C", "D", "E", "F")
Strata<-as.data.frame(Strata)
colnames(Strata)<-"Colony"
strata(snp_genid)<-Strata

#Set populations to Site first 
setPop(snp_genid)<-~Colony
popNames(snp_genid)

#Check metadata
head(snp_gen)
head(snp_genid)

snp_clone<-as.genclone(snp_genid)
info_table(snp_clone, plot=TRUE)
#make pop level all the same to ID clones 
pop(snp_clone)<-c(rep("A", 6))
Strata$Colony<-c(rep("A", 6))
strata(snp_clone)<-Strata

clonecorrect(snp_clone, strata = 1, combine = FALSE, keep = 1)
mlg(snp_clone)
# This doesn't really find any clones because it wants perfect matches of MLGs 

# Calculate individual genetic similarity instead 
snp_df<-genind2df(snp_genid)
snp_dist<-dist.gene(snp_df, method="pairwise", pairwise.deletion = FALSE, variance=FALSE)
hc=hclust(snp_dist,"ave")

#Figure S7
plot(hc,cex=0.7, main="Dendogram of Pairwise \n Genetic Distances", xlab="", ylab="" )

#########################################
# Multivariate Analyses 
Snpca_data<-scaleGen(snp_genid, NA.method="mean")
Snpca<-dudi.pca(Snpca_data, cent=FALSE, scale=FALSE,scannf =TRUE, nf=5)
col<-funky(10)
s.class(Snpca$li, pop(snp_genid), col=transp(col, 0.8), axesell = FALSE, cstar = 0, clabel=0.75, cpoint=3)
s.class(Snpca$li, pop(snp_genid), col=transp(col, 1), axesell = FALSE, cstar = 0, xax=1, yax=2, clabel=0.5, cpoint=3)
s.class(Snpca$li, pop(snp_genid), col=transp(col, 1), axesell = FALSE, cstar = 0, xax=2, yax=3, clabel=0.5, cpoint=3)

barplot(Snpca$eig[1:5],main="PCA eigenvalues", col=heat.colors(5))

#### Use dataframe from pca output to plot pca with ggplot seperately
PCA_points<-Snpca$li
PCA_points$Colony=c("A", "B", "C","D", "E", "F")
rownames(PCA_points)<-seq(1,6)

A<-ggplot(PCA_points, aes(x=Axis1, y=Axis2, colour=Colony))+
   geom_point(size=5)+theme_minimal()+xlab("PC1 (29.6%)")+ylab("PC2 (22.45%)")

B<-ggplot(PCA_df, aes(x=PC1, y=PC3, colour=Colony))+
  geom_point(size=5)+theme_minimal()+xlab("PC1 (29.6%)")+ylab("PC3 (22.14%)")


C<-ggplot(PCA_df, aes(x=PC2, y=PC3, colour=Colony))+
  geom_point(size=5)+theme_minimal()+xlab("PC2 (22.45%)")+ylab("PC3 (22.14%)")


D<-ggplot(PCA_df, aes(x=PC3, y=PC4, colour=Colony))+
  geom_point(size=5)+theme_minimal()+ylab("PC4 (20.58%)")+xlab("PC3 (22.14%)")

#Figure S8
ggarrange(A,B,C,D, ncol=2, nrow = 2,labels = c("A", "B", "C", "D"), common.legend = TRUE, legend = "bottom")

