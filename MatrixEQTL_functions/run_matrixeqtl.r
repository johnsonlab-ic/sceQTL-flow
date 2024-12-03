





suppressMessages(library(argparse))
parser <- ArgumentParser()

parser$add_argument("--exp_mat", help = "Expression matrix path")
parser$add_argument("--geno_mat", help = "Genotype dosage matrix (0,1,2) path")
parser$add_argument("--exp_loc", help = "Gene locations path")
parser$add_argument("--geno_loc", help = "SNP locations path")
parser$add_argument("--cov_file", help = "Covariate file")



##read in files
exp_mat=read.table(exp_mat)
exp_loc=read.table(exp_loc)
geno_loc=as.data.frame(data.table::fread(geno_loc))
geno_mat=as.data.frame(data.table::fread(geno_mat))


##reformat
rownames(geno_mat)<-geno_mat$snp
geno_mat$snp<-NULL
colnames(geno_mat)<-gsub("/",".",colnames(geno_mat))
colnames(geno_mat)<-gsub("-",".",colnames(geno_mat))
geno_loc<-geno_loc[,c("annot","chrom","position")]
row.names(geno_loc)<-geno_loc$annot
geno_mat<-geno_mat[rownames(geno_loc),]
geno_mat<-geno_mat[complete.cases(geno_mat),]
geno_loc<-geno_loc[rownames(geno_mat),]
row.names(geno_loc)<-rep(1:nrow(geno_loc))
geno_mat<-geno_mat[rownames(geno_mat) %in% geno_loc$annot,]



message("Keeping common samples in all matrices..")
common_names=intersect(colnames(exp_mat),colnames(geno_mat))

if(cov_file_included==TRUE){
  covmat=read.table(cov_file)
  covmat=covmat[complete.cases(covmat),]
  common_names=intersect(common_names,covmat$Individual_ID)
}

exp_mat=exp_mat[,common_names]
geno_mat=geno_mat[,common_names]
