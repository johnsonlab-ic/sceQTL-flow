

###main eqtl function 


calculate_ciseqtl=function(exp_mat,
  geno_mat,
  pcs,
  geno_pcs,
  exp_loc,
  geno_loc,
  name,
  cisDist=1e6,
  covadj=FALSE,
  covadj_pc=FALSE,
  geno_pcs_adj=TRUE,
  standardize=FALSE,
  cov_file,
  specify_genes=FALSE,
  specify_gene_list,
  covs_to_include,
  pvOutputThreshold=2e-5,
  filter_trans_FDR=FALSE,
  pvOutputThreshold_cis=5e-2,
  save_results=TRUE){

  exp_locs<-exp_loc
  geno_loc<-geno_loc

  common_names<-intersect(colnames(exp_mat),colnames(geno_mat))
  geno_mat<-geno_mat[common_names]
  exp_mat<-exp_mat[common_names]

  # this is in case the number of genes differ
  commongenes<-intersect(exp_locs$geneid,row.names(exp_mat))
  exp_mat<-exp_mat[commongenes,]

  # from MatrixEQTL manual: "The order of genes in gene locations
  # does not have to match the order of gene expression"


  cvrt=MatrixEQTL::SlicedData$new();
  if(covadj==TRUE){
    if(cov_file==""){
      covs<-as.data.frame(matrix(ncol=ncol(exp_mat),nrow=0))
      colnames(covs)<-exp_mat
    }else{
      covs<-read.table(cov_file)
      covs<-covs[covs_to_include,]
      covs<-covs[colnames(exp_mat)]
    }

    if(covadj_pc==TRUE){
      covs<-rbind(covs,pcs)
    }

    if(geno_pcs_adj==TRUE){
      covs=rbind(covs,geno_pcs)
    }


    message("Covs used: ",paste0(rownames(covs),sep=", "))
    cvrt$CreateFromMatrix(as.matrix(covs))
    saveRDS(covs,paste0(name,"covs_used.rds"))
  }

  snps=MatrixEQTL::SlicedData$new();
  snps$CreateFromMatrix(as.matrix(geno_mat))
  
  if(specify_genes==TRUE){
      exp_mat<-exp_mat[rownames(exp_mat) %in% specify_gene_list,]
  }

  if(standardize==TRUE){
    scaled<-scale(t(exp_mat),scale=T,center=F)
    exp_mat<-as.data.frame(t(scaled))
  }
  gene=MatrixEQTL::SlicedData$new();
  gene$CreateFromMatrix(as.matrix(exp_mat))
  n_indivs<-ncol(exp_mat)

  me<-suppressMessages(MatrixEQTL::Matrix_eQTL_main(snps,
                    gene,
                    cvrt = cvrt,
                    pvOutputThreshold = pvOutputThreshold,
                    useModel = MatrixEQTL::modelLINEAR,
                    errorCovariance = numeric(),
                    verbose = FALSE,
                    output_file_name=NULL,
                    output_file_name.cis =NULL,
                    pvOutputThreshold.cis =pvOutputThreshold_cis,
                    snpspos = geno_loc,
                    genepos = exp_locs,
                    cisDist = cisDist,
                    pvalue.hist = FALSE,
                    min.pv.by.genesnp = FALSE,
                    noFDRsaveMemory = FALSE))

  if(save_results==TRUE){


    if(pvOutputThreshold_cis>0){
      tmp<-me$cis$eqtls
      names(tmp)[which(names(tmp)=="statistic")]<-"t.stat"
      names(tmp)[which(names(tmp)=="pvalue")]<-"p.value"
      names(tmp)[which(names(tmp)=="snps")]<-"SNP"
      saveRDS(tmp,paste0(name,"_cis_MatrixEQTLout.rds"))
    }


    if(pvOutputThreshold>0){
      if(pvOutputThreshold_cis>0){
        tmp<-me$trans$eqtls
        names(tmp)[which(names(tmp)=="statistic")]<-"t.stat"
        names(tmp)[which(names(tmp)=="pvalue")]<-"p.value"
        names(tmp)[which(names(tmp)=="snps")]<-"SNP"
      if(filter_trans_FDR==TRUE){
          tmp<-tmp[tmp$FDR<0.2,]
        }
      saveRDS(tmp,paste0(name,"_trans_MatrixEQTLout_0.2FDR.rds"))
        }else{

          tmp<-me$all$eqtls
          if(filter_trans_FDR==TRUE){
              tmp<-tmp[tmp$FDR<0.2,]
            }

          saveRDS(tmp,paste0(name,"_trans_MatrixEQTLout_0.2FDR.rds"))
              }
      }
  }
    message(paste0("MatrixEQTL calculated for ",name, "."))
}



## PC optimisation eqtl
optimize_eqtl=function(exp_mat,
  geno_mat,
  pcs,
  geno_pcs,
  exp_loc,
  geno_loc,
  name,
  cisDist=1e6,
  filter_chr=TRUE,
  geno_pcs_adj=TRUE,
  chr="chr1",
  cov_file,
  covs_to_include,
  pvOutputThreshold=2e-5){

    
  exp_mat<-exp_mat
  geno_mat<-geno_mat
  exp_locs<-exp_loc
  geno_loc<-geno_loc

  common_names<-intersect(colnames(exp_mat),colnames(geno_mat))
  geno_mat<-geno_mat[common_names]
  exp_mat<-exp_mat[common_names]

  # this is in case the number of genes differ
  commongenes<-intersect(exp_locs$geneid,row.names(exp_mat))
  exp_mat<-exp_mat[commongenes,]


  ### get Max PCs - set to 100, less if exp_mat is small
  if(ncol(exp_mat)<100){
      
      
      down_signif <- function(x, digits = 0) {
          m <- 10^(ceiling(log(x, 10)) - digits)
          (x %/% m)*m
        }
      
      max_pcs=down_signif(ncol(exp_mat),1)/10
      max_pcs=max_pcs-1
  }else{
      max_pcs=10
      max_pcs=max_pcs-1
  }

  gene=MatrixEQTL::SlicedData$new();
  gene$CreateFromMatrix(as.matrix(exp_mat))
  n_indivs<-ncol(exp_mat)
  
  if(cov_file==""){
    covs<-as.data.frame(matrix(ncol=ncol(exp_mat),nrow=0))
    colnames(covs)<-exp_mat
  }else{
    covs<-read.table(cov_file)
    covs<-covs[covs_to_include,]
    covs<-covs[colnames(exp_mat)]
  }

  covs<-as.matrix(covs)
  ##add in pcs as provided by user
  covs<-rbind(covs,pcs)

  ##add Geno PCs
  if(geno_pcs_adj==TRUE){
      covs=rbind(covs,geno_pcs)
  }

  #filter to only contain chr1 to speed up optimisation
  if(filter_chr==TRUE){
      geno_loc=geno_loc[geno_loc$chrom %in% chr,]
      exp_loc=exp_locs[exp_locs$chr %in% chr,]
      geno_mat=geno_mat[rownames(geno_mat) %in% geno_loc$annot,]
    }
  
  snps=MatrixEQTL::SlicedData$new();
  snps$CreateFromMatrix(as.matrix(geno_mat))

  nsnps=nrow(geno_mat)
  ngenes=nrow(exp_mat)


  res_vector=vector()
  time_start=Sys.time()
  n_eqtls_results=vector()
  n_egenes_results=vector()
  pcs_vector=vector()
  
  message(paste0("A total of ",ncol(exp_mat)," individuals in exp_mat."))
  message(paste0("A total of ",ncol(geno_mat)," individuals in geno_mat."))
  message(paste0("A total of ",length(common_names)," common individuals."))



  for(i in 1:max_pcs){

        n_pcs=i*10
        pcs_vector=c(pcs_vector,n_pcs)

        covvec=c(covs_to_include,paste0(rep("PC."),1:(n_pcs)))
        tmp_covs=covs[rownames(covs) %in% covvec,]
        
        message(paste0(Sys.time()," Testing MatrixEQTL for ",n_pcs," PCs"))
        message(paste0(nsnps," SNPs and ",ngenes," genes used."))


        message("Covs used: ",paste0(rownames(tmp_covs),sep=", "))

        cvrt=MatrixEQTL::SlicedData$new();
        cvrt$CreateFromMatrix(as.matrix(tmp_covs))

  
        
        me<-suppressMessages(MatrixEQTL::Matrix_eQTL_main(snps,
                            gene,
                            cvrt=cvrt,
                            pvOutputThreshold = 0,
                            useModel = MatrixEQTL::modelLINEAR,
                            errorCovariance = numeric(),
                            verbose = FALSE,
                            output_file_name=NULL,
                            output_file_name.cis =NULL,
                            pvOutputThreshold.cis =5e-2,
                            snpspos = geno_loc,
                            genepos = exp_loc,
                            cisDist = cisDist,
                            pvalue.hist = FALSE,
                            min.pv.by.genesnp = FALSE,
                            noFDRsaveMemory = FALSE))

        eqtl_out=me$cis$eqtls
        n_eqtls=nrow(eqtl_out[eqtl_out$FDR<0.05,])
        n_egenes=length(unique(eqtl_out[eqtl_out$FDR<0.05,]$gene))
        message(paste0(n_egenes," eGenes discovered at < 5% FDR."))

        n_eqtls_results=c(n_eqtls_results,n_eqtls)
        n_egenes_results=c(n_egenes_results,n_egenes)
    }
  
  res_df=data.frame(n_eqtls_results,n_egenes_results,pcs=pcs_vector)
  res_df=res_df[order(res_df$n_egenes_results,decreasing=T),]
  res_df=res_df[1,]
  return(res_df)


}



## Function to filter the genes in the expression matrix. Takes either a minimum percentage of individuals or a minimum number of individuals
filter_pseudobulk = function(exp_mat, minimum_percentage=NULL, minimum_indivs=NULL){
  
  if(!is.null(minimum_percentage) && is.null(minimum_indivs)){

    message(paste0("filter_pseudobulk=TRUE. Keeping genes expressed in minimum ",minimum_percentage, " % individuals."))
    total_individuals = ncol(exp_mat)
    minimum_indivs = ceiling(total_individuals * (minimum_percentage / 100))
    exp_mat <- exp_mat[rowSums(exp_mat > 0) >= minimum_indivs, ]

  }else if(!is.null(minimum_indivs) && is.null(minimum_percentage)){
    exp_mat <- exp_mat[rowSums(exp_mat > 0) >= minimum_indivs, ]
     message(paste0("filter_pseudobulk=TRUE. Keeping genes expressed in minimum ",minimum_indivs, " number of individuals."))

  }else{
    stop("You must supply either a minimum_percentage or minimum_indivs value")
  }
  
  return(exp_mat)
}

## Filter genotype matrix. Retains snps with at least 2 individuals in each of the 2 genotypic categories.
filter_genotype_matrix=function(genotypemat,filter_type=c("2.2","2.3")){

  tmp<-genotypemat

  message(paste0("Filtering genotype matrix. Starting with ",nrow(genotypemat)," snps."))

  tmp$counter_0<-rowSums(tmp[1:ncol(tmp)]==0)
  tmp$counter_1<-rowSums(tmp[1:ncol(tmp)]==1)
  tmp$counter_2<-rowSums(tmp[1:ncol(tmp)]==2)

  if(filter_type=="2.2"){

    message("Retaining snps with at least 2 individuals in each of the 2 genotypic categories.")

    option1=tmp[tmp$counter_0>=2 & tmp$counter_1>=2,]
    option2=tmp[tmp$counter_0>=2 & tmp$counter_2>=2,]
    option3=tmp[tmp$counter_1>=2 & tmp$counter_2>=2,]

    intersected_rows<-unique(c(rownames(option1),rownames(option2),rownames(option3)))

  } else if (filter_type=="2.3"){

    message("Retaining snps with at least 2 individuals in each of the 3 genotypic categories.")
    tmp<-tmp[tmp$counter_0>=2 & tmp$counter_1>=2 & tmp$counter_2>=2,]
    intersected_rows<-row.names(tmp)
  }
  genotypemat<-genotypemat[intersected_rows,]
  message(paste0("Filtering complete. ",nrow(genotypemat), " snps retained."))
  return(genotypemat)
}


#function to get residuals (for now, specifically adds Random effects to Sample_Source and Diagnosis, and fixed effects to all other covariates)
get_residuals=function(exp_mat,covs_to_include,cov_file){

  covmat=read.table(cov_file)
  message(paste0("Correcting expression matrix for known covariates.\nUser-defined covariates: '",paste(covs_to_include,collapse=", "),"'"))  
  rownames(covmat)=covmat$Individual_ID
  ##
  if(length(grep(covs_to_include[1],rownames(covmat))>0)){
    covmat=as.data.frame(t(covmat))
  }

  ###relevel according to largest
  for (col_name in covs_to_include) {
    if (is.character(covmat[[col_name]])) {
      covmat[[col_name]] <- as.factor(covmat[[col_name]])
    }
    
    if (is.factor(covmat[[col_name]])) {
      # Determine the level with the highest frequency
      ref_level <- names(sort(table(covmat[[col_name]]), decreasing = TRUE))[1]
      
      # Relevel the factor with the most frequent level as the reference
      covmat[[col_name]] <- relevel(covmat[[col_name]], ref = ref_level)
    }
  }

  ### keep common names
  samples=colnames(exp_mat)
  common_samples=intersect(samples,rownames(covmat))

  covmat=covmat[common_samples,]
  exp_mat=exp_mat[,common_samples]



  lmmodel = paste0("gene ~ ", paste(covs_to_include, collapse = " + "))
  exp_mat=t(apply(exp_mat,1,function(x){
      
    lm_mat=data.frame(gene=x)
    lm_mat=cbind(lm_mat,covmat)
      # fit=lmerTest::lmer(lmmodel,lm_mat)
      fit=lm(lmmodel,lm_mat)
      resid(fit)
    
      
  }))
  return(exp_mat)

}


#function to get genotype pcs
get_geno_pcs=function(geno_mat){

  pcs=prcomp(geno_mat,scale=T)
  pcs=pcs$rotation
  pcs=pcs[,1:30]
  colnames(pcs)=paste0("geno.",colnames(pcs))
  return(t(pcs))



}

suppressMessages(library(argparse))




###########################################################
################ INPUT PARAMETERS #########################
###########################################################
parser <- ArgumentParser()

parser$add_argument("--script_dir", help = "path to scripts directory")
##### input files
###########################################################


parser$add_argument("--exp_mat", help = "Expression matrix path")
parser$add_argument("--geno_mat", help = "Genotype dosage matrix (0,1,2) path")
parser$add_argument("--exp_loc", help = "Gene locations path")
parser$add_argument("--geno_loc", help = "SNP locations path")
parser$add_argument("--cov_file", help = "Covariate file")
parser$add_argument("--covs_to_include",nargs="+",help = "")


##### processing params
###########################################################

parser$add_argument("--filter_genotype_matrix", 
help = "Filter genotype matrix",action="store_true")
parser$add_argument("--filter_pseudobulk_matrix", 
help = "Filter pseudobulk matrix",action="store_true")
parser$add_argument("--exp_mat_thresh_percent", 
help = "Filter pseudobulk matrix min indivs")
parser$add_argument("--get_residuals", 
help = "Gets residuals, correcting for covariates",action="store_true")
parser$add_argument("--specify_samples_file", 
help = "File with samples to include. Will be read with ReadLines().")
parser$add_argument("--trans_eqtls", 
help = "Whether or not to calculate trans-eQTLs 'true' or 'false' .")



args <- parser$parse_args(commandArgs(TRUE))



exp_mat=args$exp_mat
name=gsub(".*/([^/]+)_pseudobulk.csv", "\\1", exp_mat)
name=gsub("_pseudobulk.csv","",name)

geno_mat=args$geno_mat
exp_loc=args$exp_loc
geno_loc=args$geno_loc
cov_file=args$cov_file
if(length(cov_file)==0){
  cov_file=""
  cov_file_included=FALSE
}else{
  cov_file_included=TRUE
}

script_dir=args$script_dir
#source("${baseDir}/singlecell/qtl_pipeline/matrix_eqtl_funcs.r")
##source(paste0(script_dir,"/matrix_eqtl_funcs.r"))

covs=args$covs_to_include
print(covs)

##### processing params
###########################################################

trans_eqtls=args$trans_eqtls
if(trans_eqtls=="true"){
 trans_eqtls=TRUE
message("trans_eqtls set to 'true'. Will calculate trans-eQTls at nominal (p<0.05) significance.")
} else{
trans_eqtls=FALSE}

geno_filter_type="2.2"
filter_pseudobulk_thresh=as.numeric(args$exp_mat_thresh_percent)
if (!is.null(args$specify_samples_file)) {
  # Read the file and store the samples in specify_samples_vector
  specify_samples_vector <- readLines(args$specify_samples_file)
}
get_exp_residuals=args$get_residuals
ref_levels=NULL


###########################################################
###### Helper funcs
###########################################################
filter_pseudobulk = function(exp_mat, minimum_percentage=NULL, minimum_indivs=NULL){
  
  if(!is.null(minimum_percentage) && is.null(minimum_indivs)){

    message(paste0("filter_pseudobulk=TRUE. Keeping genes expressed in minimum ",minimum_percentage, " % individuals."))
    total_individuals = ncol(exp_mat)
    minimum_indivs = ceiling(total_individuals * (minimum_percentage / 100))
    exp_mat <- exp_mat[rowSums(exp_mat > 0) >= minimum_indivs, ]

  }else if(!is.null(minimum_indivs) && is.null(minimum_percentage)){
    exp_mat <- exp_mat[rowSums(exp_mat > 0) >= minimum_indivs, ]
     message(paste0("filter_pseudobulk=TRUE. Keeping genes expressed in minimum ",minimum_indivs, " number of individuals."))

  }else{
    stop("You must supply either a minimum_percentage or minimum_indivs value")
  }
  
  return(exp_mat)
}

get_residuals=function(exp_mat,covs_to_include,cov_file){

  covmat=read.table(cov_file)
  message(paste0("Correcting expression matrix for known covariates.\nUser-defined covariates: '",paste(covs_to_include,collapse=", "),"'"))  
  rownames(covmat)=covmat$Individual_ID
  ##
  if(length(grep(covs_to_include[1],rownames(covmat))>0)){
    covmat=as.data.frame(t(covmat))
  }

  ###relevel according to largest
  for (col_name in covs_to_include) {
    if (is.character(covmat[[col_name]])) {
      covmat[[col_name]] <- as.factor(covmat[[col_name]])
    }
    
    if (is.factor(covmat[[col_name]])) {
      # Determine the level with the highest frequency
      ref_level <- names(sort(table(covmat[[col_name]]), decreasing = TRUE))[1]
      
      # Relevel the factor with the most frequent level as the reference
      covmat[[col_name]] <- relevel(covmat[[col_name]], ref = ref_level)
    }
  }

  ### keep common names
  samples=colnames(exp_mat)
  common_samples=intersect(samples,rownames(covmat))

  covmat=covmat[common_samples,]
  exp_mat=exp_mat[,common_samples]



lmmodel = paste0("gene ~ ", paste(covs_to_include, collapse = " + "))

  exp_mat=t(apply(exp_mat,1,function(x){
      
    lm_mat=data.frame(gene=x)
    lm_mat=cbind(lm_mat,covmat)
      # fit=lmerTest::lmer(lmmodel,lm_mat)
      fit=lm(lmmodel,lm_mat)
      resid(fit)
    
      
  }))
  return(exp_mat)

}

##### matrixEQTL params
###########################################################


optimize_pcs=TRUE
filter_chr=FALSE
use_optimization_res=TRUE
cisDist=1e6

standardize=TRUE
pval_thresh_cis=1
if(trans_eqtls==TRUE){
 pval_thresh_trans=5e-7
pval_thresh_cis=0
} else{
 pval_thresh_trans=0
pval_thresh_cis=1
}

covadj_matrixeqtl=TRUE
covadj_pc_matrixeqtl=TRUE
geno_pc_adj=TRUE
geno_pcs_n=3



###########################################################
################ Read in files + reformat #################
###########################################################



exp_mat=read.table(exp_mat)
exp_loc=read.table(exp_loc)
geno_loc=as.data.frame(data.table::fread(geno_loc))
geno_mat=as.data.frame(data.table::fread(geno_mat))


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
message("Making sure geno_loc and geno_mat match..")
geno_mat<-geno_mat[rownames(geno_mat) %in% geno_loc$annot,]


message("Keeping common samples in all matrices..")
common_names=intersect(colnames(exp_mat),colnames(geno_mat))

if(cov_file_included==TRUE){
  covmat=read.table(cov_file)
  covmat=covmat[complete.cases(covmat),]
  # if(length(grep(covs[1],rownames(covmat))>0)){
  #   covmat=as.data.frame(t(covmat))
  # }
  common_names=intersect(common_names,covmat$Individual_ID)
}

exp_mat=exp_mat[,common_names]
geno_mat=geno_mat[,common_names]



###########################################################
###########################################################
#########        post-process / QC        #################
###########################################################
###########################################################



####################################
### specify samples if needed
####################################

if(!is.null(args$specify_samples_file)){

  message("Specify_samples=TRUE. Using a subset of the data.")
  specify_samples_vector=gsub("/",".",specify_samples_vector)
  specify_samples_vector=gsub("-",".",specify_samples_vector)
  
  common_names=intersect(colnames(exp_mat),specify_samples_vector)
  exp_mat=exp_mat[,common_names]
  geno_mat=geno_mat[,common_names]

}

message(paste0(length(common_names)," individuals retained. Double-check that this number is accurate."))

####################################
### filter expression matrix
####################################

if(args$filter_pseudobulk_matrix==TRUE){
    exp_mat=filter_pseudobulk(exp_mat,minimum_percentage=filter_pseudobulk_thresh)
}
####################################
### filter genotype matrix
####################################

if(args$filter_genotype_matrix==TRUE){
    geno_mat<-filter_genotype_matrix(geno_mat,filter_type=geno_filter_type)
    row.names(geno_loc)<-geno_loc$annot
    geno_loc<-geno_loc[rownames(geno_mat),]
    row.names(geno_loc)<-rep(1:nrow(geno_loc))
}

####################################
### STANDARDIZE EXPRESSION MATRIX
####################################
scaled<-scale(t(exp_mat),scale=T,center=F)
exp_mat<-as.data.frame(t(scaled))



####################################
### GET GENO PCS ###
####################################


geno_pcs=get_geno_pcs(geno_mat)
saveRDS(geno_pcs,paste0(name,"_geno_pcs.rds"))
geno_pcs=geno_pcs[1:geno_pcs_n,]



### GET RESIDUALS (2-STEP approach, if TRUE) ###
####################################

if(get_exp_residuals==TRUE){
 
 if(cov_file==""){
  message("You must supply a covariate file to obtain corrected expression residuals.")
 }else{
  message(paste0("Correcting expression matrix for known covariates.\nUser-defined covariates: '",paste(covs,collapse=", "),"'"))  
covs_to_include=covs
lmmodel = paste0("gene ~ ", paste(covs_to_include, collapse = " + "))

  exp_mat=suppressMessages(get_residuals(exp_mat,
  covs_to_include=covs,
  cov_file=cov_file))
  exp_mat=as.data.frame(exp_mat)
  write.table(exp_mat,paste0(name,"_residuals_pseudobulk.csv"))
  }
}

####################################
### GET PCS ###
####################################

if(ncol(exp_mat)<100){
    down_signif <- function(x, digits = 0) {
        m <- 10^(ceiling(log(x, 10)) - digits)
        (x %/% m)*m
      }
    
    max_pcs=down_signif(ncol(exp_mat),1)/10
}else{
    max_pcs=10
}
message("Getting PCs...")
pcs<-prcomp(exp_mat,scale=F,center=F)
pcs<-pcs$rotation
pcs<-pcs[,1:(max_pcs*10)]
pcs<-t(pcs)
pcs<-pcs[,colnames(exp_mat)]
rownames(pcs)<-paste0(rep("PC."),1:(max_pcs*10))
if(get_exp_residuals==TRUE){
  write.table(pcs,paste0(name,"_full_residualexpression_PCs.txt"))
}else{
  write.table(pcs,paste0(name,"_full_expression_PCs.txt"))
}


####################################
### OPTIMIZE N_PCS #################
####################################

if(optimize_pcs==TRUE && covadj_pc_matrixeqtl==TRUE){

    message("Optimize PCs = TRUE. Running Optimisation..")

    if(get_exp_residuals==TRUE){
        message("Expression matrix already corrected for covariates. Setting cov_file as null.")
        cov_file=""
        # geno_pc_adj=FALSE
      }

    best_pcs=optimize_eqtl(exp_mat=exp_mat,
        exp_loc=exp_loc,
        geno_mat=geno_mat,
        pcs=pcs,
        geno_pcs=geno_pcs,
        geno_loc=geno_loc,
        filter_chr=filter_chr,
        cisDist=cisDist,
        name=name,
        cov_file=cov_file,
        geno_pcs_adj=geno_pc_adj,
        covs_to_include=covs)

    best_pcs=best_pcs$pcs[1]


}

message(paste0("Optimization complete. ",best_pcs," PCs will be used."))

###########################################################
#########        MatrixEQTL                #################
###########################################################

nsnps=nrow(geno_mat)
ngenes=nrow(exp_mat)
nsamples=length(intersect(colnames(exp_mat),colnames(geno_mat)))
message(paste0(nsamples, " individuals will be used."))
message(paste0(nsnps," SNPs and ",ngenes," genes will be included."))

if(optimize_pcs==TRUE){
  pcs=pcs[1:best_pcs,]
}


calculate_ciseqtl(exp_mat=exp_mat,
exp_loc=exp_loc,
geno_mat=geno_mat,
geno_loc=geno_loc,
cisDist=1e6,
name=name,
specify_genes=FALSE,
specify_gene_list=specify_gene_list,
pvOutputThreshold=pval_thresh_trans,
pvOutputThreshold_cis=1,
covadj=covadj_matrixeqtl,
covadj_pc=covadj_pc_matrixeqtl,
pcs=pcs,
geno_pcs=geno_pcs,
geno_pcs_adj=geno_pc_adj,
standardize=TRUE,
cov_file=cov_file)
