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