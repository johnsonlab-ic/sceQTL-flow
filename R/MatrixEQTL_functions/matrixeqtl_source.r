


calculate_ciseqtl=function(exp_mat,
  geno_mat,
  exp_loc,
  geno_loc,
  name="celltype",
  cisDist=1e6,
  covmat=NULL,
  pvOutputThreshold=2e-5,
  filter_trans_FDR=FALSE,
  pvOutputThreshold_cis=5e-2,
  save_results=TRUE){

  if(!is.null(covmat)){
    message("Including covariates (covmat!= 'NULL').")
    covs=rownames(covmat)
    message("Covs used: ",paste0(covs,sep=", "))

    covs_meqtl=MatrixEQTL::SlicedData$new();
    covs_meqtl=covs$CreateFromMatrix(as.matrix(covs))
    saveRDS(covs,paste0(name,"_covs_used.rds"))
    
  }else{
    covs_meqtl=NULL
  }

  ##create SNP input
  geno_meqtl=MatrixEQTL::SlicedData$new();
  geno_meqtl$CreateFromMatrix(as.matrix(geno_mat))

  expr_meqtl=MatrixEQTL::SlicedData$new();
  expr_meqtl$CreateFromMatrix(as.matrix(exp_mat))
  n_indivs<-ncol(exp_mat)

  me<-suppressMessages(MatrixEQTL::Matrix_eQTL_main(geno_meqtl,
                    expr_meqtl,
                    cvrt = covs_meqtl,
                    useModel = MatrixEQTL::modelLINEAR,
                    errorCovariance = numeric(),
                    verbose = FALSE,
                    output_file_name=NULL,
                    output_file_name.cis =NULL,
                    pvOutputThreshold.cis =pvOutputThreshold_cis,
                    pvOutputThreshold = pvOutputThreshold,
                    snpspos = geno_loc,
                    genepos = exp_loc,
                    cisDist = cisDist,
                    pvalue.hist = FALSE,
                    min.pv.by.genesnp = FALSE,
                    noFDRsaveMemory = FALSE))

                        if (save_results) {
                            save_eqtls <- function(eqtls, prefix) {
                                names(eqtls)[names(eqtls) == "statistic"] <- "t.stat"
                                names(eqtls)[names(eqtls) == "pvalue"] <- "p.value"
                                names(eqtls)[names(eqtls) == "snps"] <- "SNP"
                                saveRDS(eqtls, paste0(name, "_", prefix, "_MatrixEQTLout.rds"))
                            }

                            if (pvOutputThreshold_cis > 0) {
                                save_eqtls(me$cis$eqtls, "cis")
                            }

                            if (pvOutputThreshold > 0) {
                                eqtls <- if (pvOutputThreshold_cis > 0) me$trans$eqtls else me$all$eqtls
                                if (filter_trans_FDR) {
                                    eqtls <- eqtls[eqtls$FDR < 0.2, ]
                                }
                                save_eqtls(eqtls, "trans_0.2FDR")
                            }
                        }

                        message(paste0("MatrixEQTL calculated for ", name, "."))
  }