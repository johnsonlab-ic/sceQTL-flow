process get_residuals {
    tag "${expression_mat}"
    label "process_high_memory"
    publishDir "${params.outdir}/residuals/", mode: 'copy'

    input:

    path expression_mat
    path cov_file
    path source_R

    output:
    path "*_residuals.csv", emit: residuals_results

    script:
    """
    #!/usr/bin/env Rscript
    source("$source_R")
    library(data.table)
    library(dplyr)

    exp_mat = fread("$expression_mat") %>% tibble::column_to_rownames(var="geneid")
    
    if(file.size(cov_file) > 0){
        cov_mat = fread(cov_file,head=T) %>% tibble::column_to_rownames(var="V1")

        #extract all colnames except "Individual_ID"
        covs_to_include = colnames(cov_mat)[!colnames(cov_mat) %in% c("Individual_ID")]

        exp_mat=get_residuals(exp_mat,cov_mat,covs_to_include=covs_to_include) %>% mutate(geneid=row.names(.))
    }else{
        exp_mat = exp_mat %>% mutate(geneid=row.names(.))
    }

    celltype = gsub("_pseudobulk_normalised.csv", "", "$expression_mat")
    data.table::fwrite(exp_mat, paste0(celltype, "_residuals.csv"))
    


    """

}