process optimize_pcs {
    tag "${expression_mat} and ${n_pcs} PCs"

    label "process_eqtl"
    publishDir "${params.outdir}/optimization/", mode: 'copy'

    input:
    path source_R
    path genotype_mat
    path snp_locations
    path expression_mat
    path gene_locations
    val n_pcs

    output:
    path "*egenes*.txt" , emit: egenes_results

    script:
    """
    #!/usr/bin/env Rscript
    source("$source_R")
    library(data.table)
    library(dplyr)

    exp_mat = fread("$expression_mat") %>% tibble::column_to_rownames(var="geneid")
    geno_mat = fread("$genotype_mat") %>% tibble::column_to_rownames(var="snp")
    geno_loc = fread("$snp_locations")
    exp_loc = fread("$gene_locations")

    # Update to handle both residuals and normalized file naming patterns
    if (grepl("_residuals.csv\$", basename("$expression_mat"))) {
        celltype = gsub("_residuals.csv", "", basename("$expression_mat"))
    } else if (grepl("_pseudobulk_normalised.csv\$", basename("$expression_mat"))) {
        celltype = gsub("_pseudobulk_normalised.csv", "", basename("$expression_mat"))
    } else {
        # Fallback: remove .csv extension
        celltype = gsub(".csv\$", "", basename("$expression_mat"))
    }
    common_samples = intersect(colnames(exp_mat), colnames(geno_mat))

    exp_mat = exp_mat %>% select(all_of(common_samples))
    geno_mat = geno_mat %>% select(all_of(common_samples))

    common_genes = intersect(exp_loc %>% pull(geneid), rownames(exp_mat))
    exp_mat = exp_mat %>% filter(rownames(exp_mat) %in% common_genes)

    geno_loc = geno_loc[, c("annot", "chrom", "position")] %>% tibble::column_to_rownames(var="annot")
    geno_mat = geno_mat[rownames(geno_loc), ]
    geno_mat = geno_mat[complete.cases(geno_mat), ]
    geno_loc = geno_loc[rownames(geno_mat), ]
    geno_loc = geno_loc %>% mutate(annot = rownames(geno_loc)) %>% select(annot, chrom, position)

    covmat = NULL

    exp_pcs = prcomp(t(exp_mat), scale. = TRUE)
    exp_pcs = exp_pcs\$x[, 1:${n_pcs}]
    exp_pcs = as.data.frame(exp_pcs)
    colnames(exp_pcs) = paste0("PC", 1:${n_pcs})
    exp_pcs = t(exp_pcs)

    if (!is.null(covmat)) {
        covmat = rbind(covmat, exp_pcs)
    } else {
        covmat = exp_pcs
    }

    ##finally, re-order inputs to same column order
    exp_mat = exp_mat[, common_samples]
    geno_mat = geno_mat[, common_samples]
    covmat = covmat[, common_samples]
    
    outs=calculate_ciseqtl(
        exp_mat = exp_mat,
        covmat=covmat,
        exp_loc = exp_loc,
        geno_mat = geno_mat,
        geno_loc = geno_loc,
        name = celltype,
        cisDist = ${params.cis_distance}
    )
    
    n_egenes = outs %>% filter(FDR<0.05) %>% pull(gene) %>% unique() %>% length()
    write.table(data.frame(celltype=celltype,n_pcs=${n_pcs}, n_egenes=n_egenes), file=paste0(celltype,"_egenes_vs_",${n_pcs},".txt"), sep="\t", quote=FALSE, row.names=FALSE)
    """
}