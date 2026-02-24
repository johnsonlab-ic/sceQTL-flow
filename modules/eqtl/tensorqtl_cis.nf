process run_tensorqtl {

    tag "${expression_mat} and ${fixed_pcs}"
    label "process_eqtl"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    path genotype_mat
    path snp_locations
    path expression_mat
    path gene_locations
    path fixed_pcs

    output:
    path "*_tensorqtl_done.txt", emit: run_marker
    path "*_tensorqtl.cis_qtl_pairs.*.parquet", emit: nominal_pairs
    path "*"

    script:
    """
    #!/usr/bin/env python3
    import time
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import torch
    from tensorqtl import cis


    def resolve_celltype(expression_file: str) -> str:
        name = Path(expression_file).name
        if name.endswith("_residuals.csv"):
            return name.replace("_residuals.csv", "")
        if name.endswith("_pseudobulk_normalised.csv"):
            return name.replace("_pseudobulk_normalised.csv", "")
        return Path(name).stem


    def log(msg: str) -> None:
        print(f"[tensorQTL] {msg}", flush=True)


    expression_path = Path("$expression_mat")
    genotype_path = Path("$genotype_mat")
    snp_locations_path = Path("$snp_locations")
    gene_locations_path = Path("$gene_locations")
    pcs_path = Path("$fixed_pcs")

    celltype = resolve_celltype(str(expression_path))
    prefix = f"{celltype}_tensorqtl"

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log(f"torch {torch.__version__}; cuda_available={torch.cuda.is_available()}; device={device}")

    t0 = time.perf_counter()

    log("Loading input tables")
    expression_df = pd.read_csv(expression_path)
    genotype_df = pd.read_csv(genotype_path)
    gene_locations_df = pd.read_csv(gene_locations_path)
    snp_locations_df = pd.read_csv(snp_locations_path)

    phenotype_df = (
        expression_df.set_index("geneid")
        .apply(pd.to_numeric, errors="coerce")
    )

    genotype_df = (
        genotype_df.set_index("snp")
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0)
        .astype("int8")
    )

    phenotype_pos_df = gene_locations_df.copy()
    if "geneid" not in phenotype_pos_df.columns:
        raise ValueError("gene_locations file must contain a 'geneid' column")
    if "chr" not in phenotype_pos_df.columns:
        raise ValueError("gene_locations file must contain a 'chr' column")
    if "pos" not in phenotype_pos_df.columns:
        if "left" in phenotype_pos_df.columns:
            phenotype_pos_df = phenotype_pos_df.rename(columns={"left": "pos"})
        else:
            raise ValueError("gene_locations file must contain 'pos' or 'left' column")
    phenotype_pos_df = phenotype_pos_df.set_index("geneid")[['chr', 'pos']]

    if "annot" not in snp_locations_df.columns:
        raise ValueError("snp locations file must contain an 'annot' column")
    if "position" in snp_locations_df.columns and "pos" not in snp_locations_df.columns:
        snp_locations_df = snp_locations_df.rename(columns={"position": "pos"})
    if "pos" not in snp_locations_df.columns:
        raise ValueError("snp locations file must contain 'pos' or 'position' column")
    if "chrom" not in snp_locations_df.columns:
        raise ValueError("snp locations file must contain a 'chrom' column")

    variant_df = snp_locations_df.rename(columns={"annot": "id"}).copy()
    if "index" not in variant_df.columns:
        variant_df["index"] = np.arange(variant_df.shape[0], dtype=np.int64)
    variant_df["pos"] = pd.to_numeric(variant_df["pos"], errors="coerce").astype("int32")
    variant_df = variant_df.set_index("id")[["chrom", "pos", "index"]]

    covariates_df = pd.read_csv(pcs_path, sep="\t", index_col=0).T
    covariates_df = covariates_df.apply(pd.to_numeric, errors="coerce")

    phenotype_df = phenotype_df[~phenotype_df.index.duplicated(keep="first")]
    phenotype_pos_df = phenotype_pos_df[~phenotype_pos_df.index.duplicated(keep="first")]

    common_genes = phenotype_df.index.intersection(phenotype_pos_df.index)
    phenotype_df = phenotype_df.reindex(common_genes)
    phenotype_pos_df = phenotype_pos_df.reindex(common_genes)

    common_samples = genotype_df.columns.intersection(phenotype_df.columns)
    common_samples = common_samples.intersection(covariates_df.index)
    genotype_df = genotype_df.loc[:, common_samples]
    phenotype_df = phenotype_df.loc[:, common_samples]
    covariates_df = covariates_df.loc[common_samples, :]

    common_snps = genotype_df.index.intersection(variant_df.index)
    genotype_df = genotype_df.loc[common_snps]
    variant_df = variant_df.loc[common_snps]

    log(f"Cell type: {celltype}")
    log(f"Input shapes -> genotype={genotype_df.shape}, phenotype={phenotype_df.shape}, variant={variant_df.shape}, phenotype_pos={phenotype_pos_df.shape}, covariates={covariates_df.shape}")

    if phenotype_df.shape[1] < 20:
        log(f"WARNING: only {phenotype_df.shape[1]} common samples; eQTL may be underpowered")

    log("Running cis.map_nominal")
    cis.map_nominal(
        genotype_df,
        variant_df,
        phenotype_df,
        phenotype_pos_df,
        prefix,
        covariates_df=covariates_df,
        output_dir='.',
        write_top=False,
        logger=None,
        verbose=True,
        window=int(${params.cis_distance})
    )

    elapsed = time.perf_counter() - t0
    log(f"Done in {elapsed/60:.2f} min")

    Path(f"{celltype}_tensorqtl_done.txt").write_text(
        f"celltype={celltype}\\n"
        f"samples={phenotype_df.shape[1]}\\n"
        f"phenotypes={phenotype_df.shape[0]}\\n"
        f"variants={variant_df.shape[0]}\\n"
        f"covariates={covariates_df.shape[1]}\\n"
        f"elapsed_seconds={elapsed:.2f}\\n"
    )
    """
}