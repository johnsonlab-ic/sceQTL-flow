#!/usr/bin/env python3
"""Native per-file pseudobulk for one AnnData (.h5ad).

--celltype-col accepts a comma-separated list of columns to pseudobulk on
independently (e.g. different annotation resolutions). Individual/celltype
labels are resolved to a per-cell table [cell_id -> individual, celltype...]
from EITHER:
  - an external cell-level metadata file (--cell-metadata, keyed by --id-col), OR
  - the object's own obs (when --cell-metadata is omitted).
That table is aligned to the counts by cell_id, then for each requested
celltype column, counts are summed cells->individuals within each category
via a single sparse matmul.

Emits per (celltype column x category):
  <column>-<category>__<tag>_partial.csv   geneid + one raw-count column per individual
  <column>-<category>__<tag>_ncells.csv    individual, n_cells
The "<column>-" prefix keeps categories from different columns from
colliding (and being silently summed together) downstream.
"""
import argparse
import os
import re
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

SEP = "\x1f"  # unit separator: safe group-key delimiter


def sanitize(s):
    return re.sub(r"[^A-Za-z0-9._+-]+", "_", str(s))


def resolve_base(adata, celltype_cols, args):
    """Return (meta, mask): meta is a DataFrame indexed by cell_id with columns
    [indiv_col, *celltype_cols] for cells present in the label source (in obs
    order); mask is the boolean obs selector. indiv_col is cast to str;
    celltype columns are left as-is so NaN can still be detected per-column."""
    cell_ids = pd.Index(adata.obs_names.astype(str))
    needed = [args.indiv_col] + celltype_cols
    if args.cell_metadata:
        meta = pd.read_csv(args.cell_metadata, usecols=[args.id_col] + needed)
        meta[args.id_col] = meta[args.id_col].astype(str)
        meta = meta.drop_duplicates(args.id_col).set_index(args.id_col)
    else:
        for c in needed:
            if c not in adata.obs.columns:
                raise SystemExit(
                    f"ERROR: obs column '{c}' not found and no --cell-metadata given"
                )
        meta = adata.obs[needed].copy()
        meta.index = cell_ids

    mask = cell_ids.isin(meta.index)
    n_drop = int((~mask).sum())
    if n_drop:
        print(f"[PB] {n_drop}/{len(cell_ids)} cells not in label source; dropped")
    present = cell_ids[mask]
    if len(present) == 0:
        raise SystemExit("ERROR: no cells matched between counts and label source")
    meta = meta.loc[present, needed].copy()
    meta[args.indiv_col] = meta[args.indiv_col].astype(str)
    return meta, mask


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--celltype-col", required=True,
                    help="comma-separated list of celltype columns to pseudobulk on")
    ap.add_argument("--indiv-col", required=True)
    ap.add_argument("--cell-metadata", default=None,
                    help="optional cell-level metadata CSV/.gz; if omitted, labels read from obs")
    ap.add_argument("--id-col", default="cell_id",
                    help="cell-id column in --cell-metadata (matches obs_names)")
    ap.add_argument("--counts-layer", default="counts")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()

    tag = args.tag or os.path.splitext(os.path.basename(args.h5ad))[0]
    os.makedirs(args.outdir, exist_ok=True)
    celltype_cols = [c.strip() for c in args.celltype_col.split(",") if c.strip()]

    adata = ad.read_h5ad(args.h5ad)
    meta, mask = resolve_base(adata, celltype_cols, args)

    if args.counts_layer and args.counts_layer in adata.layers:
        X = adata.layers[args.counts_layer]
    else:
        print(f"[PB] layer '{args.counts_layer}' not found; using X")
        X = adata.X
    X = sp.csr_matrix(X)[mask]  # cells x genes, restricted to resolved cells (obs order)

    genes = np.asarray(adata.var_names)
    ind_all = meta[args.indiv_col].to_numpy()

    for col in celltype_cols:
        valid = meta[col].notna().to_numpy()
        n_invalid = int((~valid).sum())
        if n_invalid:
            print(f"[PB] [{col}] {n_invalid} cells with missing celltype label; dropped")

        ct = meta[col].to_numpy()[valid].astype(str)
        ind = ind_all[valid]
        Xc = X[valid]

        key = pd.Series(ct) + SEP + pd.Series(ind)
        codes, uniques = pd.factorize(key, sort=True)
        design = sp.csr_matrix(
            (np.ones(len(codes)), (np.arange(len(codes)), codes)),
            shape=(len(codes), len(uniques)),
        )
        summed = np.rint(np.asarray((Xc.T @ design).todense())).astype(np.int64)
        ncells = np.asarray(design.sum(axis=0)).ravel().astype(np.int64)

        grp_ct = np.array([u.split(SEP)[0] for u in uniques])
        grp_ind = np.array([u.split(SEP)[1] for u in uniques])

        for ctname in pd.unique(grp_ct):
            m = grp_ct == ctname
            inds = grp_ind[m]
            df = pd.DataFrame(summed[:, m], columns=inds)
            df.insert(0, "geneid", genes)
            label = sanitize(f"{col}-{ctname}")
            df.to_csv(os.path.join(args.outdir, f"{label}__{tag}_partial.csv"), index=False)
            pd.DataFrame({"individual": inds, "n_cells": ncells[m]}).to_csv(
                os.path.join(args.outdir, f"{label}__{tag}_ncells.csv"), index=False
            )
            print(f"[PB] [{col}] {ctname}: {summed.shape[0]} genes x {len(inds)} individuals (tag={tag})")


if __name__ == "__main__":
    main()
