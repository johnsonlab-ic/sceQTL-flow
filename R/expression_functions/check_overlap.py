#!/usr/bin/env python3
"""Genotype <-> single-cell sample overlap gate.

Lists genotype individuals and single-cell individuals, reports overlap.
  - WARN if shared < --warn-frac of single-cell individuals.
  - FAIL (exit 1) if there is zero overlap.

Genotype IDs are the column names of a genotype matrix CSV. If the matrix uses
IDs that differ from the single-cell individual labels (e.g. chip_id vs caseid),
pass --id-map / --map-from / --map-to to translate before comparing.
"""
import argparse
import sys
import pandas as pd


def read_genotype_ids(path):
    cols = list(pd.read_csv(path, nrows=0).columns)
    drop = {"", "snp", "Unnamed: 0"}
    return [c for c in cols if c not in drop]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genotype-matrix", required=True)
    ap.add_argument("--cell-metadata", required=True)
    ap.add_argument("--indiv-col", required=True, help="individual column in --cell-metadata")
    ap.add_argument("--id-map", default=None, help="optional CSV mapping genotype IDs -> individual")
    ap.add_argument("--map-from", default=None, help="column in --id-map matching genotype IDs")
    ap.add_argument("--map-to", default=None, help="column in --id-map giving the individual label")
    ap.add_argument("--warn-frac", type=float, default=0.5)
    args = ap.parse_args()

    geno = read_genotype_ids(args.genotype_matrix)
    if args.id_map:
        m = pd.read_csv(args.id_map)
        lut = dict(zip(m[args.map_from].astype(str), m[args.map_to].astype(str)))
        mapped = [(g, lut.get(str(g))) for g in geno]
        unmapped = [g for g, v in mapped if v is None]
        geno = [v for _, v in mapped if v is not None]
        if unmapped:
            print(f"[OVERLAP] {len(unmapped)} genotype IDs not in id-map (e.g. {unmapped[:3]})")

    geno_set = set(map(str, geno))
    sc = pd.read_csv(args.cell_metadata, usecols=[args.indiv_col])[args.indiv_col].astype(str)
    sc_set = set(sc.unique())

    inter = geno_set & sc_set
    print(f"[OVERLAP] genotype individuals: {len(geno_set)}  e.g. {sorted(geno_set)[:3]}")
    print(f"[OVERLAP] single-cell individuals: {len(sc_set)}  e.g. {sorted(sc_set)[:3]}")
    pct = 100 * len(inter) / max(len(sc_set), 1)
    print(f"[OVERLAP] shared: {len(inter)}  ({pct:.1f}% of single-cell individuals)")
    only_sc = sorted(sc_set - geno_set)
    only_geno = sorted(geno_set - sc_set)
    if only_sc:
        print(f"[OVERLAP] single-cell without genotype: {len(only_sc)}  e.g. {only_sc[:5]}")
    if only_geno:
        print(f"[OVERLAP] genotype without single-cell: {len(only_geno)}  e.g. {only_geno[:5]}")

    if len(inter) == 0:
        print("[OVERLAP] FAIL: no overlap between genotype and single-cell individuals.")
        sys.exit(1)
    if len(inter) < args.warn_frac * len(sc_set):
        print(f"[OVERLAP] WARNING: overlap below {100*args.warn_frac:.0f}% of single-cell individuals.")
    print("[OVERLAP] PASS")


if __name__ == "__main__":
    main()
