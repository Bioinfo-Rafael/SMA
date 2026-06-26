#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import scanpy as sc
import pandas as pd

warnings.simplefilter("ignore", category=FutureWarning)

# =========================================================
# input / output
# =========================================================
INPUT_H5AD = Path(
    #"/home/suzuki/Learn/SMA/v2/results/08_classical_full_inner_microglia_reclustering/02_full_clustering/full_inner_classical_clustered.h5ad"
    "/home/suzuki/Learn/SMA/v2/results/08_classical_full_inner_microglia_reclustering/04_microglia_reclustering/microglia_classical_reclustered.h5ad"
)

OUTDIR = Path(
    #"/home/suzuki/Learn/SMA/v2/results/08_classical_full_inner_microglia_reclustering/02_full_clustering/additional_marker_panels"
    "/home/suzuki/Learn/SMA/v2/results/08_classical_full_inner_microglia_reclustering/04_microglia_reclustering/marker_genes"
)

CLUSTER_COL ="microglia_leiden_r1_5" #"leiden_r1_5"
LAYER = "logexpr_for_clustering"

# =========================================================
# marker groups
# =========================================================
marker_groups = {
    "Microglia_resident_microglia": [
        "P2ry12", "Tmem119", "Cx3cr1", "Sall1", "Hexb", "Fcrls", "Olfml3",
        "Gpr34", "P2ry13", "Siglech", "Slc2a5", "Csf1r", "Aif1", "C1qa",
        "C1qb", "C1qc", "Ctss"
    ],
    "DAM_activated_microglia": [
        "Apoe", "Tyrobp", "Trem2", "Gpnmb", "Lpl", "Cst7", "Cd68", "Itgax",
        "Axl", "Clec7a", "Cd9", "Cd63", "Spp1", "Lgals3", "Ctsb", "Ctsd", "Ctsz"
    ],
    "Neuron": [
        "Snap25", "Syt1", "Rbfox1", "Rbfox2", "Rbfox3", "Tubb3", "Map2",
        "Nefl", "Nefm", "Stmn2", "Slc17a7", "Slc17a6", "Gad1", "Gad2", "Slc32a1"
    ],
    "Astrocyte": [
        "Aqp4", "Aldh1l1", "Slc1a2", "Slc1a3", "Gja1", "Sox9", "S100b",
        "Gfap", "Agt", "Sparcl1", "Clu"
    ],
    "Reactive_astrocyte": [
        "Gfap", "Vim", "Serpina3n", "Lcn2", "C3", "Hif3a"
    ],
    "Oligodendrocyte": [
        "Plp1", "Mbp", "Mog", "Mag", "Mobp", "Cnp", "Mal", "Opalin",
        "Car2", "Ugt8a", "Ermn"
    ],
    "OPC_oligodendrocyte_precursor": [
        "Pdgfra", "Cspg4", "Vcan", "Olig1", "Olig2", "Sox10", "Tnr",
        "Bcas1", "Enpp6", "Tcf7l2"
    ],
    "Endothelial": [
        "Pecam1", "Cldn5", "Kdr", "Flt1", "Tek", "Ly6c1", "Slco1a4",
        "Bsg", "Esam", "Vwf"
    ],
    "Pericyte_vascular_mural": [
        "Pdgfrb", "Rgs5", "Kcnj8", "Abcc9", "Notch3", "Vtn", "Des",
        "Acta2", "Tagln", "Myh11"
    ],
    "Ependymal": [
        "Foxj1", "Tmem212", "Dnah5", "Dnah12", "Cfap126", "Nnat",
        "Pifo", "Rsph1"
    ],
    "Schwann_cell_peripheral_nerve_contamination": [
        "Mpz", "Pmp22", "Mbp", "Plp1", "Sox10", "Ncmap", "Prx"
    ],
    "Meningeal_fibroblast_leptomeningeal": [
        "Dcn", "Col1a1", "Col1a2", "Lum", "Pdgfra", "Fn1", "Cxcl12", "Pi16"
    ],
    "Monocyte_macrophage_contamination": [
        "Ccr2", "Ly6c2", "Lyz2", "S100a8", "S100a9", "Ms4a7", "Fcgr1",
        "Itgam", "Cd14", "Mrc1", "Cd163", "Pf4", "Lyve1", "Folr2"
    ],
    "Tcell_NK": [
        "Cd3d", "Cd3e", "Trac", "Lck", "Nkg7", "Gzma", "Gzmb", "Klrb1c"
    ],
    "Bcell_plasma": [
        "Cd79a", "Cd79b", "Ms4a1", "Cd19", "Jchain", "Mzb1", "Xbp1"
    ],
    "Neutrophil": [
        "S100a8", "S100a9", "Mpo", "Elane", "Lcn2", "Retnlg", "Cxcr2"
    ],
}

# =========================================================
# utility
# =========================================================
def log(msg):
    print(msg, flush=True)

def save_current_figure(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close("all")

def build_upper_map(var_names):
    d = {}
    for g in var_names:
        d.setdefault(str(g).upper(), str(g))
    return d

def resolve_present_genes(adata, groups):
    upper = build_upper_map(adata.var_names)
    resolved = {}
    records = []

    for group, genes in groups.items():
        present = []
        for g in genes:
            matched = upper.get(str(g).upper())
            records.append({
                "group": group,
                "requested_gene": g,
                "present": matched is not None,
                "matched_var_name": matched if matched is not None else ""
            })
            if matched is not None:
                present.append(matched)

        # 同一group内重複除去
        seen = set()
        present_unique = []
        for x in present:
            if x not in seen:
                seen.add(x)
                present_unique.append(x)

        if len(present_unique) > 0:
            resolved[group] = present_unique

    return resolved, pd.DataFrame(records)

# =========================================================
# main
# =========================================================
def main():
    if not INPUT_H5AD.exists():
        raise FileNotFoundError(f"Input h5ad not found: {INPUT_H5AD}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "combined").mkdir(parents=True, exist_ok=True)
    (OUTDIR / "by_group").mkdir(parents=True, exist_ok=True)

    log(f"Reading: {INPUT_H5AD}")
    adata = sc.read_h5ad(INPUT_H5AD)

    if CLUSTER_COL not in adata.obs.columns:
        raise ValueError(
            f"{CLUSTER_COL} not found in adata.obs.\n"
            f"Available obs columns:\n{list(adata.obs.columns)}"
        )

    # plotting用に logexpr layer を .X に入れた copy を使う
    adata_plot = adata.copy()
    if LAYER in adata.layers:
        log(f"Using layer '{LAYER}' for plotting.")
        adata_plot.X = adata.layers[LAYER].copy()
    else:
        log(f"Layer '{LAYER}' not found. Using adata.X as-is for plotting.")

    resolved_groups, presence_df = resolve_present_genes(adata_plot, marker_groups)

    presence_path = OUTDIR / "marker_presence_additional_panels.csv"
    presence_df.to_csv(presence_path, index=False)
    log(f"Saved marker presence table: {presence_path}")

    if len(resolved_groups) == 0:
        raise ValueError("None of the requested marker genes were found in adata.var_names.")

    # group summary
    summary_rows = []
    for group, genes in resolved_groups.items():
        summary_rows.append({
            "group": group,
            "n_present_genes": len(genes),
            "present_genes": ", ".join(genes)
        })
    pd.DataFrame(summary_rows).to_csv(
        OUTDIR / "marker_group_summary.csv", index=False
    )

    # -------------------------
    # combined dotplot
    # -------------------------
    log("Creating combined dotplot...")
    sc.pl.dotplot(
        adata_plot,
        var_names=resolved_groups,
        groupby=CLUSTER_COL,
        standard_scale="var",
        show=False
    )
    save_current_figure(OUTDIR / "combined" / f"dotplot_all_markers_by_{CLUSTER_COL}.png")

    # -------------------------
    # combined tracksplot
    # -------------------------
    log("Creating combined tracksplot...")
    sc.pl.tracksplot(
        adata_plot,
        var_names=resolved_groups,
        groupby=CLUSTER_COL,
        show=False
    )
    save_current_figure(OUTDIR / "combined" / f"tracksplot_all_markers_by_{CLUSTER_COL}.png")

    # -------------------------
    # per-group dotplot / tracksplot
    # -------------------------
    for group, genes in resolved_groups.items():
        safe_group = group.replace("/", "_").replace(" ", "_")

        log(f"Creating plots for group: {group}")

        group_dir = OUTDIR / "by_group" / safe_group
        group_dir.mkdir(parents=True, exist_ok=True)

        # dotplot
        sc.pl.dotplot(
            adata_plot,
            var_names=genes,
            groupby=CLUSTER_COL,
            standard_scale="var",
            show=False
        )
        save_current_figure(group_dir / f"dotplot_{safe_group}_by_{CLUSTER_COL}.png")

        # tracksplot
        sc.pl.tracksplot(
            adata_plot,
            var_names=genes,
            groupby=CLUSTER_COL,
            show=False
        )
        save_current_figure(group_dir / f"tracksplot_{safe_group}_by_{CLUSTER_COL}.png")

    log("Done.")

if __name__ == "__main__":
    main()