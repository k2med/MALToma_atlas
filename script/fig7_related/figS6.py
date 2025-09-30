# ---- packages ----
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ---- global rcParams (editor-friendly text in vector outputs) ----
rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42
rcParams["svg.fonttype"] = "none"
rcParams["text.usetex"] = False
rcParams["font.sans-serif"] = ["Arial"]
rcParams["font.size"] = 15
rcParams["axes.titlesize"] = 15
rcParams["axes.labelsize"] = 15
rcParams["xtick.labelsize"] = 15
rcParams["ytick.labelsize"] = 15

# ---- paths & inputs ----
sp_data_folder = "../data/spatial_rnaseq/"
sample_names = ["S31_ST", "S32_ST", "S71_ST", "S01_ST", "S02_ST", "S89_ST", "S04_ST", "S06_ST", "S87_ST"]

# ---- plotting params ----
marker_genes = ["MS4A1", "CD3D", "CD4", "CD8A", "CXCL13", "PECAM1", "VEGFA", "COL1A1"]
metric_name = "n_genes_by_counts"  # or "total_counts"
size_val = 1.3                     # spot size
img_key = "hires"
alpha_img = 1.0

# ---- helpers ----
def get_library_id(adata):
    """Return the (only) Visium library_id stored in adata.uns['spatial']."""
    if "spatial" not in adata.uns or not isinstance(adata.uns["spatial"], dict) or len(adata.uns["spatial"]) == 0:
        raise ValueError("adata.uns['spatial'] is missing or empty.")
    # after normalization we store {sample_name: original_spatial_dict}
    return list(adata.uns["spatial"].keys())[0]


def read_and_qc(sample_name, path=sp_data_folder):
    """
    Read one 10x Visium sample and perform basic QC annotations.

    Parameters
    ----------
    sample_name : str
        Sample ID (folder under `path`).
    path : str
        Base folder containing the sample.

    Returns
    -------
    adata : anndata.AnnData
    """
    adata = sc.read_visium(
        path + str(sample_name),
        count_file="filtered_feature_bc_matrix.h5",
        load_images=True,
    )
    adata.var_names_make_unique()

    # annotate sample and gene symbols
    adata.obs["sample"] = str(sample_name)
    adata.var["SYMBOL"] = adata.var_names

    # mitochondrial flags and QC metrics (works with sparse matrices)
    adata.var["MT"] = adata.var["SYMBOL"].str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["MT"], log1p=False, inplace=True)
    # MT fraction is available as 'pct_counts_MT' (percent); keep also a ratio if needed
    adata.obs["MT_frac"] = adata.obs["total_counts_MT"] / adata.obs["total_counts"]

    # prefix obs_names by sample to avoid collisions across slides
    adata.obs_names = adata.obs["sample"].astype(str) + "_" + adata.obs_names.astype(str)
    adata.obs.index.name = "spot_id"

    # normalize the structure of adata.uns['spatial'] to use sample_name as the key
    if "spatial" in adata.uns and isinstance(adata.uns["spatial"], dict) and len(adata.uns["spatial"]) > 0:
        original_key = list(adata.uns["spatial"].keys())[0]
        adata.uns["spatial"] = {sample_name: adata.uns["spatial"][original_key]}

    return adata


def remove_colorbar_axes(fig, aspect_threshold=0.08):
    """
    Remove colorbar-like axes from a matplotlib Figure based on aspect ratio.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    aspect_threshold : float
        Axes with width/height lower than this are considered colorbars and removed.
    """
    to_remove = []
    for ax in fig.axes:
        bbox = ax.get_position()
        w, h = bbox.width, bbox.height
        if h > 0 and (w / h) < aspect_threshold:
            to_remove.append(ax)
    for ax in to_remove:
        fig.delaxes(ax)


# ---- read data ----
slides = [read_and_qc(sn, path=sp_data_folder) for sn in sample_names]
name_to_slide = {sn: ad for sn, ad in zip(sample_names, slides)}
slides_ordered = [name_to_slide[sn] for sn in sample_names]

# ---- figure layout (width = 600 mm) ----
n_rows = 2 + len(marker_genes)   # 1 row for H&E background, 1 row for metric, rest for genes
n_cols = len(sample_names)

FIG_WIDTH_MM = 600.0
width_in = FIG_WIDTH_MM / 25.4
cell_w_in = width_in / n_cols
cell_h_in = cell_w_in * 1.05
height_in = cell_h_in * n_rows

fig, axes = plt.subplots(n_rows, n_cols, figsize=(width_in, height_in))
plt.subplots_adjust(wspace=0.2, hspace=1.0)

# ---- row 1: H&E background (no spots) ----
for j, (ad, sid) in enumerate(zip(slides_ordered, sample_names)):
    lib_id = get_library_id(ad)
    ax = axes[0, j]
    sc.pl.spatial(
        ad,
        library_id=lib_id,
        color=None,
        size=0,
        alpha_img=alpha_img,
        img_key=img_key,
        show=False,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(sid)

# ---- row 2: QC metric per sample (individual scaling) ----
for j, (ad, sid) in enumerate(zip(slides_ordered, sample_names)):
    lib_id = get_library_id(ad)
    ax = axes[1, j]
    sc.pl.spatial(
        ad,
        library_id=lib_id,
        color=metric_name,
        size=size_val,
        alpha_img=alpha_img,
        img_key=img_key,
        vmin=0,
        vmax="p99.2",
        show=False,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"{sid} - {metric_name}")

# ---- remaining rows: one gene per row (per-sample scaling) ----
for i, gene in enumerate(marker_genes, start=2):
    for j, (ad, sid) in enumerate(zip(slides_ordered, sample_names)):
        lib_id = get_library_id(ad)
        ax = axes[i, j]
        if gene in ad.var_names:
            sc.pl.spatial(
                ad,
                library_id=lib_id,
                color=gene,
                size=size_val,
                alpha_img=alpha_img,
                img_key=img_key,
                vmin=0,
                vmax="p99.2",
                show=False,
                ax=ax,
            )
            ax.set_title(f"{sid} - {gene}")
        else:
            sc.pl.spatial(
                ad,
                library_id=lib_id,
                color=None,
                size=0,
                alpha_img=alpha_img,
                img_key=img_key,
                show=False,
                ax=ax,
            )
            ax.set_title(f"{sid} - {gene} (absent)")
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticks([])
        ax.set_yticks([])

# ---- clean up and export ----
remove_colorbar_axes(fig)
fig.savefig("figS6.pdf", bbox_inches="tight")
