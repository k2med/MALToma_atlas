# ---- packages ----
import os
from pathlib import Path
import logging

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import cell2location

import matplotlib as mpl
import matplotlib.pyplot as plt

from pathlib import Path
from cell2location.utils.filtering import filter_genes
from cell2location import run_colocation
from cell2location.plt import plot_spatial
from matplotlib import rcParams

# ---- global config ----
rcParams["pdf.fonttype"] = 42
scvi.settings.seed = 0
sc.settings.verbosity = 2

DATA_DIR = Path("../../data")
SP_DATA_DIR = DATA_DIR / "spatial_rnaseq"
REF_H5AD = Path("./snrna_reference.h5ad") # obtained from fig2e.R
REF_MODEL_DIR = Path("./reference_signatures")
ST_OUT_DIR = Path("./cell2location_map")

SAMPLE_NAMES = [
    "S01_ST", "S02_ST", "S04_ST", "S06_ST",
    "S31_ST", "S32_ST", "S71_ST", "S87_ST", "S89_ST"
]

REF_MODEL_DIR.mkdir(parents=True, exist_ok=True)
ST_OUT_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s - %(message)s")
logger = logging.getLogger("c2l")


# ---- utilities ----
def read_and_qc(sample_name: str, base_path: Path = SP_DATA_DIR):
    """Read one 10x Visium experiment into an AnnData object and compute QC metrics."""
    sample_dir = base_path / sample_name
    adata = sc.read_visium(sample_dir, count_file="filtered_feature_bc_matrix.h5", load_images=True)
    adata.var_names_make_unique()

    adata.obs["sample"] = str(sample_name)
    adata.var["SYMBOL"] = adata.var_names

    sc.pp.calculate_qc_metrics(adata, inplace=True)

    adata.var["MT"] = [g.startswith("MT-") for g in adata.var["SYMBOL"]]
    mt_mask = adata.var["MT"].values
    mt_sum = np.array(adata[:, mt_mask].X.sum(axis=1)).ravel()
    adata.obs["MT_frac"] = mt_sum / adata.obs["total_counts"].replace(0, np.nan)

    adata.obs.index = adata.obs["sample"].astype(str) + "_" + adata.obs_names.astype(str)
    adata.obs.index.name = "spot_id"

    if "spatial" in adata.uns and len(adata.uns["spatial"]) > 0:
        first_key = list(adata.uns["spatial"].keys())[0]
        adata.uns["spatial"] = {sample_name: adata.uns["spatial"][first_key]}

    return adata


def select_slide(adata, sample, s_col="sample"):
    """Select one slide from a combined spatial AnnData object."""
    slide = adata[adata.obs[s_col].isin([sample]), :]
    keys = list(slide.uns["spatial"].keys())
    spatial_key = np.array(keys)[[sample in k for k in keys]][0]
    slide.uns["spatial"] = {spatial_key: slide.uns["spatial"][spatial_key]}
    return slide


# ---- reference model ----
adata_ref = sc.read_h5ad(REF_H5AD)

selected = filter_genes(
    adata_ref,
    cell_count_cutoff=5,
    cell_percentage_cutoff2=0.03,
    nonz_mean_cutoff=1.12,
)
adata_ref = adata_ref[:, selected].copy()

cell2location.models.RegressionModel.setup_anndata(
    adata=adata_ref,
    batch_key="orig.ident",
    labels_key="fine_cell_type",
)

from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)
mod.view_anndata_setup()
mod.train(max_epochs=250)

adata_ref = mod.export_posterior(
    adata_ref,
    sample_kwargs={"num_samples": 1000, "batch_size": 2500, "use_gpu": True},
)
adata_ref = mod.export_posterior(
    adata_ref,
    use_quantiles=True,
    add_to_varm=["q05", "q50", "q95", "q0001"],
    sample_kwargs={"batch_size": 2500},
)

mod.save(REF_MODEL_DIR, overwrite=True)

adata_ref.__dict__["_raw"].__dict__["_var"] = adata_ref.__dict__["_raw"].__dict__["_var"].rename(
    columns={"_index": "features"}
)

if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
    inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
        [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
    ].copy()
else:
    inf_aver = adata_ref.var[
        [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
    ].copy()
inf_aver.columns = adata_ref.uns["mod"]["factor_names"]


# ---- spatial data ----
slides = [read_and_qc(name, base_path=SP_DATA_DIR) for name in SAMPLE_NAMES]

adata_total = slides[0].concatenate(
    slides[1:],
    batch_key="sample",
    uns_merge="unique",
    batch_categories=SAMPLE_NAMES,
    index_unique=None,
)

adata_total.obsm["MT"] = adata_total[:, adata_total.var["MT"].values].X.toarray()
adata_total = adata_total[:, ~adata_total.var["MT"].values]

intersect = np.intersect1d(adata_total.var_names, inf_aver.index)
adata = adata_total[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

cell2location.models.Cell2location.setup_anndata(adata=adata, batch_key="sample")

mod = cell2location.models.Cell2location(
    adata,
    cell_state_df=inf_aver,
    N_cells_per_location=30,
    detection_alpha=20,
)
mod.view_anndata_setup()
mod.train(
    max_epochs=30000,
    batch_size=None,
    train_size=1,
    use_gpu=True,
)

mod.plot_history(1000)
plt.legend(labels=["full data training"])

adata = mod.export_posterior(
    adata,
    sample_kwargs={"num_samples": 1000, "batch_size": mod.adata.n_obs, "use_gpu": True},
)

pd.DataFrame(adata.obsm['q05_cell_abundance_w_sf']).to_csv(ST_OUT_DIR / "st_cell2location_res.csv")
adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']

adata.write(ST_OUT_DIR / "st_cell2location_res_q05.h5ad")

FIG_DIR = Path(".")
FIG_DIR.mkdir(parents=True, exist_ok=True)

# colocation
RESULTS_DIR = Path("./colocation")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
run_name = RESULTS_DIR / "cell2location_map"

res_dict, adata = run_colocation(
    adata=adata,
    model_name="CoLocatedGroupsSklearnNMF",
    train_args={
        "n_fact": np.arange(2, 10),
        "sample_name_col": "sample",
        "n_restarts": 3,
    },
    model_kwargs={"alpha": 0.01, "init": "random", "nmf_kwd_args": {"tol": 1e-6}},
    export_args={"path": str(run_name / "CoLocatedComb")},
)

# Fig. 7c (HE background only)
slide = select_slide(adata, "S01_ST")
sc.pl.spatial(
    adata=slide,
    img_key="hires",
    color=None,
    show=False,
    size=1.5,
)
plt.savefig(FIG_DIR / "fig7c_HE.pdf", dpi=300, bbox_inches="tight")
# Other panels in fig7c were obtained from:
# "./colocation/cell2location_map/CoLocatedComb/CoLocatedGroupsSklearnNMF_30515locations_25factors/spatial"

# fig7e (left/right): S06_ST, up to 6 clusters
clust_labels = ["B_Mal", "CD8T_HAVCR2", "CD4T_FOXP3", "cDC_XCR1", "cDC_CLEC10A", "Fibro_CCL19"]
clust_col = [str(x) for x in clust_labels]  # keep compatibility if column names differ
slide = select_slide(adata, "S06_ST")

coords = slide.obsm["spatial"]
sample_id = next(iter(slide.uns["spatial"]))

scalef = slide.uns["spatial"][sample_id]["scalefactors"]["tissue_hires_scalef"]
img = slide.uns["spatial"][sample_id]["images"]["hires"]
h, w = img.shape[:2]

crop_x = (
    max(0, coords[:, 0].min() * scalef) - 10,
    min(coords[:, 0].max() * scalef, w) + 10,
)
crop_y = (
    max(0, coords[:, 1].min() * scalef) - 10,
    min(coords[:, 1].max() * scalef, h) + 10,
)

plt.figure()
with mpl.rc_context({"figure.figsize": (10, 10)}):
    _ = plot_spatial(
        adata=slide,
        color=clust_col,
        labels=clust_labels,
        style="fast",
        max_color_quantile=0.992,
        circle_diameter=6,
        crop_x=crop_x,
        crop_y=crop_y,
        colorbar_position="right",
    )
plt.savefig(FIG_DIR / "fig7e_left.pdf", dpi=300, bbox_inches="tight")
plt.close()

plt.figure()
with mpl.rc_context({"axes.facecolor": "black", "figure.figsize": [4.5, 5]}):
    sc.pl.spatial(
        slide,
        cmap="magma",
        color=clust_labels,
        ncols=4,
        size=1.3,
        img_key="hires",
        vmin=0,
        vmax="p99.2",
        show=False,
        alpha_img=0,
    )
plt.savefig(FIG_DIR / "fig7e_right.pdf", dpi=300, bbox_inches="tight")
plt.close()

# fig7f (left/right): S31_ST, up to 6 clusters
clust_labels = ["B_Mal", "Neu_G0S2", "Mac_MARCO", "Endo_CXCL12", "Endo_ACKR1", "Peri_PDGFRB"]
clust_col = [str(x) for x in clust_labels]
slide = select_slide(adata, "S31_ST")

coords = slide.obsm["spatial"]
sample_id = next(iter(slide.uns["spatial"]))

scalef = slide.uns["spatial"][sample_id]["scalefactors"]["tissue_hires_scalef"]
img = slide.uns["spatial"][sample_id]["images"]["hires"]
h, w = img.shape[:2]

crop_x = (
    max(0, coords[:, 0].min() * scalef) - 10,
    min(coords[:, 0].max() * scalef, w) + 10,
)
crop_y = (
    max(0, coords[:, 1].min() * scalef) - 10,
    min(coords[:, 1].max() * scalef, h) + 10,
)

plt.figure()
with mpl.rc_context({"figure.figsize": (10, 10)}):
    _ = plot_spatial(
        adata=slide,
        color=clust_col,
        labels=clust_labels,
        style="fast",
        max_color_quantile=0.992,
        circle_diameter=6,
        crop_x=crop_x,
        crop_y=crop_y,
        colorbar_position="right",
    )
plt.savefig(FIG_DIR / "fig7f_left.pdf", dpi=300, bbox_inches="tight")
plt.close()

plt.figure()
with mpl.rc_context({"axes.facecolor": "black", "figure.figsize": [4.5, 5]}):
    sc.pl.spatial(
        slide,
        cmap="magma",
        color=clust_labels,
        ncols=4,
        size=1.3,
        img_key="hires",
        vmin=0,
        vmax="p99.2",
        show=False,
        alpha_img=0,
    )
plt.savefig(FIG_DIR / "fig7f_right.pdf", dpi=300, bbox_inches="tight")
plt.close()