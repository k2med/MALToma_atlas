# ---- packages ----
suppressPackageStartupMessages({
  library(RColorBrewer)
  library(viridis)
})

# ---- subtype palettes ----
subtype_col1 <- c(
  "LME_1" = "#fadcb4",
  "LME_2" = "#dea3a2",
  "LME_3" = "#a7b9d7"
)

subtype_col2 <- c(
  "LME_1" = "#e3b87f",
  "LME_2" = "#b57979",
  "LME_3" = "#576fa0"
)

# ---- NMF cluster palettes ----
NMF_cluster_col1 <- c(
  "#fadcb4", "#dea3a2", "#a7b9d7",
  "#c0686c", "#c7b4a2", "#6f9bb9", "#ac9db8"
)

NMF_cluster_col2 <- c(
  "#6cc1a4", "#f78c65", "#8f9ec9",
  "#e588c1", "#a4d955", "#f8db3f", "#f8db3f"
)

# ---- cohort / group / clinical palettes ----
cohort_col <- c(
  "WCH"       = "#49c1ad",
  "Shi_2022"  = "#bf8ca9",
  "Wolf_2023" = "#e3bc7b"
)

group_col <- c(
  "Tumor"  = "#cdb783",
  "Normal" = "#52a09e"
)

sex_col <- c(
  "Female" = "#c5e4ae",
  "Male"   = "#f5d8a8",
  "NA"     = "#d0d0d0"
)

age_col <- c(
  ">=80" = "#006837",
  "70-79" = "#41ab5d",
  "60-69" = "#78c679",
  "50-59" = "#addd8e",
  "40-49" = "#d9f0a3",
  "<40"   = "#ffffe5",
  "NA"    = "#d0d0d0"
)

age_group_col <- c(
  "<60"  = "#d9f0a3",
  ">=60" = "#006837"
)

site_col <- c(
  "Ocular_adnexa"   = "#CDBED2",
  "Stomach"         = "#EACFB9",
  "Lung"            = "#49C1AD",
  "Intestine"       = "#b5d2e2",
  "Nasopharynx"     = "#E7B922",
  "Hypopharynx"     = "#fdb462",
  "Paranasal_sinus" = "#b3de69",
  "Liver"           = "#80b1d3",
  "Salivary_gland"  = "#B1E4C2",
  "Lymph_node"      = "#F4D8A8",
  "Others"          = "#a8b4d6"
)

stage_col <- c(
  "I"  = "#efedf5",
  "II" = "#bcbddc",
  "III"= "#807dba",
  "IV" = "#6a51a3",
  "NA" = "#d0d0d0"
)

stage_group_col <- c(
  "Stage I_III" = "#efedf5",
  "Stage IV"    = "#6a51a3"
)

fusion_col <- c(
  "Neg" = "#80b1d3",
  "Pos" = "#fb8072",
  "NA"  = "#d0d0d0"
)

# ---- continuous scales ----
scale_col1 <- brewer.pal(n = 11, name = "RdYlBu")
scale_col2 <- rev(hcl.colors(n = 10, palette = "Earth"))

# ---- heatmaps ----
heatmap_col1 <- colorRampPalette(c(
  "#38629D", "#357BA2", "#3492A8", "#38AAAC", "#49C1AD", "#78D6AE",
  "#B1E4C2", "#DEF5E5", "white", "#F0EAF9", "#E8CDE3", "#DFABC9", "#D485AA",
  "#BC6892", "#93658F", "#71608C", "#534C7A"
))(20)

heatmap_col2 <- colorRampPalette(rev(brewer.pal(5, "RdYlBu")))(20)

heatmap_col3 <- c(
  "black", "#006aff", "#0068fe", "#0066fc", "#0063fa", "#0061f8", "#005ff6", "#005df4",
  "#005af2", "#0058f0", "#0056ee", "#0054ec", "#0052ea", "#004fe8", "#004de6", "#004be4",
  "#0049e2", "#0046e0", "#0045de", "#0044dd", "#0042db", "#0041d9", "#0040d8", "#003fd6",
  "#003ed4", "#003dd3", "#003cd1", "#003bcf", "#003ace", "#0039cc", "#0038ca", "#0037c9",
  "#0036c7", "#0034c5", "#0033c4", "#0131ba", "#0c2eae", "#132ca1", "#172995", "#192688",
  "#1b247c", "#1c2170", "#1c1f65", "#1b1c59", "#1a1a4e", "#191743", "#171539", "#15122e",
  "#130e25", "#100a1b", "#080511", "#000000", "#121205", "#1d1e0a", "#26290e", "#30360f",
  "#3b4210", "#464f0f", "#515c0e", "#5c6a0c", "#677807", "#738601", "#7f9500", "#8ba400",
  "#97b300", "#a4c200", "#b1d200", "#bde100", "#c6ec00", "#c8ee00", "#c9ef00", "#caf000",
  "#cbf100", "#ccf300", "#cef400", "#cff500", "#d0f600", "#d1f800", "#d2f900", "#d4fa00",
  "#d5fc00", "#d6fd00", "#d7fe00", "#d8ff00", "#daff00", "#dbff00", "#ddff00", "#dfff00",
  "#e0ff00", "#e2ff00", "#e4ff00", "#e6ff00", "#e7ff00", "#e9ff00", "#ebff00", "#edff02",
  "#eeff07", "#f0ff0b", "#f2ff0e", "#f4ff12", "#f5ff14", "#f7ff17"
)

heatmap_col4 <- c(
  colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10),
  rev(magma(323, begin = 0.18))
)

heatmap_col5 <- colorRampPalette(rev(brewer.pal(5, "RdBu")))(20)

# ---- cell-type palettes ----
coarse_cell_type_col <- c(
  "B_cell"        = "#9bbcb7",
  "T_cell"        = "#9fccdd",
  "Myeloid_cell"  = "#e4d173",
  "Stromal_cell"  = "#c4a9db",
  "Hepatocyte"    = "#d7e4f4",
  "Alveolar_cell" = "#e8d7ba"
)

fine_cell_type_col <- c(
  "B_Mal"        = "#d9dec5",
  "B_Nor"        = "#8acac1",
  
  "CD8T_CCL5"    = "#ec8c70",
  "CD8T_GNLY"    = "#dbd08c",
  "CD8T_HAVCR2"  = "#e9b77a",
  "CD4T_FOXP3"   = "#4690a9",
  "CD4T_CXCL13"  = "#bebada",
  "CD4T_IL7R"    = "#a0d09d",
  
  "Mac_SPARC"    = "#dfd179",
  "Mac_STAB1"    = "#bda988",
  "Mac_MARCO"    = "#c1c7b1",
  "Mac_CHIT1"    = "#a77668",
  "cDC_XCR1"     = "#b7d0db",
  "cDC_CLEC10A"  = "#72c0c7",
  "pDC_CLEC4C"   = "#87b67f",
  "Neu_G0S2"     = "#e0eab6",
  "Mast_CPA3"    = "#f8d0cc",
  
  "Endo_ACKR1"   = "#dfdfc1",
  "Endo_CXCL12"  = "#aed491",
  "Endo_CCL21"   = "#f4e4b1",
  "Endo_HPGD"    = "#a2988f",
  "Peri_PDGFRB"  = "#646669",
  "Fibro_PI16"   = "#7fa2b7",
  "Fibro_CCL19"  = "#8cc0e1",
  "Fibro_LRRC15" = "#b2c1c7",
  "Fibro_TCF21"  = "#d6d2e9",
  "FDC_CR2"      = "#deafbf",
  
  "Hepa_ALB"     = "#d7e4f4",
  "AT1_AGER"     = "#c8b88a",
  "AT2_SFTPC"    = "#59888f"
)