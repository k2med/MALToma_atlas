# ---- run_misty_seurat: one-shot wrapper to build & run a MISTy workflow ----
run_misty_seurat <- function(
    visium.slide,      # Seurat object with Visium data
    view.assays,       # named list: assay per view
    view.features = NULL, # named list: features per view (NULL = all)
    view.types,        # named list: view type per assay ("intra"/"para"/"juxta")
    view.params,       # named list: parameter per view (NULL or value)
    spot.ids = NULL,   # vector of spot IDs (NULL = all)
    out.alias = "results" # output folder alias
) {
  mistyR::clear_cache()
  
  # geometry (Visium tissue positions)
  geometry <- GetTissueCoordinates(
    visium.slide, cols = c("row", "col"), scale = NULL
  )
  
  # extract per-view data from Seurat
  view.data <- map(
    view.assays, extract_seurat_data,
    geometry = geometry, visium.slide = visium.slide
  )
  
  # build & run pipeline
  build_misty_pipeline(
    view.data     = view.data,
    view.features = view.features,
    view.types    = view.types,
    view.params   = view.params,
    geometry      = geometry,
    spot.ids      = spot.ids,
    out.alias     = out.alias
  )
}

# ---- extract_seurat_data: matrix from assay, aligned to geometry ----
extract_seurat_data <- function(visium.slide, assay, geometry) {
  print(assay)
  data <- GetAssayData(visium.slide, assay = assay) %>%
    as.matrix() %>% t() %>% as_tibble(rownames = NA)
  
  return(
    data %>% slice(match(rownames(.), rownames(geometry)))
  )
}

# ---- filter_data_features: keep only selected features (or all if NULL) ----
filter_data_features <- function(data, features) {
  if (is.null(features)) features <- colnames(data)
  
  return(
    data %>%
      rownames_to_column() %>%
      select(rowname, all_of(features)) %>%
      rename_with(make.names) %>%
      column_to_rownames()
  )
}

# ---- create_default_views: build a view by type/param and subset spots ----
create_default_views <- function(
    data, view.type, view.param, view.name, spot.ids, geometry
) {
  mistyR::clear_cache()
  view.data.init <- create_initial_view(data)
  
  if (!(view.type %in% c("intra", "para", "juxta"))) view.type <- "intra"
  
  if (view.type == "intra") {
    data.red <- view.data.init[["intraview"]]$data %>%
      rownames_to_column() %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
    
  } else if (view.type == "para") {
    view.data.tmp <- view.data.init %>% add_paraview(geometry, l = view.param)
    data.ix  <- paste0("paraview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
    
  } else if (view.type == "juxta") {
    view.data.tmp <- view.data.init %>% add_juxtaview(
      positions = geometry, neighbor.thr = view.param
    )
    data.ix  <- paste0("juxtaview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  }
  
  misty.view <- if (is.null(view.param)) {
    create_view(paste0(view.name), data.red)
  } else {
    create_view(paste0(view.name, "_", view.param), data.red)
  }
  
  return(misty.view)
}

# ---- build_misty_pipeline: assemble views and run MISTy ----
build_misty_pipeline <- function(
    view.data, view.features, view.types, view.params,
    geometry, spot.ids = NULL, out.alias = "default"
) {
  if (is.null(spot.ids)) spot.ids <- rownames(view.data[[1]])
  
  # feature filtering per view
  view.data.filt <- map2(view.data, view.features, filter_data_features)
  
  # main (intra) view
  views.main <- create_initial_view(
    view.data.filt[[1]] %>%
      rownames_to_column() %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  )
  
  # additional views (para/juxta/extra intra)
  view.names <- names(view.data.filt)
  all.views <- pmap(
    list(
      view.data.filt[-1],
      view.types[-1],
      view.params[-1],
      view.names[-1]
    ),
    create_default_views,
    spot.ids = spot.ids,
    geometry = geometry
  )
  
  pline.views <- add_views(views.main, unlist(all.views, recursive = FALSE))
  
  # run MISTy
  run_misty(pline.views, out.alias, cached = FALSE)
}