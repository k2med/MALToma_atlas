distance_comparison <- function(distances, proportions, topic1, topic2, t1_thresh, t2_thresh) {
  
  # ---- filter topic1 spots above threshold ----
  filtered_props_1 <- proportions[proportions[, topic1] > t1_thresh, ]
  rownames_holder  <- rownames(filtered_props_1)
  filtered_props_1 <- data.frame(filtered_props_1[, c(topic1)])
  rownames(filtered_props_1) <- rownames_holder
  rm(rownames_holder)
  colnames(filtered_props_1) <- c(topic1)
  print(filtered_props_1)
  
  # ---- filter topic2 spots above threshold ----
  filtered_props_2 <- proportions[proportions[, topic2] > t2_thresh, ]
  rownames_holder  <- rownames(filtered_props_2)
  filtered_props_2 <- data.frame(filtered_props_2[, c(topic2)])
  rownames(filtered_props_2) <- rownames_holder
  rm(rownames_holder)
  colnames(filtered_props_2) <- c(topic2)
  print(filtered_props_2)
  
  # ---- restrict distance matrix to topic1 (rows) vs topic2 (columns) ----
  dist_filtered <- distances
  dist_filtered <- dist_filtered[rownames(dist_filtered) %in% rownames(filtered_props_1), ]
  dist_filtered <- dist_filtered[, colnames(dist_filtered) %in% rownames(filtered_props_2)]
  print(dist_filtered)
  
  # ---- compute minimum distance from each topic1 spot to topic2 set ----
  min_df <- apply(dist_filtered, 1, function(x) min(x))
  
  # ---- format as dataframe ----
  min_df <- as.data.frame(min_df)
  colnames(min_df) <- c("min")
  
  comp_str <- paste0(topic1, " vs. ", topic2)
  min_df$barcodes <- comp_str
  print(min_df)
  
  return(min_df)
}


min_distance_to_topic <- function(distances, proportions, topic, threshold) {
  
  # ---- identify spots with topic proportion above threshold ----
  topic_high_spots <- rownames(proportions[proportions[, topic] > threshold, , drop = FALSE])
  
  # ---- restrict distance matrix to only include topic-high spots (columns) ----
  distances <- distances[, colnames(distances) %in% topic_high_spots, drop = FALSE]
  
  # ---- compute minimum distance from each spot to the set of topic-high spots ----
  min_dist <- apply(distances, 1, function(x) {
    if (length(x) == 0) return(NA)
    return(min(x, na.rm = TRUE))
  })
  
  # ---- format as dataframe ----
  min_dist_df <- data.frame(
    spot = names(min_dist),
    min_distance_to_topic = min_dist,
    reference_topic = topic
  )
  
  return(min_dist_df)
}
