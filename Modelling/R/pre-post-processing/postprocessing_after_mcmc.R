###############################################################################
# R postprocessing, between step 1 (training) and 2 (predictions) in TreePPL
###############################################################################

###########################################
# Functions needed
###########################################

post_processing_single_mcmc <- function(post_data) {
  samples <- post_data$samples$`__data__` #  samples <- post_data$samples[[1]]$`__data__`
  n <- length(samples$k)  # number of MCMC samples
  post_table <- vector("list", n)
  
  for (i in 1:n) {
    # extract k value (only one)
    k_val <- samples$k[[i]]
    
    # extract c values (length 15 vector)
    c_vals <- samples$c[[i]]
    names(c_vals) <- paste0("c", 1:15)
    
    # extract theta 
    theta_val <- samples$theta[[i]]
    
    # create single row data frame
    row_df <- data.frame(
      theta_dist = theta_val,
      k_bio = k_val,
      t(c_vals),
      weights = 1 / n,
      row.names = NULL,
      check.names = FALSE
    )
    
    post_table[[i]] <- row_df
  }
  
  # Combine into one data framesource
  post_df <- do.call(rbind, post_table)

  # Determine thinning indices
  N <- 500 
  if (n < N) {
    stop(paste0("Not enough samples to thin to ", N))
  }
  keep_idx <- sample(unique(seq(1, n, length.out=n)), N)

  # Thin 
  post_df_thinned <- post_df[keep_idx, ]
  
  return(post_df_thinned)
}

post_processing_comb_single_mcmc <- function(post_data) {
  samples <- post_data$samples$`__data__` #  samples <- post_data$samples[[1]]$`__data__`
  n <- length(samples$k_bio)  # number of MCMC samples
  post_table <- vector("list", n)
  
  for (i in 1:n) {
    # extract k value (only one)
    k_val <- samples$k_bio[[i]]
    
    # extract c values (length 15 vector)
    c_vals <- samples$c[[i]]
    names(c_vals) <- paste0("c", 1:15)
    
    # extract theta 
    theta_val <- samples$theta[[i]]
    
    # create single row data frame
    row_df <- data.frame(
      theta_dist = theta_val,
      k_bio = k_val,
      t(c_vals),
      weights = 1 / n,
      row.names = NULL,
      check.names = FALSE
    )
    
    post_table[[i]] <- row_df
  }
  
  # Combine into one data frame
  post_df <- do.call(rbind, post_table)
  
  # Determine thinning indices
  N <- 500 
  if (n < N) {
    stop(paste0("Not enough samples to thin to ", N))
  }
  keep_idx <- sample(unique(seq(1, n, length.out=n)), N)

  # Thin 
  post_df_thinned <- post_df[keep_idx, ]
  
  return(post_df_thinned)
}

post_processing_comb_6k_mcmc <- function(post_data) {
  samples <- post_data$samples$`__data__` #  samples <- post_data$samples[[1]]$`__data__`
  n <- length(samples$k)  # number of MCMC samples
  post_table <- vector("list", n)
  
  for (i in 1:n) {
    # extract k values (length 6 vector)
    k_vals <- samples$k[[i]]
    names(k_vals) <- paste0("k", 1:6)
    
    # extract c values (length 15 vector)
    c_vals <- samples$c[[i]]
    names(c_vals) <- paste0("c", 1:15)
    
    # extract theta (you might need to adjust this key)
    theta_val <- samples$theta[[i]]
    
    # create single row data frame
    row_df <- data.frame(
      theta_dist = theta_val,
      t(k_vals),
      t(c_vals),
      weights = 1 / n,
      row.names = NULL,
      check.names = FALSE
    )
    
    post_table[[i]] <- row_df
  }
  
  # Combine into one data frame
  post_df <- do.call(rbind, post_table)
  
  # Determine thinning indices
  N <- 500 
  if (n < N) {
    stop(paste0("Not enough samples to thin to ", N))
  }
  keep_idx <- sample(unique(seq(1, n, length.out=n)), N)

  # Thin 
  post_df_thinned <- post_df[keep_idx, ]
  
  return(post_df_thinned)
}

post_processing_comb_8k_mcmc <- function(post_data) {
  samples <- post_data$samples$`__data__` #  samples <- post_data$samples[[1]]$`__data__`
  n <- length(samples$k)  # number of MCMC samples
  post_table <- vector("list", n)
  
  for (i in 1:n) {
    # Get 6 selected k values (3rd to 8th)
    k_vals <- samples$k[[i]][3:8]
    names(k_vals) <- paste0("k", 1:6)
    
    # Get c values
    c_vals <- samples$c[[i]]
    names(c_vals) <- paste0("c", 1:15)
    
    # Extract theta value (adjust key if needed)
    theta_val <- samples$theta[[i]]  # Replace with correct param if different
    
    # Combine into a row
    row_df <- data.frame(
      theta_dist = theta_val,
      t(k_vals),
      t(c_vals),
      weights = 1 / n,
      row.names = NULL,
      check.names = FALSE
    )
    
    post_table[[i]] <- row_df
  }

  # Combine into one data frame
  post_df <- do.call(rbind, post_table)
  
  # Determine thinning indices
  N <- 500 
  if (n < N) {
    stop(paste0("Not enough samples to thin to ", N))
  }
  keep_idx <- sample(unique(seq(1, n, length.out=n)), N)

  # Thin 
  post_df_thinned <- post_df[keep_idx, ]
  
  return(post_df_thinned)
}

output_for_pred_single_mcmc <- function(post_table, prior_data, output_file) {
  # Extract theta_list
  theta_list <- as.numeric(post_table$theta_dist)
  
  # Extract and transpose matrices: [param x sample]
  k_bio_dist <- as.numeric(post_table$k_bio)
  c_matrix <- t(as.matrix(post_table[, 3:17]))     # c1–c15

  # Convert c to list of vectors, then remove row names to avoid "1", "2", etc.
  c_dist <- unname(split(c_matrix, row(c_matrix)))

  # Extract weights
  weights <- as.numeric(post_table$weights)

  # Clean up prior data
  prior_data <- as.list(prior_data)
  prior_data[c("k_1", "k_2", "theta_1", "theta_2", "tau_1")] <- NULL

  # Construct final list
  combined_list <- list(
    k_bio_dist = k_bio_dist,
    c_dist = c_dist,
    theta_list = theta_list,
    weights = weights
  )
  combined_list <- c(combined_list, prior_data)

  # Write to JSON
  json_str <- jsonlite::toJSON(combined_list, pretty = FALSE, auto_unbox = TRUE, digits = 16)
  writeLines(json_str, con = output_file)
}

output_for_pred_comb_mcmc <- function(post_table, prior_data, output_file) {
  # Extract theta_list
  theta_list <- as.numeric(post_table$theta_dist)
  
  # Extract and transpose matrices: [param x sample]
  k_matrix <- t(as.matrix(post_table[, 2:7]))      # k1–k6
  c_matrix <- t(as.matrix(post_table[, 8:22]))     # c1–c15

  # Convert to list of vectors, then remove row names to avoid "1", "2", etc.
  k_bio_dist <- unname(split(k_matrix, row(k_matrix)))
  c_dist <- unname(split(c_matrix, row(c_matrix)))

  # Extract weights
  weights <- as.numeric(post_table$weights)

  # Clean up prior data
  prior_data <- as.list(prior_data)
  prior_data[c("k_1", "k_2", "theta_1", "theta_2", "tau_1")] <- NULL

  # Construct final list
  combined_list <- list(
    k_bio_dist = k_bio_dist,
    c_dist = c_dist,
    theta_list = theta_list,
    weights = weights
  )
  combined_list <- c(combined_list, prior_data)

  # Write to JSON
  json_str <- jsonlite::toJSON(combined_list, pretty = FALSE, auto_unbox = TRUE, digits = 16)
  writeLines(json_str, con = output_file)
}