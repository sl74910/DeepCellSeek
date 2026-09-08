#' Cell Type Annotation using LLMs
#' 
#' @param input Input data for cell type annotation (data frame or named list)
#' @param tissuename Tissue name (e.g., "Brain", "Lung")
#' @param species Species name (e.g., "Human", "Mouse") 
#' @param model Model to use for annotation
#' @param topgenenumber Number of top genes to consider
#' @param api_key API key for the LLM service
#' @param wait_indefinitely Whether to wait indefinitely for the LLM response
#'   instead of using the default request timeout (default: FALSE)
#' @param allowed_cell_types Optional character vector, data frame with a
#'   `cell_type` column, or RDS path containing the only labels the model may use
#' @return Cell type annotation results
#' @export
llm_celltype <- function(input, tissuename = NULL, species = "Human", model = "deepseek-v4-flash",
                         topgenenumber = 10, api_key = NULL,
                         allowed_cell_types = NULL, wait_indefinitely = FALSE) {

  if (!requireNamespace("glue", quietly = TRUE)) {
    stop("Package 'glue' is required. Please install it with: install.packages('glue')")
  }

  if (!is_model_supported(model)) {
    stop("❌ Model not supported: ", model, "\nSupported models: ",
         paste(unlist(get_supported_models()), collapse = ", "))
  }

  if (!is.logical(wait_indefinitely) || length(wait_indefinitely) != 1L ||
      is.na(wait_indefinitely)) {
    stop("'wait_indefinitely' must be a single TRUE or FALSE value")
  }

  config <- get_model_config(model)
  if (is.null(api_key)) {
    api_key <- resolve_api_key(config)
    if (api_key == "") {
      cat("📝 Note:", config$provider, "API key not found: returning the prompt itself.\n")
      API.flag <- 0
    } else {
      cat("🔑 Note:", config$provider, "API key found: proceeding with annotation.\n")
      API.flag <- 1
    }
  } else {
    API.flag <- 1
  }

  processed_input <- process_input_data(input, topgenenumber)
  allowed_cell_types <- resolve_allowed_cell_types(allowed_cell_types)
  if (!API.flag) {
    marker_data <- paste0(names(processed_input), ':', unlist(processed_input), collapse = "\n")
    message <- glue::glue("Identify cell types of {species} {tissuename} cells using the following markers separately for each row.
{annotation_instruction}

{marker_data}",
                         species = species,
                         tissuename = tissuename,
                         num_clusters = length(processed_input),
                         annotation_instruction = format_annotation_instruction(
                           allowed_cell_types, length(processed_input), subtype = FALSE
                         ),
                         marker_data = marker_data)
    return(message)
  }

  timeout_seconds <- if (wait_indefinitely) NULL else 1000

  cutnum <- ceiling(length(processed_input) / 30)
  if (cutnum > 1) {
    cid <- as.numeric(cut(1:length(processed_input), cutnum))
  } else {
    cid <- rep(1, length(processed_input))
  }
  
  allres <- sapply(1:cutnum, function(i) {
    id <- which(cid == i)
    retry_count <- 0
    max_retries <- 3
    
    while (retry_count < max_retries) {
      marker_data <- paste(names(processed_input)[id], ':', processed_input[id], collapse = '\n')
      
      prompt <- glue::glue("Identify cell types of {species} {tissuename} cells using the following markers separately for each row.
{annotation_instruction}

{marker_data}",
                          species = species,
                          tissuename = tissuename,
                          num_clusters = length(id),
                          annotation_instruction = format_annotation_instruction(
                            allowed_cell_types, length(id), subtype = FALSE
                          ),
                          marker_data = marker_data)

      result <- call_llm_api(model, prompt, temperature = NULL,
                             timeout_seconds = timeout_seconds, api_key = api_key)

      if (!is.null(result)) {
        res <- trimws(unlist(strsplit(result, '\n')))
        res <- res[nzchar(res)]

        if (length(res) == length(id)) {
          names(res) <- names(processed_input)[id]
          return(res)
        } else {
          cat("⚠️ Response length mismatch for batch", i, ". Expected:", length(id), "Got:", length(res), "\n")
          
          if (retry_count >= max_retries - 1) {
            valid_results <- min(length(res), length(id))
            final_res <- rep("Failed to annotate", length(id))
            if (valid_results > 0) {
              final_res[1:valid_results] <- res[1:valid_results]
            }
            names(final_res) <- names(processed_input)[id]
            return(final_res)
          }
        }
      }
      
      retry_count <- retry_count + 1
      cat("🔄 Retrying batch", i, "(attempt", retry_count + 1, ")\n")
      Sys.sleep(2)
    }
    
    cat("❌ Failed to get proper response for batch", i, "after", max_retries, "retries\n")
    res <- rep("Failed to annotate", length(id))
    names(res) <- names(processed_input)[id]
    return(res)
  }, simplify = FALSE)
  
  cat("✅ Cell type annotation completed\n")
  cat("📝 Note: It is recommended to check results for potential AI hallucinations before downstream analysis\n")
  
  result <- gsub(',$', '', unlist(allres))
  result <- sub("^[^a-zA-Z]*([a-zA-Z].*)", "\\1", result)
  result <- restrict_to_allowed_cell_types(result, allowed_cell_types)
  
  return(result)
}

#' Cell Subtype Annotation using LLMs
#' 
#' @param input Input data for cell subtype annotation (data frame or named list)
#' @param tissuename Tissue name (e.g., "Brain", "Lung")
#' @param species Species name (e.g., "Human", "Mouse")
#' @param celltypename Parent cell type name
#' @param model Model to use for annotation
#' @param topgenenumber Number of top genes to consider
#' @param api_key API key for the LLM service
#' @param allowed_cell_types Optional character vector, data frame with a
#'   `cell_type` column, or RDS path containing the only labels the model may use
#' @return Cell subtype annotation results
#' @export
llm_subcelltype <- function(input, tissuename = NULL, species = "Human", celltypename = NULL,
                           model = "deepseek-v4-flash", topgenenumber = 10, api_key = NULL,
                           allowed_cell_types = NULL) {

  if (!requireNamespace("glue", quietly = TRUE)) {
    stop("Package 'glue' is required. Please install it with: install.packages('glue')")
  }

  if (!is_model_supported(model)) {
    stop("❌ Model not supported: ", model, "\nSupported models: ", 
         paste(unlist(get_supported_models()), collapse = ", "))
  }

  config <- get_model_config(model)
  if (is.null(api_key)) {
    api_key <- resolve_api_key(config)
    if (api_key == "") {
      cat("📝 Note:", config$provider, "API key not found: returning the prompt itself.\n")
      API.flag <- 0
    } else {
      cat("🔑 Note:", config$provider, "API key found: proceeding with annotation.\n")
      API.flag <- 1
    }
  } else {
    API.flag <- 1
  }

  processed_input <- process_input_data(input, topgenenumber)
  allowed_cell_types <- resolve_allowed_cell_types(allowed_cell_types)
  if (!API.flag) {
    marker_data <- paste0(names(processed_input), ':', unlist(processed_input), collapse = "\n")
    message <- glue::glue("Identify detailed cell subtypes for {celltypename} in {species} {tissuename} cells using the following markers, provided separately for each row.
{annotation_instruction}

{marker_data}",
                         celltypename = celltypename,
                         species = species,
                         tissuename = tissuename,
                         num_clusters = length(processed_input),
                         annotation_instruction = format_annotation_instruction(
                           allowed_cell_types, length(processed_input), subtype = TRUE
                         ),
                         marker_data = marker_data)
    return(message)
  }

  cutnum <- ceiling(length(processed_input) / 30)
  if (cutnum > 1) {
    cid <- as.numeric(cut(1:length(processed_input), cutnum))
  } else {
    cid <- rep(1, length(processed_input))
  }
  
  allres <- sapply(1:cutnum, function(i) {
    id <- which(cid == i)
    retry_count <- 0
    max_retries <- 3
    
    while (retry_count < max_retries) {
      marker_data <- paste(names(processed_input)[id], ':', processed_input[id], collapse = '\n')
      
      prompt <- glue::glue("Identify detailed cell subtypes for {celltypename} in {species} {tissuename} cells using the following markers, provided separately for each row.
{annotation_instruction}

{marker_data}",
                          celltypename = celltypename,
                          species = species,
                          tissuename = tissuename,
                          num_clusters = length(id),
                          annotation_instruction = format_annotation_instruction(
                            allowed_cell_types, length(id), subtype = TRUE
                          ),
                          marker_data = marker_data)
      
      result <- call_llm_api(model, prompt, temperature = NULL, timeout_seconds = 1000, api_key = api_key)
      
      if (!is.null(result)) {
        res <- trimws(unlist(strsplit(result, '\n')))
        res <- res[nzchar(res)]
        
        if (length(res) == length(id)) {
          names(res) <- names(processed_input)[id]
          return(res)
        } else {
          cat("⚠️ Response length mismatch for batch", i, ". Expected:", length(id), "Got:", length(res), "\n")
          
          if (retry_count >= max_retries - 1) {
            valid_results <- min(length(res), length(id))
            final_res <- rep("Failed to annotate", length(id))
            if (valid_results > 0) {
              final_res[1:valid_results] <- res[1:valid_results]
            }
            names(final_res) <- names(processed_input)[id]
            return(final_res)
          }
        }
      }
      
      retry_count <- retry_count + 1
      cat("🔄 Retrying batch", i, "(attempt", retry_count + 1, ")\n")
      Sys.sleep(2)
    }
    
    cat("❌ Failed to get proper response for batch", i, "after", max_retries, "retries\n")
    res <- rep("Failed to annotate", length(id))
    names(res) <- names(processed_input)[id]
    return(res)
  }, simplify = FALSE)
  
  cat("✅ Cell subtype annotation completed\n")
  cat("📝 Note: It is recommended to check results for potential AI hallucinations before downstream analysis\n")
  
  result <- gsub(',$', '', unlist(allres))
  result <- sub("^[^a-zA-Z]*([a-zA-Z].*)", "\\1", result)
  result <- restrict_to_allowed_cell_types(result, allowed_cell_types)
  
  return(result)
}
