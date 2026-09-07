#' Ensemble Voting for Cell Type Annotation
#' 
#' This module implements ensemble voting functionality using multiple LLM models
#' with intelligent arbitration for improved accuracy.

#' Cell Type Ensemble Annotation using Multiple LLMs
#' 
#' @param input Input data for cell type annotation (data frame or named list)
#' @param tissuename Tissue name (e.g., "Brain", "Lung")
#' @param species Species name (e.g., "Human", "Mouse") 
#' @param elite_models Vector of model names to use for initial predictions
#' @param arbitrator_model Model name to use for final arbitration
#' @param topgenenumber Number of top genes to consider per cluster
#' @param allowed_cell_types Optional character vector, data frame with a
#'   `cell_type` column, or RDS path containing the only labels the model may use
#' @param parallel Whether to call models in parallel (default: TRUE)
#' @param timeout_seconds Timeout for each API call (default: 1000)
#' @return Named vector of cell type annotations
#' @export
llm_celltype_ensemble <- function(input, 
                                 tissuename, 
                                 species = "Human",
                                 elite_models = c("deepseek-v4-flash", "deepseek-v4-pro", "kimi-k2.6"),
                                 arbitrator_model = "kimi-k2.6",
                                 topgenenumber = 10,
                                 parallel = TRUE,
                                 timeout_seconds = 1000,
                                 allowed_cell_types = NULL) {
  
  cat("🗳️ Starting ensemble voting for", species, tissuename, "cell type annotation\n")
  cat("🤖 Elite models:", paste(elite_models, collapse = ", "), "\n")
  cat("⚖️ Arbitrator model:", arbitrator_model, "\n")

  processed_input <- process_input_data(input, topgenenumber)
  allowed_cell_types <- resolve_allowed_cell_types(allowed_cell_types)
  num_clusters <- length(processed_input)
  
  cat("📊 Processing", num_clusters, "cell clusters\n")

  elite_predictions <- get_elite_predictions(
    processed_input, 
    tissuename, 
    species,
    elite_models, 
    parallel, 
    timeout_seconds,
    allowed_cell_types
  )
  
  if (length(elite_predictions) == 0) {
    stop("❌ Failed to get predictions from any elite model")
  }

  arbitration_data <- prepare_arbitration_data(elite_predictions, processed_input)

  final_result <- get_arbitrator_decision(
    arbitration_data,
    tissuename,
    species, 
    arbitrator_model,
    num_clusters,
    timeout_seconds,
    names(processed_input),
    allowed_cell_types
  )
  
  if (is.null(final_result)) {
    stop("❌ Arbitrator failed to provide final decision")
  }
  
  cat("✅ Ensemble voting completed successfully\n")
  return(final_result)
}

#' Cell Subtype Ensemble Annotation using Multiple LLMs
#' 
#' @param input Input data for cell subtype annotation
#' @param tissuename Tissue name
#' @param species Species name
#' @param celltypename Parent cell type name
#' @param elite_models Vector of model names for initial predictions
#' @param arbitrator_model Model for final arbitration
#' @param topgenenumber Number of top genes to consider
#' @param allowed_cell_types Optional character vector, data frame with a
#'   `cell_type` column, or RDS path containing the only labels the model may use
#' @param parallel Whether to use parallel processing
#' @param timeout_seconds API timeout
#' @return Named vector of cell subtype annotations
#' @export
llm_subcelltype_ensemble <- function(input,
                                    tissuename,
                                    species = "Human", 
                                    celltypename,
                                    elite_models = c("deepseek-v4-flash", "deepseek-v4-pro", "kimi-k2.6"),
                                    arbitrator_model = "kimi-k2.6",
                                    topgenenumber = 10,
                                    parallel = TRUE,
                                    timeout_seconds = 1000,
                                    allowed_cell_types = NULL) {
  
  cat("🗳️ Starting ensemble voting for", species, tissuename, celltypename, "subtype annotation\n")

  processed_input <- process_input_data(input, topgenenumber)
  allowed_cell_types <- resolve_allowed_cell_types(allowed_cell_types)
  num_clusters <- length(processed_input)

  elite_predictions <- get_elite_subtype_predictions(
    processed_input,
    tissuename,
    species,
    celltypename,
    elite_models,
    parallel,
    timeout_seconds,
    allowed_cell_types
  )
  
  if (length(elite_predictions) == 0) {
    stop("❌ Failed to get predictions from any elite model")
  }

  arbitration_data <- prepare_arbitration_data(elite_predictions, processed_input)

  final_result <- get_arbitrator_subtype_decision(
    arbitration_data,
    tissuename,
    species,
    celltypename,
    arbitrator_model,
    num_clusters,
    timeout_seconds,
    names(processed_input),
    allowed_cell_types
  )
  
  if (is.null(final_result)) {
    stop("❌ Arbitrator failed to provide final decision")
  }
  
  cat("✅ Ensemble subtype voting completed successfully\n")
  return(final_result)
}

get_elite_predictions <- function(processed_input, tissuename, species, elite_models, parallel,
                                  timeout_seconds, allowed_cell_types = NULL) {

  marker_data <- paste(names(processed_input), ':', processed_input, collapse = '\n')
  
  prompt <- glue::glue("Identify cell types of {species} {tissuename} cells using the following markers separately for each row.
{annotation_instruction}

{marker_data}",
                      species = species,
                      tissuename = tissuename,
                      num_clusters = length(processed_input),
                      annotation_instruction = format_annotation_instruction(
                        allowed_cell_types, length(processed_input), subtype = FALSE
                      ),
                      marker_data = marker_data)

  processed_predictions <- list()
  
  for (model in elite_models) {
    cat("🤖 Calling elite model:", model, "\n")
    
    result <- call_llm_api(model, prompt, temperature = NULL, timeout_seconds = timeout_seconds)
    
    if (!is.null(result) && result != "") {
      lines <- trimws(unlist(strsplit(result, "\n")))
      lines <- lines[lines != ""]
      
      if (length(lines) == length(processed_input)) {
        lines <- restrict_to_allowed_cell_types(lines, allowed_cell_types)
        processed_predictions[[model]] <- lines
        cat("✅", model, "provided", length(lines), "predictions\n")
      } else {
        cat("⚠️", model, "returned", length(lines), "predictions, expected", length(processed_input), "\n")
      }
    } else {
      cat("❌", model, "failed to provide predictions\n")
    }
  }
  
  return(processed_predictions)
}

get_elite_subtype_predictions <- function(processed_input, tissuename, species, celltypename,
                                          elite_models, parallel, timeout_seconds,
                                          allowed_cell_types = NULL) {
  
  marker_data <- paste(names(processed_input), ':', processed_input, collapse = '\n')
  
  prompt <- glue::glue("Identify detailed cell subtypes for {celltypename} in {species} {tissuename} cells using the following markers, provided separately for each row.
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

  processed_predictions <- list()
  
  for (model in elite_models) {
    cat("🤖 Calling elite model for subtypes:", model, "\n")
    
    result <- call_llm_api(model, prompt, temperature = NULL, timeout_seconds = timeout_seconds)
    
    if (!is.null(result) && result != "") {
      lines <- trimws(unlist(strsplit(result, "\n")))
      lines <- lines[lines != ""]
      
      if (length(lines) == length(processed_input)) {
        lines <- restrict_to_allowed_cell_types(lines, allowed_cell_types)
        processed_predictions[[model]] <- lines
        cat("✅", model, "provided", length(lines), "subtype predictions\n")
      } else {
        cat("⚠️", model, "returned", length(lines), "predictions, expected", length(processed_input), "\n")
      }
    } else {
      cat("❌", model, "failed to provide subtype predictions\n")
    }
  }
  
  return(processed_predictions)
}

prepare_arbitration_data <- function(elite_predictions, processed_input) {
  
  cluster_names <- names(processed_input)
  prediction_lines <- c()
  
  for (i in 1:length(cluster_names)) {
    cluster_name <- cluster_names[i]

    model_predictions <- list()
    for (j in 1:length(elite_predictions)) {
      model_name <- paste0("LLM", j)
      model_results <- elite_predictions[[j]]
      
      if (i <= length(model_results)) {
        model_predictions[[model_name]] <- model_results[i]
      } else {
        model_predictions[[model_name]] <- "unknown"
      }
    }

    pred_parts <- paste0(names(model_predictions), ":", model_predictions)
    cluster_line <- paste0(i, ": ", paste(pred_parts, collapse = ", "))
    prediction_lines <- c(prediction_lines, cluster_line)
  }
  
  return(paste(prediction_lines, collapse = "\n"))
}

get_arbitrator_decision <- function(arbitration_data, tissuename, species, arbitrator_model,
                                    num_clusters, timeout_seconds, cluster_names,
                                    allowed_cell_types = NULL) {
  ontology_arbitration_instruction <- if (is.null(allowed_cell_types)) {
    paste(
      "Do not consolidate predictions that have a hierarchical relationship",
      "(e.g., parent-child classes in Cell Ontology) into a single entity.",
      paste0(
        "Before voting, filter out any prediction whose Cell Ontology ID is not a",
        " descendant of the declared lineage root for ", species, " ", tissuename, " cells."
      )
    )
  } else {
    ""
  }

  arbitration_prompt <- glue::glue("
Integrate multiple AI model predictions to determine the final cell type for {species} {tissuename} cells. Predictions are provided separately for each cell type.
{ontology_arbitration_instruction}
Determine the final annotation by majority vote. If this vote results in a tie, select the prediction from the first model listed.
{annotation_instruction}

{arbitration_data}",
                                 species = species,
                                 tissuename = tissuename,
                                 num_clusters = num_clusters,
                                 ontology_arbitration_instruction = ontology_arbitration_instruction,
                                 annotation_instruction = format_annotation_instruction(
                                   allowed_cell_types, num_clusters, subtype = FALSE
                                 ),
                                 arbitration_data = arbitration_data)
  
  cat("⚖️ Calling arbitrator model for final decision...\n")
  
  result <- call_llm_api(arbitrator_model, arbitration_prompt, temperature = NULL, timeout_seconds = timeout_seconds)
  
  if (is.null(result) || result == "") {
    return(NULL)
  }

  result_lines <- trimws(unlist(strsplit(result, "\n")))
  result_lines <- result_lines[result_lines != ""]
  
  if (length(result_lines) != num_clusters) {
    cat("⚠️ Arbitrator returned", length(result_lines), "results, expected", num_clusters, "\n")
    return(NULL)
  }

  final_results <- sub("^[^a-zA-Z]*([a-zA-Z].*)", "\\1", result_lines)
  names(final_results) <- cluster_names
  final_results <- restrict_to_allowed_cell_types(final_results, allowed_cell_types)
  
  return(final_results)
}

get_arbitrator_subtype_decision <- function(arbitration_data, tissuename, species, celltypename,
                                            arbitrator_model, num_clusters, timeout_seconds,
                                            cluster_names, allowed_cell_types = NULL) {
  ontology_arbitration_instruction <- if (is.null(allowed_cell_types)) {
    paste(
      "Do not consolidate predictions that have a hierarchical relationship",
      "(e.g., parent-child classes in Cell Ontology) into a single entity.",
      paste0(
        "Before voting, filter out any prediction whose Cell Ontology ID is not a",
        " descendant of the declared lineage root for ", celltypename, " in ", species,
        " ", tissuename, " cells."
      )
    )
  } else {
    ""
  }

  arbitration_prompt <- glue::glue("
Integrate multiple AI model predictions to determine the final cell subtype for {celltypename} in {species} {tissuename} cells. Predictions are provided separately for each cell subtype.
{ontology_arbitration_instruction}
Determine the final annotation by majority vote. If this vote results in a tie, select the prediction from the first model listed.
{annotation_instruction}

{arbitration_data}",
                                 celltypename = celltypename,
                                 species = species,
                                 tissuename = tissuename,
                                 num_clusters = num_clusters,
                                 ontology_arbitration_instruction = ontology_arbitration_instruction,
                                 annotation_instruction = format_annotation_instruction(
                                   allowed_cell_types, num_clusters, subtype = TRUE
                                 ),
                                 arbitration_data = arbitration_data)
  
  cat("⚖️ Calling arbitrator model for final subtype decision...\n")
  
  result <- call_llm_api(arbitrator_model, arbitration_prompt, temperature = NULL, timeout_seconds = timeout_seconds)
  
  if (is.null(result) || result == "") {
    return(NULL)
  }

  result_lines <- trimws(unlist(strsplit(result, "\n")))
  result_lines <- result_lines[result_lines != ""]
  
  if (length(result_lines) != num_clusters) {
    cat("⚠️ Arbitrator returned", length(result_lines), "results, expected", num_clusters, "\n")
    return(NULL)
  }

  final_results <- sub("^[^a-zA-Z]*([a-zA-Z].*)", "\\1", result_lines)
  names(final_results) <- cluster_names
  final_results <- restrict_to_allowed_cell_types(final_results, allowed_cell_types)
  
  return(final_results)
}
