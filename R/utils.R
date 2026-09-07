#' Data Processing Utilities for DeepCellSeek
#' 
#' This module provides common utility functions used across different
#' components of the DeepCellSeek package.

#' Process input data for cell annotation
#' 
#' This function standardizes different input formats (data frames, lists, vectors)
#' into a consistent format for downstream processing by both local API calls 
#' and DeepCellSeek gateway API calls.
#' 
#' @importFrom magrittr %>%
#' @param input Input data in various formats:
#'   - data.frame: Seurat FindAllMarkers output with columns: gene, cluster, avg_log2FC, p_val_adj
#'   - list: Named list where each element contains gene names for a cluster
#'   - character: Already processed gene list
#' @param topgenenumber Number of top genes to select per cluster (default: 10)
#' @return Named vector where names are cluster identifiers and values are 
#'   comma-separated gene lists
#' @export
process_input_data <- function(input, topgenenumber = 10) {
  if (is.character(input) || (is.list(input) && !is.data.frame(input))) {
    if (is.character(input)) {
      processed <- input
    } else {
      processed <- sapply(input, paste, collapse = ',')
    }
  } 
  else if (is.data.frame(input)) {
    required_cols <- c("gene", "cluster", "avg_log2FC", "p_val_adj")
    missing_cols <- setdiff(required_cols, names(input))
    if (length(missing_cols) > 0) {
      stop("❌ Input data.frame is missing required columns: ", 
           paste(missing_cols, collapse = ", "))
    }

    if (!requireNamespace("dplyr", quietly = TRUE)) {
      input <- input[input$avg_log2FC > 0, , drop = FALSE]
      processed <- tapply(input$gene, list(input$cluster), 
                         function(i) paste0(i[1:min(length(i), topgenenumber)], collapse = ','))
    } else {
      top_genes <- input %>%
        dplyr::filter(p_val_adj < 2.2e-16, avg_log2FC > 0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        dplyr::slice_head(n = topgenenumber) %>%
        dplyr::ungroup()
      
      processed <- tapply(top_genes$gene, list(top_genes$cluster), 
                         function(i) paste0(i, collapse = ','))
    }
  } 
  else {
    stop("❌ Input must be either a data frame with gene, cluster, avg_log2FC, p_val_adj columns or a processed gene list")
  }
  
  return(processed)
}

# Resolve the optional set of labels an annotation is permitted to use.  Keeping
# this separate from prompt construction lets every annotation entry point use
# the same validation behaviour.
resolve_allowed_cell_types <- function(allowed_cell_types) {
  if (is.null(allowed_cell_types)) {
    return(NULL)
  }

  source_description <- "the supplied value"
  if (is.character(allowed_cell_types) && length(allowed_cell_types) == 1L &&
      file.exists(allowed_cell_types)) {
    source_description <- paste0("RDS file '", allowed_cell_types, "'")
    allowed_cell_types <- readRDS(allowed_cell_types)
  }

  if (is.data.frame(allowed_cell_types)) {
    if (!"cell_type" %in% names(allowed_cell_types)) {
      stop(
        "❌ Allowed cell types supplied as a data frame must contain a 'cell_type' column"
      )
    }
    allowed_cell_types <- allowed_cell_types$cell_type
  }

  if (is.factor(allowed_cell_types)) {
    allowed_cell_types <- as.character(allowed_cell_types)
  }
  if (!is.character(allowed_cell_types)) {
    stop(
      "❌ allowed_cell_types must be a character vector, an RDS path, or a data frame with a 'cell_type' column"
    )
  }

  allowed_cell_types <- trimws(allowed_cell_types)
  allowed_cell_types <- unique(allowed_cell_types[!is.na(allowed_cell_types) & nzchar(allowed_cell_types)])
  if (!length(allowed_cell_types)) {
    stop("❌ No usable allowed cell types were found in ", source_description)
  }

  allowed_cell_types
}

format_annotation_instruction <- function(allowed_cell_types, num_clusters,
                                           subtype = FALSE) {
  if (is.null(allowed_cell_types)) {
    label_instruction <- if (subtype) {
      paste0(
        "You MUST use standardized cell subtype names from the Cell Ontology (CL). ",
        "If no exact CL term exists, use the most specific CL term available and ",
        "add descriptive modifiers."
      )
    } else {
      paste0(
        "You MUST use standardized cell type names from the Cell Ontology (CL). ",
        "If no exact CL term exists, use the most specific CL term available and ",
        "add descriptive modifiers."
      )
    }
    mixture_instruction <- if (subtype) {
      "\nNote: Some rows may represent a mixture of multiple subtypes."
    } else {
      "\nSome can be a mixture of multiple cell types."
    }
    return(paste0(
      label_instruction,
      "\nOnly ", if (subtype) "output" else "provide",
      " the cell ", if (subtype) "subtype" else "type",
      " name. Do not include any numbers or extra annotations before the name. ",
      "Do not include any explanatory text, introductory phrases, or descriptions.",
      "\nIMPORTANT: Return exactly ", num_clusters,
      " lines, one for each row.",
      mixture_instruction
    ))
  }

  paste0(
    "You MUST choose every annotation only from the permitted cell-type labels below. ",
    "Return one or more selected labels exactly as written. If a cluster contains ",
    "multiple cell types, join the labels with ' + '. Every component must be one ",
    "of the listed labels; do not use synonyms, modifiers, or any label outside ",
    "this list.\n",
    "This label list is authoritative and takes precedence over all other label rules.\n",
    paste0("- ", allowed_cell_types, collapse = "\n"),
    "\nIMPORTANT: Return exactly ", num_clusters,
    " lines, one for each row.\n",
    "Only output the selected label(s), with no numbers, explanations, or extra text."
  )
}

# Canonicalise model output back to the spelling supplied by the caller.  An
# unrecognised output becomes a status value rather than an out-of-vocabulary
# cell type, which keeps the constraint strict even if a model ignores it.
restrict_to_allowed_cell_types <- function(results, allowed_cell_types) {
  if (is.null(allowed_cell_types)) {
    return(results)
  }

  normalise <- function(x) {
    tolower(gsub("\\s+", " ", trimws(x)))
  }

  allowed_normalised <- normalise(allowed_cell_types)
  normalise_composite <- function(x) {
    # Delimiters require whitespace around '+', so marker names such as
    # "CD4+ T cell" remain a single candidate label.
    parts <- trimws(unlist(strsplit(x, "\\s+\\+\\s+")))
    if (!length(parts) || any(!nzchar(parts))) {
      return(NA_character_)
    }
    parts <- normalise(parts)
    if (anyDuplicated(parts)) {
      return(NA_character_)
    }
    matched_parts <- match(parts, allowed_normalised)
    if (anyNA(matched_parts)) {
      return(NA_character_)
    }
    paste(matched_parts, collapse = "+")
  }

  canonical_keys <- vapply(results, normalise_composite, character(1))
  invalid <- is.na(canonical_keys)

  if (any(invalid)) {
    warning(
      sum(invalid),
      " annotation(s) were outside allowed_cell_types and were marked as failed.",
      call. = FALSE
    )
  }

  constrained <- results
  if (any(!invalid)) {
    constrained[!invalid] <- vapply(strsplit(canonical_keys[!invalid], "+", fixed = TRUE),
                                    function(indices) paste(allowed_cell_types[as.integer(indices)], collapse = " + "),
                                    character(1))
  }
  constrained[invalid] <- "Failed to annotate: outside allowed cell types"
  constrained
}
