#' Unified API Client for LLM Services
#' 
#' This module provides a unified interface for calling different LLM APIs
#' with standardized error handling and response processing.

#' Get effective temperature for a specific model
#' 
#' @param model Model name
#' @param requested_temperature The temperature requested by user
#' @return Effective temperature to use for this model
get_effective_temperature <- function(model, requested_temperature) {
  if (grepl("gpt-5", model, ignore.case = TRUE)) {
    return(1.0)
  }

  return(0.0)
}

#' Get model provider from model name
#' 
#' @param model Model name
#' @return Provider name
get_model_provider <- function(model) {
  config <- get_model_config(model)
  if (is.null(config)) {
    return(NULL)
  }
  return(config$provider)
}

# Resolve a key from the configured primary and fallback environment variables.
resolve_api_key <- function(config, api_key = NULL) {
  if (!is.null(api_key)) {
    return(api_key)
  }

  env_vars <- config$env_var
  if (!is.null(config$fallback_env_var)) {
    env_vars <- c(env_vars, config$fallback_env_var)
  }
  env_values <- Sys.getenv(env_vars, unset = "")
  api_key <- env_values[which(nzchar(env_values))[1]]
  if (is.na(api_key)) {
    api_key <- ""
  }
  api_key
}

# libcurl uses a zero timeout to mean no total request timeout. httr2 has the
# same behaviour when no timeout option is added to a request.
httr_timeout <- function(timeout_seconds) {
  if (is.null(timeout_seconds)) {
    return(httr::config(timeout = 0))
  }
  httr::timeout(timeout_seconds)
}

httr2_timeout <- function(req, timeout_seconds) {
  if (is.null(timeout_seconds)) {
    return(req)
  }
  httr2::req_timeout(req, timeout_seconds)
}

#' Make API call to LLM service
#' 
#' @param model Model name to use
#' @param prompt Text prompt to send
#' @param temperature Sampling temperature (default: 0.0)
#' @param timeout_seconds Request timeout in seconds (default: 1000)
#' @param api_key Optional API key (if not provided, will use environment variable)
#' @return API response content or NULL if failed
call_llm_api <- function(model, prompt, temperature = NULL, timeout_seconds = 1000, api_key = NULL) {

  config <- get_model_config(model)
  if (is.null(config)) {
    stop("❌ Model not supported: ", model)
  }
  
  if (is.null(temperature)) {
    temperature <- 0.0
  }

  effective_temperature <- get_effective_temperature(model, temperature)

  if (effective_temperature != temperature) {
    cat("🌡️ Model", model, "using temperature", effective_temperature, "instead of requested", temperature, "\n")
  }
  
  model_type <- get_model_provider(model)

  if (is.null(api_key)) {
    env_vars <- c(config$env_var, config$fallback_env_var)
    env_vars <- env_vars[!is.na(env_vars) & nzchar(env_vars)]
    api_key <- resolve_api_key(config)
    if (api_key == "") {
      cat("❌ API key not found for", config$provider, "(checked", paste(env_vars, collapse = ", "), ")\n")
      return(NULL)
    }
  }
  
  cat("🔌 Calling", config$provider, "API with model:", model, "\n")
  
  tryCatch({
    
    if (config$provider == "external_openai") {
      result <- call_external_openai_api(config, model, prompt, effective_temperature, timeout_seconds, api_key)

    } else if (config$provider == "openai") {
      result <- call_openai_api(config, model, prompt, effective_temperature, timeout_seconds, api_key)
      
    } else if (config$provider == "deepseek") {
      result <- call_deepseek_api(config, model, prompt, effective_temperature, timeout_seconds, api_key)
      
    } else if (config$provider == "claude") {
      result <- call_claude_api(config, model, prompt, effective_temperature, timeout_seconds, api_key)
      
    } else if (config$provider == "gemini") {
      result <- call_gemini_api(config, model, prompt, effective_temperature, timeout_seconds, api_key)
      
    } else if (config$provider == "grok") {
      result <- call_grok_api(config, model, prompt, effective_temperature, timeout_seconds, api_key)
      
    } else if (config$provider == "kimi") {
      result <- call_kimi_api(config, model, prompt, effective_temperature, timeout_seconds, api_key)
      
    } else if (config$provider == "doubao") {
      result <- call_doubao_api(config, model, prompt, effective_temperature, timeout_seconds, api_key)
      
    } else {
      stop("❌ Provider not implemented: ", config$provider)
    }
    
    cat("✅ API call successful\n")
    return(result)
    
  }, error = function(e) {
    cat("❌ API call failed for", config$provider, ":", e$message, "\n")
    return(NULL)
  })
}

call_openai_api <- function(config, model, prompt, temperature, timeout_seconds, api_key) {
  if (!is.null(timeout_seconds)) {
    timeout_seconds <- if (!is.null(config$api$timeout)) config$api$timeout else 360
  }
  
  response <- httr::POST(
    "https://api.openai.com/v1/chat/completions",
    httr::add_headers(
      "Content-Type" = "application/json",
      "Authorization" = paste("Bearer", api_key)
    ),
    body = jsonlite::toJSON(list(
      model = model,
      messages = list(
        list(role = "user", content = prompt)
      ),
      temperature = temperature
    ), auto_unbox = TRUE),
    httr_timeout(timeout_seconds),
    httr::config(connecttimeout = if (is.null(timeout_seconds)) 0 else 360),
    encode = "raw"
  )
  
  if (httr::status_code(response) == 200) {
    resp_content <- httr::content(response, "text", encoding = "UTF-8")
    resp_json <- jsonlite::fromJSON(resp_content, flatten = TRUE)
    return(resp_json$choices$message.content[1])
  } else {
    stop("HTTP error: ", httr::status_code(response))
  }
}

call_openai_compatible_api <- function(config, model, prompt, temperature, timeout_seconds, api_key) {
  request_body <- list(
    model = model,
    messages = list(list(role = "user", content = prompt))
  )
  if (!is.null(temperature)) {
    request_body$temperature <- temperature
  }

  response <- httr::POST(
    config$endpoint,
    httr::add_headers(
      "Content-Type" = "application/json",
      "Authorization" = paste(config$auth_prefix, api_key)
    ),
    body = jsonlite::toJSON(request_body, auto_unbox = TRUE),
    httr_timeout(timeout_seconds),
    encode = "raw"
  )

  if (httr::status_code(response) != 200) {
    stop("HTTP error: ", httr::status_code(response))
  }

  resp <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"), flatten = TRUE)
  return(resp$choices$message.content[1])
}

# Call an OpenAI Responses-compatible endpoint (for example, a relay service).
# Responses reasoning models do not use chat-completions' temperature field.
call_external_openai_api <- function(config, model, prompt, temperature, timeout_seconds, api_key) {
  base_url <- Sys.getenv(config$base_url_env_var, unset = config$base_url)
  base_url <- sub("/+$", "", base_url)
  endpoint_path <- Sys.getenv(
    config$endpoint_path_env_var,
    unset = config$endpoint_path
  )
  endpoint <- if (grepl("/responses$", base_url, ignore.case = TRUE)) {
    base_url
  } else if (grepl("/v1$", base_url, ignore.case = TRUE)) {
    paste0(base_url, "/responses")
  } else {
    paste0(base_url, endpoint_path)
  }

  reasoning_effort <- Sys.getenv(
    "DEEPCELLSEEK_REASONING_EFFORT",
    unset = if (!is.null(config$reasoning_effort)) config$reasoning_effort else "max"
  )
  request_body <- list(
    model = model,
    input = prompt,
    reasoning = list(effort = reasoning_effort),
    # Streaming keeps the relay connection active while max reasoning runs,
    # avoiding the relay's 120-second idle proxy timeout.
    stream = TRUE
  )

  response <- httr::POST(
    endpoint,
    httr::add_headers(
      "Content-Type" = "application/json",
      "Authorization" = paste(config$auth_prefix, api_key)
    ),
    body = jsonlite::toJSON(request_body, auto_unbox = TRUE, null = "null"),
    httr_timeout(timeout_seconds),
    httr::config(connecttimeout = if (is.null(timeout_seconds)) 0 else min(timeout_seconds, 360)),
    encode = "raw"
  )

  status <- httr::status_code(response)
  if (status < 200 || status >= 300) {
    error_body <- httr::content(response, "text", encoding = "UTF-8")
    stop("HTTP error: ", status, if (nzchar(error_body)) paste0(" - ", error_body))
  }

  response_body <- httr::content(response, "text", encoding = "UTF-8")
  extract_stream_response_text(response_body)
}

# Parse an SSE Responses response and extract the completed response payload.
extract_stream_response_text <- function(response_body) {
  lines <- strsplit(response_body, "\\r?\\n", perl = TRUE)[[1]]
  data_lines <- sub("^data:\\s*", "", lines[grepl("^data:\\s*", lines)])
  data_lines <- data_lines[nzchar(data_lines) & data_lines != "[DONE]"]

  if (!length(data_lines)) {
    # Some relays ignore stream=TRUE and return one ordinary JSON response.
    resp <- jsonlite::fromJSON(response_body, simplifyVector = FALSE)
    return(extract_responses_text(resp))
  }

  for (data in rev(data_lines)) {
    event <- tryCatch(
      jsonlite::fromJSON(data, simplifyVector = FALSE),
      error = function(e) NULL
    )
    if (!is.null(event$response)) {
      return(extract_responses_text(event$response))
    }
  }

  stop("Responses stream completed without output text")
}

# Responses APIs may expose output text directly or nest it in output/content.
extract_responses_text <- function(resp) {
  if (!is.null(resp$output_text) && length(resp$output_text) > 0) {
    return(paste(unlist(resp$output_text), collapse = ""))
  }

  texts <- character()
  output <- if (!is.null(resp$output)) resp$output else list()
  for (item in output) {
    content <- if (!is.null(item$content)) item$content else list()
    for (part in content) {
      if (!is.null(part$text)) {
        texts <- c(texts, as.character(part$text))
      }
    }
  }

  if (!length(texts)) {
    stop("Responses API returned no output text")
  }
  paste(texts, collapse = "")
}

call_deepseek_api <- function(config, model, prompt, temperature, timeout_seconds, api_key) {
  call_openai_compatible_api(config, model, prompt, temperature, timeout_seconds, api_key)
}

call_claude_api <- function(config, model, prompt, temperature, timeout_seconds, api_key) {
  req <- httr2::request(config$endpoint) |>
    httr2::req_headers(
      "Content-Type" = "application/json",
      !!config$auth_header := api_key,
      "anthropic-version" = "2023-06-01"
    ) |>
    httr2::req_body_json(list(
      model = model,
      messages = list(list(role = "user", content = prompt)),
      temperature = temperature,
      max_tokens = config$max_tokens %||% 4096
    )) |>
    httr2_timeout(timeout_seconds) |>
    httr2::req_perform()
  
  resp <- httr2::resp_body_json(req)
  return(resp$content[[1]]$text)
}

call_gemini_api <- function(config, model, prompt, temperature, timeout_seconds, api_key) {
  endpoint <- paste0(config$endpoint, "/", model, ":generateContent?key=", api_key)
  
  req <- httr2::request(endpoint) |>
    httr2::req_headers("Content-Type" = "application/json") |>
    httr2::req_body_json(list(
      contents = list(list(parts = list(list(text = prompt)))),
      generationConfig = list(temperature = temperature)
    )) |>
    httr2_timeout(timeout_seconds) |>
    httr2::req_perform()
  
  resp <- httr2::resp_body_json(req)
  return(resp$candidates[[1]]$content$parts[[1]]$text)
}

call_grok_api <- function(config, model, prompt, temperature, timeout_seconds, api_key) {
  req <- httr2::request(config$endpoint) |>
    httr2::req_headers(
      "Content-Type" = "application/json",
      !!config$auth_header := paste(config$auth_prefix, api_key)
    ) |>
    httr2::req_body_json(list(
      model = model,
      messages = list(list(role = "user", content = prompt)),
      temperature = temperature
    )) |>
    httr2_timeout(timeout_seconds) |>
    httr2::req_perform()
  
  resp <- httr2::resp_body_json(req)
  return(resp$choices[[1]]$message$content)
}

call_kimi_api <- function(config, model, prompt, temperature, timeout_seconds, api_key) {
  if (model %in% c("kimi-k2.6", "kimi-k2.5")) {
    temperature <- NULL
  }
  call_openai_compatible_api(config, model, prompt, temperature, timeout_seconds, api_key)
}

call_doubao_api <- function(config, model, prompt, temperature, timeout_seconds, api_key) {
  req <- httr2::request(config$endpoint) |>
    httr2::req_headers(
      "Content-Type" = "application/json",
      !!config$auth_header := paste(config$auth_prefix, api_key)
    ) |>
    httr2::req_body_json(list(
      model = model,
      messages = list(list(role = "user", content = prompt)),
      temperature = temperature
    )) |>
    httr2_timeout(timeout_seconds) |>
    httr2::req_perform()
  
  resp <- httr2::resp_body_json(req)
  return(resp$choices[[1]]$message$content)
}

#' Call multiple models in parallel
#' 
#' @param models Vector of model names
#' @param prompt Text prompt to send to all models
#' @param temperature Sampling temperature (will be adjusted per model)
#' @param timeout_seconds Request timeout
#' @return Named list of results (model_name -> response)
call_models_parallel <- function(models, prompt, temperature = NULL, timeout_seconds = 1000) {
  
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' is required for parallel processing. Please install it.")
  }
  
  if (is.null(temperature)) {
    temperature <- 0.0
  }

  future::plan(future::multisession, workers = min(length(models), 4))

  futures <- list()
  for (model in models) {
    futures[[model]] <- future::future({
      call_llm_api(model, prompt, temperature, timeout_seconds)
    }, seed = TRUE)
  }

  results <- list()
  for (model in names(futures)) {
    tryCatch({
      results[[model]] <- future::value(futures[[model]])
    }, error = function(e) {
      cat("❌ Failed to get result from", model, ":", e$message, "\n")
      results[[model]] <- NULL
    })
  }

  future::plan(future::sequential)

  successful_results <- results[!sapply(results, is.null)]
  
  return(successful_results)
}
