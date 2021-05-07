
plco.request <- function(resource,
                         query = list(),
                         api_root = getOption(
                           "gwas.explorer.plco.api_root",
                           default = "https://exploregwas.cancer.gov/plco-atlas/api"
                         ),
                         api_key = getOption(
                           "gwas.explorer.plco.api_key",
                           default = NULL
                         ),
                         is_verbose = getOption(
                           "gwas.explorer.plco.is_verbose",
                           default = F
                         )) {
  query_defaults <- list(
    api_key = api_key
  )

  response <- httr::GET(
    url = paste(api_root, resource, sep = "/"),
    query = utils::modifyList(query, query_defaults),
    config = httr::config(verbose = is_verbose)
  )

  text <- httr::content(response, as = "text")
  results <- jsonlite::fromJSON(text, flatten = T)

  if (query$raw == "true" &&
    hasName(results, "data") &&
    hasName(results, "columns")) {
    colnames(results$data) <- results$columns
    results <- as.data.frame(results$data)
  }

  results
}



plco.get_metadata <- function(phenotype_id = NULL,
                              sex = NULL,
                              ancestry = NULL,
                              ...) {
  resource <- "metadata"
  params <- list(
    phenotype_id = phenotype_id,
    sex = sex,
    ancestry = ancestry
  )

  plco.request(resource, params, ...)
}

plco.get_pca <- function(phenotype_id, platform, pc_x, pc_y, ...) {
  resource <- "pca"
  params <- list(
    phenotype_id = phenotype_id,
    platform = platform,
    pc_x = pc_x,
    pc_y = pc_y,
    raw = "true"
  )

  plco.request(resource, params, ...)
}

plco.get_phenotypes <- function(query = NULL,
                                ...) {
  resource <- "phenotypes"
  params <- list(
    q = query
  )

  plco.request(resource, params, ...)
}


plco.get_phenotype <- function(phenotype_id,
                               type = "frequency",
                               ...) {
  resource <- "phenotype"
  params <- list(
    id = phenotype_id,
    type = type
  )

  plco.request(resource, params, ...)
}

plco.get_variants <- function(phenotype_id,
                              sex,
                              ancestry,
                              chromosome = NULL,
                              columns = NULL,
                              snp = NULL,
                              position = NULL,
                              position_min = NULL,
                              position_max = NULL,
                              p_value_min = NULL,
                              p_value_max = NULL,
                              order_by = "p_value",
                              order = "asc",
                              offset = 0,
                              limit = 1000,
                              ...) {
  resource <- "variants"
  params <- list(
    phenotype_id = phenotype_id,
    sex = sex,
    ancestry = ancestry,
    chromosome = chromosome,
    columns = columns,
    snp = snp,
    position = position,
    position_min = position_min,
    position_max = position_max,
    p_value_min = p_value_min,
    p_value_max = p_value_max,
    order_by = order_by,
    order = order,
    offset = offset,
    limit = limit
  )

  plco.request(resource, params, ...)
}

plco.download_dataset <- function(phenotype_id,
                                  filepath) {
  download_root <- plco.request("config", list(key = "downloadRoot"))
}

plco.get_summary <- function(phenotype_id,
                             sex,
                             ancestry) {
  resource <- "summary"
  params <- list(
    phenotype_id = phenotype_id,
    sex = sex,
    ancestry = ancestry
  )

  plco.request(resource, params, ...)
}

plco.get_points <- function(phenotype_id,
                            sex,
                            ancestry) {
  resource <- "points"
  params <- list(
    phenotype_id = phenotype_id,
    sex = sex,
    ancestry = ancestry
  )

  plco.request(resource, params, ...)
}
