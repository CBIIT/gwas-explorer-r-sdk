#' gwas.explorer: A package for qyerying GWAS results from studies supported by
#' the National Cancer Institute
#'
#' The GWAS Explorer SDK enables users to query GWAS results and
#' metadata from studies that were supported by the National Cancer Institute.
#' The SDK currently provides GWAS results for the Prostate, Lung, Colorectal,
#' and Ovarian Cancer Screening Trial (PLCO).
#'
#'
#' @docType package
#' @name gwas.explorer
NULL


#' Performs a request against the PLCO API
#'
#' @param resource A resource name
#' @param query A list of query parameters
#' @param api_root The API root of the request
#' @param api_key An API key to apply to the request
#' @param is_verbose If TRUE, enables verbose mode
#' @return The results of the qpi request.
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

  # return non-json responses as-is
  if (!grepl("^application/json", response$header$`content-type`, ignore.case = T)) {
    return(text)
  }

  # otherwise, parse and flatten json responses
  results <- jsonlite::fromJSON(text, flatten = T)

  if (hasName(results, "error")) {
    stop(paste(results$error, results$message, sep = ": "))
  }

  if (hasName(query, "raw") &&
    query$raw == "true" &&
    hasName(results, "data") &&
    hasName(results, "columns")) {
    data <- results$data
    colnames(data) <- results$columns
    return (data)
  }

  results
}

#' Retrieves variant metadata for the specified phenotypes
#'
#' @param phenotype_id Optional. A numeric phenotype id
#' @param sex Optional. A character vector specifying sexes to retrieve data for
#' @param ancestry Optional. A character vector specifying ancestries to retrieve data for
#' @return A dataframe containing phenotype metadata
#' @examples
#' plco.get_metadata(1010, "female", "european")
plco.get_metadata <- function(phenotype_id = NULL,
                              sex = NULL,
                              ancestry = NULL,
                              ...) {
  resource <- "metadata"
  params <- list(
    phenotype_id = phenotype_id,
    sex = if (is.null(sex) || length(sex) == 0) NULL
          else paste0(sex, collapse = ','),
    ancestry = if (is.null(ancestry) || length(ancestry) == 0) NULL
          else paste0(ancestry, collapse = ',')
  )

  plco.request(resource, params, ...)
}

#' Retrieve PCA coordinates for the specified phenotype and platform
#'
#' @param phenotype_id A numeric phenotype id
#' @param platform A character vector specifying the platform to retrieve data for
#' @param pc_x A numeric value (1-20) specifying the x axis's principal component
#' @param pc_y A numeric value (1-20) specifying the y axis's principal component
#' @return A dataframe containing pca coordinates
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

#' Retrieve phenotypes
#'
#' @param query A character vector specifying a term to search for
#' @return If query is specified, returns a list of phenotypes which contain the query term. Otherwise, returns all phenotypes. Note that phenotypes without an import_count do not contain associated gwas summary results
#' @examples
#' plco.get_phenotypes("cancer")
plco.get_phenotypes <- function(query = NULL,
                                ...) {
  resource <- "phenotypes"
  params <- list(
    q = query
  )

  plco.request(resource, params, ...)
}

#' Retrieves aggregate counts for participants. Aggregate counts under 10 are returned as "< 10".
#'
#' @param columns A character vector specifying properties for which to retrieve counts for. Valid properties are: value, ancestry, genetic_ancestry, sex, and age
#' @param precision Optional. For continuous phenotypes, a numeric value specifying the -log10(precision) to which values should be rounded to.
#' @examples
#' plco.get_participants(1010)
#' plco.get_participants(1010, c("value", "age"))
#' plco.get_participants(2250, c("value"), -1)
plco.get_participants <- function(phenotype_id,
                                  columns=c("value"),
                                  precision=0,
                                  ...) {
  resource <- "participants"
  params <- list(
    phenotype_id = phenotype_id,
    columns = paste0(columns, collapse = ','),
    precision = precision,
    raw = "true"
  )

  plco.request(resource, params, ...)
}

#' Retrieves variants for a specified phenotype, sex, and ancestry
#'
#' @param phenotype_id A numeric phenotype id
#' @param sex A character vector specifying a sex to retrieve data for
#' @param ancestry A character vector specifying ancestries to retrieve data for
#' @param chromosome Optional. A character vector specifying a chromosome to retrieve data for
#' @param columns Optional. A character vector specifying properties to include for each variant
#' @param snp Optional. A character vector specifying snps to retrieve
#' @param position Optional. A numeric value specifying the chromosome position to retrieve data for
#' @param position_min Optional. A numeric value specifying the minimum chromosome position for variants
#' @param position_max Optional. A numeric value specifying the maximum chromosome position  for variants
#' @param p_value_min Optional. A numeric value specifying the minimum p-value for variants
#' @param p_value_max Optional. A numeric value specifying the maximum p-value for variants
#' @param order_by Optional. A character vector specifying the property to sort variants by
#' @param order Optional. A character vector specifying the order to sort variants by. Either "asc" or "desc"
#' @param offset Optional. A numeric value specifying the offset from which to begin retrieving variants
#' @param limit Optional. A numeric value specifying the maximum number of variants to retrieve
#' @return A dataframe containing variants
#' @examples
#' plco.get_variants(1010, "female", "european", limit=100)
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
    limit = limit,
    raw = "true"
  )

  plco.request(resource, params, ...)
}

#' Downloads the original association results in tsv.gz format to the specified path
#'
#' @param phenotype_id A numeric phenotype id
#' @param sex A character vector specifying the path to save the file to
#' @examples
#' plco.download_dataset(1010, "download.tsv.gz")
plco.download_dataset <- function(phenotype_id,
                                  filepath,
                                  ...) {
  resource <- "download"
  params <- list(
    phenotype_id = phenotype_id,
    get_link_only = "true"
  )

  url <- plco.request(resource, params, ...)

  httr::GET(url, httr::write_disk(filepath, overwrite=TRUE))
}

#' Retrieves aggregate variants for the specified phenotype, sex, and ancestry
#'
#' @param phenotype_id A numeric phenotype id
#' @param sex A character vector specifying a sex to retrieve data for
#' @param ancestry A character vector specifying ancestries to retrieve data for
#' @param p_value_nlog_min Optional. A numeric value specifying the minimum p-value for aggregate variants
#' @return A dataframe containing aggregated variants
#' @examples
#' plco.get_summary(1010, "female", "european")
plco.get_summary <- function(phenotype_id,
                             sex,
                             ancestry,
                             p_value_nlog_min = 2,
                             ...) {
  resource <- "summary"
  params <- list(
    phenotype_id = phenotype_id,
    sex = sex,
    ancestry = ancestry,
    p_value_nlog_min = p_value_nlog_min,
    raw = "true"
  )

  plco.request(resource, params, ...)
}

#' Retrieves sampled variants suitable for visualizing a QQ plot for the specified phenotype, sex, and ancestry
#'
#' @param phenotype_id A numeric phenotype id
#' @param sex A character vector specifying a sex to retrieve data for
#' @param ancestry A character vector specifying ancestries to retrieve data for
#' @return A dataframe containing variants
#' @examples
#' plco.get_points(1010, "female", "european")
plco.get_points <- function(phenotype_id,
                            sex,
                            ancestry,
                            ...) {
  resource <- "points"
  params <- list(
    phenotype_id = phenotype_id,
    sex = sex,
    ancestry = ancestry,
    raw = "true"
  )

  plco.request(resource, params, ...)
}
