% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plco.R
\name{plco.get_metadata}
\alias{plco.get_metadata}
\title{Retrieves variant metadata for the specified phenotypes}
\usage{
plco.get_metadata(phenotype_id = NULL, sex = NULL, ancestry = NULL, ...)
}
\arguments{
\item{phenotype_id}{Optional. A numeric phenotype id}

\item{sex}{Optional. A character vector specifying sexes to retrieve data for}

\item{ancestry}{Optional. A character vector specifying ancestries to retrieve data for}
}
\value{
A dataframe containing phenotype metadata
}
\description{
Retrieves variant metadata for the specified phenotypes
}
\examples{
plco.get_metadata(1010, "female", "european")
}
