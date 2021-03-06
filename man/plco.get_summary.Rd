% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plco.R
\name{plco.get_summary}
\alias{plco.get_summary}
\title{Retrieves aggregate variants for the specified phenotype, sex, and ancestry}
\usage{
plco.get_summary(phenotype_id, sex, ancestry, p_value_nlog_min = 2, ...)
}
\arguments{
\item{phenotype_id}{A numeric phenotype id}

\item{sex}{A character vector specifying a sex to retrieve data for}

\item{ancestry}{A character vector specifying ancestries to retrieve data for}

\item{p_value_nlog_min}{Optional. A numeric value specifying the minimum p-value for aggregate variants}
}
\value{
A dataframe containing aggregated variants
}
\description{
Retrieves aggregate variants for the specified phenotype, sex, and ancestry
}
\examples{
plco.get_summary(1010, "female", "european")
}
