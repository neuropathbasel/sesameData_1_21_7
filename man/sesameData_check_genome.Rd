% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genome.R
\name{sesameData_check_genome}
\alias{sesameData_check_genome}
\title{Find genome assembly version(s) supported for a platform}
\usage{
sesameData_check_genome(genome, platform)
}
\arguments{
\item{genome}{mm10, hg38, ..., or NULL}

\item{platform}{HM27, HM450, EPIC, EPICv2, MSA, ...}
}
\value{
genome as string
}
\description{
Find genome assembly version(s) supported for a platform
}
\examples{
sesameData_check_genome(NULL, "Mammal40")
}
