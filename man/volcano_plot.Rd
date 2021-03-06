% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volcano_plot.R
\name{volcano_plot}
\alias{volcano_plot}
\title{Function to draw a Volcano Plot}
\usage{
volcano_plot(
  pval,
  fc,
  names,
  size = 2,
  tFC = 2,
  tPV = -log10(0.001),
  show.labels = TRUE,
  show.effect = FALSE
)
}
\arguments{
\item{pval}{numeric vector of P.Values}

\item{fc}{numeric vector of fold change}

\item{names}{character vector with the feature's names.}

\item{size}{(default \code{2}) Sice of the labels in case they are
placed.}

\item{tFC}{(default \code{2}) fold change threshold. It can be set to
\code{NULL} to not filter.}

\item{tPV}{(default \code{-log10(0.001)}) P-Value threshold. It can be set
to \code{NULL} to not filter.}

\item{show.labels}{(default \code{TRUE}) If set to \code{TRUE}, features are 
labelled.}

\item{show.effect}{(default \code{FALSE}) If set to \code{TRUE}, the X-axis
will should \code{2^logFC} instead to the default \code{logFC}.}
}
\value{
A \code{ggplot} object
}
\description{
Function that takes two numeric vectors (P-Value and fold change)
and draws a volcano plot using \code{ggplot2}
}
\examples{
data(rset)
w1 <- getAssociation(rset, rid = 1, fNames = NULL)
volcano_plot(w1$P.Value, w1$logFC, rownames(w1))
}
