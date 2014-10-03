\name{getIdeo}
\alias{getIdeo}
\title{Download ideogram table from UCSC}
\description{Download table containing chromosomal extent and band locations from the UCSC genome browser}
\usage{getIdeo(ideoSource)}
\arguments{
  \item{ideoSource}{(character) Genome build for data (e.g. mm10).}
}
\details{Uses \code{rtracklayer} to retrieve the \emph{cytoBandIdeo} table from the UCSC genome browser. The \emph{cytoBandIdeo} table  contains chromosomal ideogram information and is used to graph the chromosomal bands in \code{plotOnIdeo()}. This table is provided as input to \code{plotOnIdeo()}. In the case where the user bins the data, the output of this function can also be used as input to generate bin coordinates for binning the data (see \code{avgByBin()}). }
\value{(data.frame) ideogram table}
\examples{getIdeo("mm9")}
\seealso{\code{avgByBin()},\code{getBins()}}
\author{Shraddha Pai <Shraddha.Pai@camh.ca>, Jingliang Ren}





