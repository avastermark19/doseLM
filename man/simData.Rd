% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simData.R
\name{simData}
\alias{simData}
\title{simData Function to make simData.}
\usage{
simData(FA, MA, FZ, MZ)
}
\arguments{
\item{FA}{FA fold change parameter.}

\item{MA}{MA fold change parameter.}

\item{FZ}{FZ fold change parameter.}

\item{MZ}{MZ fold change parameter.}
}
\value{
Returns simData.
}
\description{
This function generates simData.
}
\details{
This function makes simData.
}
\examples{
se <- simData(FA=100)
str(se)
}
\references{
The 'doseLM' package, 2019 (in press).
}
\author{
AJ Vaestermark, JR Walters.
}
