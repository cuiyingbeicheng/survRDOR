\name{sim_data}
\alias{sim_data}
\docType{data}
\title{
Data example
}
\description{
This is a simulated data example.
}
\usage{data("sim_data")}
\format{
  A data frame with 200 observations on the following 4 variables.
  \describe{
  \item{\code{y1}}{Time to response (month)}
  \item{\code{y2}}{Progression-free survival time (month)}
  \item{\code{d1}}{Status for response (=1 for responders; =0 for non-responders)}
  \item{\code{d2}}{Status for progression-free survival (=1 for event; =0 for censored)}
  }
}
\examples{
data(sim_data)
str(sim_data)
}
\keyword{datasets}
