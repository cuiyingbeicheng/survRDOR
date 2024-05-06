\name{get.SRDOR_IPW.naive}
\alias{get.SRDOR_IPW.naive}
\title{
Naive IPW adjusted estimator for the distribution of restricted duration of response.
}
\description{
Function to estimate the survival function of the restricted duration of response using the naive IPW adjusted estimator.
}
\usage{
get.SRDOR_IPW.naive(X1,X2,delta1,delta2,tau)
}
\arguments{
  \item{X1}{Time to response (month)}
  \item{X2}{Progression-free survival time (month)}
  \item{delta1}{Status for response (=1 for responders; =0 for non-responders)}
  \item{delta2}{Status for progression-free survival (=1 for event; =0 for censored)}
  \item{tau}{Truncation time for defining the time-window of interest.}
}
\value{
\item{td}{Time points of interest for the survival function of the duration of response.}
\item{est.curve}{The distribution of DOR during the time window [0, tau]}
\item{est.mu}{The mean response-duration during the time window [0, tau]}
\item{se.curve}{The standard error of est.curve}
\item{se.mu}{The standard error of est.mu}
}

\examples{
\donttest{
res <- get.SRDOR_IPW.naive(X1=sim_data$y1,X2=sim_data$y2,delta1=sim_data$d1,delta2=sim_data$d2,
tau = 30)
}
}

\keyword{inference functions}
