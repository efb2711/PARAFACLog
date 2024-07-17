\name{X3way}
\alias{X3way}
\docType{data}
\title{Three way data array}
\description{The data contains the reactions of a group of people to certain situations}
\usage{data(X3way)}
\format{
A data array with a hipotetical 8 persons x 7 emotions x 6 situations:
\subsection{$ The situations are the following:}{
  \describe{
    \item{\code{S1}}{Quarrelling with someone}
    \item{\code{S2}}{Partner leaves you}
    \item{\code{S3}}{Someone is telling lies about you}
    \item{\code{S4}}{Giving a bad speech}
    \item{\code{S5}}{Failing a test}
    \item{\code{S6}}{Writing a bad paper}
  }}
  \subsection{$ The responses are the following:}{
  \describe{
    \item{\code{R1}}{Other anger}
    \item{\code{R2}}{Shame}
    \item{\code{R3}}{Love}
    \item{\code{R4}}{Sorrow}
    \item{\code{R5}}{Fear}
    \item{\code{R6}}{Guilt}
    \item{\code{R7}}{Self anger}
  }}
}
\examples{
data(X3way)
str(X3way)
}
\keyword{datasets}
