\name{MunsellHVC}
\alias{MunsellHVC}
\title{Munsell string to HVC}
\description{\code{MunsellHVC} Converts a Munsell string to HVC}
\usage{MunsellHVC(MunIn) }
\arguments{
  \item{MunIn}{ A Munsell specification in Munsell format }
}
\value{
    \item{H}{ Munsell hue}
    \item{V}{ Munsell value}
    \item{C}{ Munsell chroma}
  \item{Status.ind}{ status of the operation}
}
\author{Jose Gama}
\examples{
MunsellHVC('7.6P 8.9/2.2')
}
\keyword{datasets}

