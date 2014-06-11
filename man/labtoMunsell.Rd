\name{labtoMunsell}
\alias{labtoMunsell}
\title{Convert CIE Lab coordinates into a Munsell specification}
\description{\code{labtoMunsell} Converts CIE Lab coordinates into a Munsell 
specification, by interpolating over the extrapolated Munsell renotation data
}
\usage{labtoMunsell(Lab) }
\arguments{
  \item{Lab}{ CIE Lab coordinates of MunsellSpec illuminated by Illuminant C}
}
\value{
  \item{MunsellSpec}{ A Munsell specification in Munsell format }
  \item{MunsellVec}{ A Munsell specification, in ColorLab format}
  \item{Status.ind}{ status of the operation}
}
\source{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{http://www.99main.com/~centore/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}
\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{http://www.99main.com/~centore/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}
\author{Jose Gama}
\examples{
labtoMunsell(c(0.310897, 0.306510, 74.613450))
}
\keyword{datasets}

