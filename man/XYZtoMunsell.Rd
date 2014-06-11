\name{XYZtoMunsell}
\alias{XYZtoMunsell}
\title{Convert XYZ coordinates into a Munsell specification}
\description{\code{XYZtoMunsell} Converts XYZ coordinates into a Munsell 
specification, by interpolating over the extrapolated Munsell renotation data
}
\usage{XYZtoMunsell(XYZ) }
\arguments{
  \item{XYZ}{ CIE XYZ coordinates of MunsellSpec illuminated by Illuminant C}
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
XYZtoMunsell(c(0.310897, 0.306510, 74.613450))
}
\keyword{datasets}

