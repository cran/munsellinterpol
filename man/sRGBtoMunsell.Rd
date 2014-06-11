\name{sRGBtoMunsell}
\alias{sRGBtoMunsell}
\title{Convert sRGB coordinates into a Munsell specification}
\description{\code{sRGBtoMunsell} Converts sRGB coordinates into a Munsell 
specification, by interpolating over the extrapolated Munsell renotation data
}
\usage{sRGBtoMunsell(sRGB) }
\arguments{
  \item{sRGB}{ sRGB coordinates}
}
\value{
  \item{MunsellSpec}{ standard Munsell specification or a Munsell vector in 
  ColorLab format }
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
#sRGBtoMunsell(c(23,45,67))
}
\keyword{datasets}

