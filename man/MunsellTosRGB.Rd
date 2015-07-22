\name{MunsellTosRGB}
\alias{MunsellTosRGB}
\title{Convert a Munsell specification into sRGB coordinates}
\description{\code{MunsellTosRGB} Converts a 
Munsell specification into sRGB coordinates, by interpolating over the 
extrapolated Munsell renotation data
}
\usage{MunsellTosRGB(MunsellSpec) }
\arguments{
  \item{MunsellSpec}{ standard Munsell specification or a Munsell vector in 
  ColorLab format }
}
\value{
  \item{sRGB}{ sRGB coordinates}
  \item{Status.ind}{ status of the operation}
  \item{OutOfGamutFlag}{ logical, TRUE = result is out of gamut}
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
MunsellTosRGB('7.6P 8.9/2.2')
}
\keyword{datasets}

