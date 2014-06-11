\name{MunsellToLuv}
\alias{MunsellToLuv}
\title{Convert a Munsell specification into CIE Luv coordinates}
\description{\code{MunsellToLuv} Converts a 
Munsell specification into CIE Luv coordinates, by interpolating over the 
extrapolated Munsell renotation data
}
\usage{MunsellToLuv(MunsellSpec, InterpolateByLuminanceFactor=TRUE) }
\arguments{
  \item{MunsellSpec}{ standard Munsell specification or a Munsell vector in 
  ColorLab format }
  \item{InterpolateByLuminanceFactor}{ logical, TRUE = interpolation by 
  luminance factor }
}
\value{
  \item{L, u, v}{ CIE Lab coordinates of MunsellSpec illuminated by Illuminant C}
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
MunsellToLuv('7.6P 8.9/2.2')
}
\keyword{datasets}

