\name{MunsellToxyY}
\alias{MunsellToxyY}
\title{Convert a Munsell specification into xyY coordinates}
\description{\code{MunsellToxyY} Converts a 
Munsell specification into xyY coordinates, by interpolating over the 
extrapolated Munsell renotation data
}
\usage{MunsellToxyY(MunsellSpec, InterpolateByLuminanceFactor=TRUE) }
\arguments{
  \item{MunsellSpec}{ standard Munsell specification or a Munsell vector in 
  ColorLab format }
  \item{InterpolateByLuminanceFactor}{ logical, TRUE = interpolation by 
  luminance factor }
}
\value{
  \item{x, y, Y}{ CIE xyY coordinates of MunsellSpec illuminated by Illuminant C}
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
MunsellToxyY('7.6P 8.9/2.2')
}
\keyword{datasets}

