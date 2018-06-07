\name{luv2xyz}
\alias{luv2xyz}
\title{Convert CIE Luv coordinates to XYZ}
\description{\code{luv2xyz} Converts CIE Luv coordinates to XYZ.
}

\usage{luv2xyz( Luv, white=c(95.047,100,108.883) ) }
\arguments{
  \item{Luv}{CIE Luv coordinates. An Nx3 matrix, or a vector that can be converted to such a matrix. Each row of the matrix has Luv.}
  \item{white}{ XYZ for the reference white.  The default is Illuminant=D65, 2 observer }
}
\value{
An Nx3 matrix with the CIE XYZ coordinates. Each row of the matrix has XYZ.
The rownames are set to those of Luv.
}

\references{
CIE 15: Technical Report: Colorimetry.
CIE 15:2004.
3rd edition.
ISBN 3901906339.

Luv Color Space.
\url{https://en.wikipedia.org/wiki/Luv_color_space}
}

\author{Jose Gama and Glenn Davis}

\examples{
luv2xyz(c(55, -5.06510, 14.613450))
}

\keyword{datasets}

