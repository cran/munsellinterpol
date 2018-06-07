\name{srgb2xyz}
\alias{srgb2xyz}
\title{Convert sRGB coordinates to XYZ}

\description{Convert sRGB coordinates to XYZ}

\usage{srgb2xyz( RGBmatrix, maxValue=255 ) }
\arguments{
\item{RGBmatrix}{ a numeric Nx3 matrix with sRGB coordinates in the rows, or a vector that can be converted to such a matrix, by row.
These are non-linear display values, but they are not required to be integers.
}
\item{maxValue}{maximum value of sRGB for display.
Other popular values are 1, 1023, and 65535.
Even when 1, they are still taken to be non-linear display values.}
}

\value{
a numeric Nx3 matrix with XYZ coordinates in the rows.
These are for viewing under Illuminant D65 (from the sRGB standard), with Y=100.
}

\references{
Wikipedia. \bold{sRGB}.
\url{https://en.wikipedia.org/wiki/SRGB}.

Logicol S.r.l., 2014
EasyRGB color search engine
\url{http://www.easyrgb.com/}
}

\author{Jose Gama and Glenn Davis}


\seealso{
\code{\link{xyz2srgb}()}
}

\examples{
srgb2xyz( c(128,200,255, 0,0,0,  255,255,255) )
##              X         Y        Z
##  [1,] 47.60334  53.11601 102.3549
##  [2,]  0.00000   0.00000   0.0000
##  [3,] 95.04559 100.00000 108.9058
}


\keyword{sRGB}

