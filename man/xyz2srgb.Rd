\name{xyz2srgb}
\alias{xyz2srgb}
\title{Convert XYZ coordinates to sRGB}

\description{Convert XYZ coordinates to sRGB}

\usage{xyz2srgb( XYZ, maxValue=255 ) }
\arguments{
\item{XYZ}{a numeric Nx3 matrix with CIE XYZ coordinates in the rows, or a vector that can be converted to such a matrix, by row.
The XYZ are for viewing in an environment with Illuminant D65 (from the sRGB standard)
with Y=100.}

\item{maxValue}{maximum value of sRGB for display.
Other popular values are 1, 1023, and 65535..
Even when 1, they are still taken to be non-linear display values.}
}

\value{
a \code{data.frame} with N rows and these columns
\item{sRGB}{non-linear sRGB coordinates. All values are clamped to the appropriate cube, e.g. \eqn{[0,255]^3}.
Values are not rounded.}
\item{OutOfGamutFlag}{logical vector, TRUE means the result was out of gamut (the cube)
before clamping it}
}


\references{
Wikipedia. \bold{sRGB}.
\url{https://en.wikipedia.org/wiki/SRGB}

Logicol S.r.l., 2014
EasyRGB color search engine
\url{http://www.easyrgb.com/}
}

\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{srgb2xyz}()}
}

\examples{
xyz2srgb( c(80.310897, 90.306510, 84.613450) )
##      sRGB.R   sRGB.G   sRGB.B OutOfGamutFlag
##  1 230.1676 249.4122 225.2472          FALSE
}

\keyword{sRGB}

