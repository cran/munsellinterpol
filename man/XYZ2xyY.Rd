\name{XYZ2xyY}
\alias{XYZ2xyY}
\title{Convert XYZ coordinates to xyY}
\description{\code{XYZ2xyY} Converts XYZ coordinates to xyY.
}

\usage{XYZ2xyY( XYZ ) }

\arguments{
  \item{XYZ}{ a numeric Nx3 matrix with CIE XYZ coordinates in the rows, or a vector that can be converted to such a matrix, by row.}
}

\value{
a numeric Nx3 matrix with xyY coordinates in the rows.
If the sum of XYZ is not positive, x and y are set to \code{NA}.
The rownames are copied from input to output.  
}


\references{
Wikipedia.
CIE 1931 color space.
\url{https://en.wikipedia.org/wiki/CIE_1931_color_space}

Logicol S.r.l., 2014
EasyRGB color search engine
\url{http://www.easyrgb.com/}
}

\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{xyY2XYZ}()}
}


\examples{
XYZ2xyY(c(0.310897, 0.306510, 0.3))
}

\keyword{XYZ}

