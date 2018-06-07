\name{xyY2XYZ}
\alias{xyY2XYZ}
\title{Convert xyY coordinates to XYZ}

\description{\code{xyY2XYZ} Converts xyY coordinates to XYZ.
}
\usage{xyY2XYZ( xyY ) }

\arguments{
  \item{xyY}{ a numeric Nx3 matrix with CIE xyY coordinates in the rows, or a vector that can be converted to such a matrix, by row.}
}

\value{
a numeric Nx3 matrix with XYZ coordinates in the rows.
The rownames are copied from input to output.
}

\note{
An input Y=0 is treated as a special case.
In this case x and y are ignored, and all XYZ are set to 0.
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
\code{\link{XYZ2xyY}()}
}

\examples{
xyY2XYZ(c(0.310897, 0.306510, 74.613450))
}

\keyword{XYZ}

