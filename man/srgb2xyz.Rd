\name{srgb2xyz}
\alias{srgb2xyz}
\title{Convert sRGB coordinates to XYZ}
\description{\code{srgb2xyz} Converts sRGB coordinates to XYZ.
}
\usage{srgb2xyz(RGBmatrix) }
\arguments{
  \item{RGBmatrix}{ sRGB coordinates}
}
\value{
  CIE XYZ coordinates
}
\source{
Logicol S.r.l., 2014
EasyRGB color search engine
\url{http://www.easyrgb.com/}
}
\references{
Logicol S.r.l., 2014
EasyRGB color search engine
\url{http://www.easyrgb.com/}
}
\author{Jose Gama}
\examples{
srgb2xyz(c(0.310897, 0.306510, 74.613450))
}
\keyword{datasets}

