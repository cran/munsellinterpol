\name{luv2xyz}
\alias{luv2xyz}
\title{Convert CIE Luv coordinates to XYZ}
\description{\code{luv2xyz} Converts CIE Luv coordinates to XYZ.
}
\usage{luv2xyz(C,bbb) }
\arguments{
  \item{C}{ CIE Luv coordinates}
  \item{bbb}{ Reference white, by default Observer = 2deg, Illuminant = D65}
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
luv2xyz(c(0.310897, 0.306510, 74.613450))
}
\keyword{datasets}

