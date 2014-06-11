\name{xyz2luv}
\alias{xyz2luv}
\title{Convert XYZ coordinates to CIE Luv}
\description{\code{xyz2luv} Converts XYZ coordinates to CIE Luv.
}
\usage{xyz2luv(XYZ,bbb) }
\arguments{
  \item{XYZ}{ XYZ coordinates}
  \item{bbb}{ Reference white, by default Observer = 2deg, Illuminant = D65}
}
\value{
  CIE Luv coordinates
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
xyz2luv(c(0.310897, 0.306510, 74.613450))
}
\keyword{datasets}

