\name{xyz2lab}
\alias{xyz2lab}
\title{Convert XYZ coordinates to CIE lab}
\description{\code{xyz2lab} Converts XYZ coordinates to CIE lab.
}
\usage{xyz2lab(C,bbb) }
\arguments{
  \item{C}{ XYZ coordinates}
  \item{bbb}{ Reference white, by default Observer = 2deg, Illuminant = D65}
}
\value{
  CIE lab coordinates
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
xyz2lab(c(0.310897, 0.306510, 74.613450))
}
\keyword{datasets}

