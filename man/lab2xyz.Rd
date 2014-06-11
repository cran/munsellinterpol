\name{lab2xyz}
\alias{lab2xyz}
\title{Convert CIE Lab coordinates to XYZ}
\description{\code{lab2xyz} Converts CIE Lab coordinates to XYZ.
}
\usage{lab2xyz(C,bbb) }
\arguments{
  \item{C}{ CIE Lab coordinates}
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
lab2xyz(c(17,12,12))
}
\keyword{datasets}

