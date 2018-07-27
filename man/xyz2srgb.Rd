\name{xyz2srgb}
\alias{xyz2srgb}
\title{Convert XYZ coordinates to sRGB}
\description{\code{xyz2srgb} Converts XYZ coordinates to sRGB.
This function is here only because package \pkg{colorscience} requires it.
A better function - \code{RGBfromXYZ()} - is now in \pkg{spacesRGB}.
}
\usage{xyz2srgb(XYZ) }
\arguments{
  \item{XYZ}{ XYZ coordinates}
}
\value{
  sRGB coordinates
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

\seealso{
\code{\link{RGBfromXYZ}()}
}

\examples{
xyz2srgb(c(0.310897, 0.306510, 74.613450))
}
\keyword{datasets}

