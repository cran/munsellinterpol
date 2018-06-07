\name{xyz2luv}
\alias{xyz2luv}
\title{Convert XYZ coordinates to CIE Luv}

\description{Convert XYZ coordinates to CIE Luv}

\usage{xyz2luv( XYZ, white=c(95.047,100,108.883) ) }
\arguments{
\item{XYZ}{a numeric Nx3 matrix with CIE XYZ coordinates in the rows, or a vector that can be converted to such a matrix, by row.}
\item{white}{XYZ for the reference white.  The default is Illuminant D65, 2-degree observer}
}
\value{
An Nx3 matrix with the CIE Luv coordinates in the rows.
The rownames of the input are copied to the output.
}
 
\references{
CIE 15: Technical Report: Colorimetry.
CIE 15:2004.
3rd edition.
ISBN 3901906339.

Luv Color Space.
\url{https://en.wikipedia.org/wiki/Luv_color_space}
}


\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{luv2xyz}()}
}

\examples{
xyz2luv(  c(90.310897, 99.306510, 74.613450) )
}

\keyword{XYZ}

