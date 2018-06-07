\name{xyz2lab}
\alias{xyz2lab}
\title{Convert XYZ coordinates to CIE Lab}

\description{Convert XYZ coordinates to CIE Lab}

\usage{xyz2lab( XYZ, white=c(95.047,100,108.883) ) }
\arguments{
\item{XYZ}{a numeric Nx3 matrix with CIE XYZ coordinates in the rows, or a vector that can be converted to such a matrix, by row.}
\item{white}{XYZ for the reference white.  The default is Illuminant D65, 2-degree observer}
}

\value{
An Nx3 matrix with the CIE Lab coordinates in the rows.
The rownames of the input are copied to the output.
}


\references{
CIE 15: Technical Report: Colorimetry.
CIE 15:2004.
3rd edition.
ISBN 3901906339.

Lab Color Space.
\url{https://en.wikipedia.org/wiki/Lab_color_space}
}

\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{lab2xyz}()}
}

\examples{
xyz2lab( c(90.310897, 99.306510, 74.613450) )
}

\keyword{datasets}

