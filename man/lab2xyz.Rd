\name{lab2xyz}
\alias{lab2xyz}
\title{Convert CIE Lab coordinates to XYZ}
\description{\code{lab2xyz} Converts CIE Lab coordinates to XYZ.
}
\usage{lab2xyz( Lab, white=c(95.047,100,108.883) ) }
\arguments{
\item{Lab}{a numeric Nx3 matrix with CIE Lab coordinates in the rows, or a vector that can be converted to such a matrix, by row.}
\item{white}{ XYZ for the reference white.  The default is Illuminant=D65, 2 observer }
}
\value{
An Nx3 matrix with the CIE XYZ coordinates in the rows.
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
\code{\link{xyz2lab}()}
}

\examples{
lab2xyz( c(17,-12,12) )
}
\keyword{datasets}

