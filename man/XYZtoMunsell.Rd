\name{XYZtoMunsell}
\alias{XYZtoMunsell}

\title{Convert XYZ coordinates to a Munsell specification}

\description{\code{XYZtoMunsell} Convert XYZ coordinates to a Munsell 
specification, by interpolating over the extrapolated Munsell renotation data
}

\usage{XYZtoMunsell( XYZ, ... ) }
\arguments{
\item{XYZ}{a numeric Nx3 matrix with CIE XYZ coordinates in the rows,
or a vector that can be converted to such a matrix, by row.
The XYZ are for viewing in an environment with Illuminant C.}

\item{...}{other parameters passed to \code{\link{xyYtoMunsell}()} }
}

\value{
an Nx3 matrix with Munsell HVC in the rows.
The rownames are copied from input to output.
}

\details{
the function calls \code{XYZ2xyY()} and \code{xyYtoMunsell()}.
}

\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{https://www.munsellcolourscienceforpainters.com/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{MunsellToXYZ}()}
}

\examples{
XYZtoMunsell(c(0.310897, 0.306510, 74.613450))
}

\keyword{XYZ}

