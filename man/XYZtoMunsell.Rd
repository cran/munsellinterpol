\name{XYZtoMunsell}
\alias{XYZtoMunsell}

\title{Convert XYZ coordinates to Munsell HVC}

\description{\code{XYZtoMunsell} Convert XYZ coordinates to Munsell HVC, by interpolating over the  Munsell renotation data
}

\usage{XYZtoMunsell( XYZ, ... ) }
\arguments{
\item{XYZ}{a numeric Nx3 matrix with CIE XYZ coordinates in the rows,
or a vector that can be converted to such a matrix, by row.
The XYZ are for viewing in an environment with Illuminant C.
XYZ must be scaled so that Y=100 for the perfect reflecting diffuser.
}

\item{...}{other parameters passed to \code{\link{xyYtoMunsell}()} }
}

\value{
an Nx3 matrix with Munsell HVC in the rows.
The rownames are copied from input to output, unless the rownames are \code{NULL}
when they are set to the Munsell notations for HVC.\cr
In case of error, it returns \code{NULL}.
}

\details{
the function calls \code{XYZ2xyY()} and \code{xyYtoMunsell()}.

The conversion is done in these steps:
\itemize{
\item XYZ  \rarrow  xyY using \code{\link[spacesXYZ:xyYfromXYZ]{spacesXYZ::xyYfromXYZ}()} 
\item xyY  \rarrow  HVC using \code{\link{xyYtoMunsell}()} and \code{...}
}
}

\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{https://www.munsellcolourscienceforpainters.com/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{MunsellToXYZ}()},
\code{\link[spacesXYZ:xyYfromXYZ]{spacesXYZ::xyYfromXYZ}()} 
}

\examples{
XYZtoMunsell( c(31.0897, 30.6510, 74.613450) )
##                      H        V        C
## 5.4PB 6.1/9.4 75.43786 6.117631 9.416488
}

\keyword{XYZ}

