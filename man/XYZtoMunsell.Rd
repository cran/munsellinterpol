\name{XYZtoMunsell}
\alias{XYZtoMunsell}

\title{Convert XYZ coordinates to Munsell HVC}

\description{\code{XYZtoMunsell} Convert XYZ coordinates to Munsell HVC, by interpolating over the  Munsell renotation data
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
XYZtoMunsell( c(0.310897, 0.306510, 74.613450) )
##                         H         V        C
##  0.51PB 0.27/9.9 70.50874 0.2671142 9.916039
}

\keyword{XYZ}

