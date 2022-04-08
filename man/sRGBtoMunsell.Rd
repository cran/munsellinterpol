\name{sRGBtoMunsell}
\alias{sRGBtoMunsell}
\title{Convert sRGB coordinates to Munsell HVC}

\description{Converts non-linear sRGB coordinates to Munsell HVC, 
by interpolating over the extrapolated Munsell renotation data
}

\usage{sRGBtoMunsell( sRGB, maxSignal=255, ... ) }

\arguments{
\item{sRGB}{ a numeric Nx3 matrix with signal sRGB coordinates in the rows,
or a vector that can be converted to such a matrix, by row.
These are non-linear signal values, but they are not required to be integers.}

\item{maxSignal}{maximum value of signal sRGB.
Other popular values are 1, 100, 1023, and 65535.
Even when 1, they are still taken to be non-linear signal values.}

\item{...}{other parameters passed to \code{\link{XYZtoMunsell}()}    }
}

\value{
a numeric Nx3 matrix with HVC coordinates in the rows.
The rownames are copied from input to output, unless the rownames are \code{NULL}
when they are set to the Munsell notations for HVC.\cr
In case of error, it returns \code{NULL}.
}


\details{
The conversion is done in these steps:
\itemize{
\item sRGB  \rarrow  XYZ using \code{\link[spacesRGB:XYZfromRGB]{spacesXYZ::XYZfromRGB}()} with the given \code{maxSignal}

\item XYZ is adapted from Illuminant D65 (from the sRGB standard) to Illuminant C,
using \code{\link[spacesXYZ:adaptXYZ]{spacesXYZ::adaptXYZ}()} with the \emph{Bradford} chromatic adaptation method

\item XYZ  \rarrow  HVC using \code{\link{XYZtoMunsell}()} and the given \code{...}
}
}

\note{
The more general function \code{\link{RGBtoMunsell}()} also performs this conversion.
This function has the advantage that it takes a little less time,
since the CAT (using the \emph{Bradford} method) is precomputed during
\code{\link[base:.onLoad]{base::.onLoad}()}.
But it has the disadvantage that the chromatic adaptation method cannot be changed.
}


\references{
Wikipedia. \bold{sRGB}.
\url{https://en.wikipedia.org/wiki/SRGB}.

Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{https://www.munsellcolourscienceforpainters.com/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{RGBtoMunsell}()},
\code{\link{XYZtoMunsell}()},
\code{\link[spacesXYZ:CAT]{spacesXYZ::CAT}()},
\code{\link[spacesXYZ:adaptXYZ]{spacesXYZ::adaptXYZ}()},
\code{\link[spacesRGB:XYZfromRGB]{spacesXYZ::XYZfromRGB}()}
}

\examples{
sRGBtoMunsell( c(255,45,67) )
##                     H        V        C
##  5.4R 5.5/18 5.401135 5.477315 18.01984

sRGBtoMunsell( c(1,0,1), maxSignal=1 )
##                 H        V        C
##  8P 6/26 87.98251 5.981297 25.64534
}

\keyword{sRGB}

