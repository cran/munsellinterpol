\name{MunsellTosRGB}
\alias{MunsellTosRGB}

\title{Convert a Munsell specification to sRGB coordinates}

\description{\code{MunsellTosRGB} Converts a 
Munsell specification to non-linear sRGB coordinates, by interpolating over the 
extrapolated Munsell renotation data
}

\usage{MunsellTosRGB( MunsellSpec, maxSignal=255, ... ) }
\arguments{
  \item{MunsellSpec}{ a numeric Nx3 matrix with HVC values in the rows, or a vector that can be converted to such a matrix, by row.\cr
\code{MunsellSpec} can also be a character vector with Munsell Notations,
which is converted to an Nx3 matrix using \code{\link{HVCfromMunsellName}()}. }
  \item{maxSignal}{maximum of signal sRGB; Other popular values are 1, 1023, and 65535}  

  \item{...}{other parameters passed to \code{\link{MunsellToxyY}()} }       
}

\value{
a \code{data.frame} with these columns
\item{xyY}{ an Nx3 matrix with xyY values in the rows, that are adapted to Illuminant C.
This is an intermediate result that is sometimes useful, e.g. it can be passed
to \code{\link{IsWithinMacAdamLimits}()}. }
\item{RGB}{an Nx3 matrix with non-linear signal sRGB values in the rows. All values are clamped to the appropriate cube, e.g. \eqn{[0,255]^3} }
\item{OutOfGamut}{logical vector, TRUE means the result was out of gamut (the cube)
before clamping it}
}


\details{
The conversion is done in these steps.
\itemize{
\item HVC  \rarrow  xyY using \code{\link{MunsellToxyY}()}. This xyY is for Illuminant C.
\item xyY is adapted from Illuminant C to Illuminant D65 (from the sRGB standard) using
\code{\link{adaptxyY}()} and the \emph{Bradford} chromatic adaptation method (CAT)
\item xyY  \rarrow  XYZ  using \code{\link[spacesXYZ:toXYZ]{XYZfromxyY}()} 
\item XYZ  \rarrow  sRGB using \code{\link[spacesRGB]{RGBfromXYZ}()} with the given \code{maxSignal}
}
}

\note{
The function \code{\link{MunsellToRGB}()} also performs this conversion.
The only reason to use this one is that it takes a little less time,
since the CAT (using the \emph{Bradford} method) is precomputed.
}

\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{http://centore.isletech.net/~centore/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{MunsellToXYZ}()},
\code{\link{MunsellToRGB}()},
\code{\link[spacesRGB]{RGBfromXYZ}()},
\code{\link[spacesXYZ:toXYZ]{XYZfromxyY}()},
\code{\link[spacesXYZ:adaptation]{CAT}()},
\code{\link{IsWithinMacAdamLimits}()}
}


\examples{
MunsellTosRGB('7.6P 8.9/2.2')
##                SAMPLE_NAME      xyY.x      xyY.y      xyY.Y    RGB.R    RGB.G    RGB.B OutOfGamut
##  7.6P 8.9/2.2 7.6P 8.9/2.2  0.3109520  0.3068719 74.6134498 231.3575 221.1421 230.3501      FALSE
}
\keyword{datasets}

