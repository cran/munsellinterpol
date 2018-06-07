\name{MunsellTosRGB}
\alias{MunsellTosRGB}
\title{Convert a Munsell specification to sRGB coordinates}
\description{\code{MunsellTosRGB} Converts a 
Munsell specification to sRGB coordinates, by interpolating over the 
extrapolated Munsell renotation data
}

\usage{MunsellTosRGB( MunsellSpec, maxValue=255, adaption='bradford', ... ) }
\arguments{
  \item{MunsellSpec}{ a numeric Nx3 matrix with HVC values in the rows, or a vector that can be converted to such a matrix, by row.\cr
\code{MunsellSpec} can also be a character vector with Munsell Notations,
which is converted to an Nx3 matrix using \code{\link{HVCfromMunsellName}()}. }
  \item{maxValue}{maximum of sRGB for display; Other popular values are 1, 1023, and 65535}  
  \item{adaption}{method for chromatic adaption, see \code{\link{adaption}} for valid values.
  Also see \bold{Details}. }

  \item{...}{other parameters passed to \code{\link{MunsellToxyY}()} }       
}

\value{
a \code{data.frame} with these columns
\item{xyY}{ an Nx3 matrix with xyY values in the rows, that are adapted to Illuminant C.
This is an intermediate result that is sometimes useful, e.g. it can be passed
to \code{\link{IsWithinMacAdamLimits}()}. }
\item{sRGB}{an Nx3 matrix with non-linear sRGB display values in the rows. All values are clamped to the appropriate cube, e.g. \eqn{[0,255]^3} }
\item{OutOfGamutFlag}{logical vector, TRUE means the result was out of gamut (the cube)
before clamping it}
}


\details{
The conversion is done in these steps.
\itemize{
\item HVC  \rarrow  xyY using \code{\link{MunsellToxyY}()}. This xyY is for Illuminant C.
\item xyY is adapted from Illuminant C to Illuminant D65 (from the sRGB standard) using
\code{\link{adaptxyY}()} and the given \code{\link{adaption}} method
\item xyY  \rarrow  XYZ using \code{\link{xyY2XYZ}()} 
\item XYZ  \rarrow  sRGB using \code{\link{xyz2srgb}()} with the given \code{maxValue}
}
}

\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{http://www.99main.com/~centore/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{labtoMunsell}()},
\code{\link{MunsellToXYZ}()},
\code{\link{xyz2srgb}()},
\code{\link{xyY2XYZ}()},
\code{\link{IsWithinMacAdamLimits}()}
}


\examples{
MunsellTosRGB('7.6P 8.9/2.2')
}
\keyword{datasets}

