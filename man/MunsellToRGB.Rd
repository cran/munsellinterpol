\name{MunsellToRGB}
\alias{MunsellToRGB}
\title{Convert a Munsell specification to RGB coordinates}
\description{\code{MunsellToRGB} Converts a 
Munsell specification to RGB coordinates, by interpolating over the 
extrapolated Munsell renotation data
}

\usage{MunsellToRGB( MunsellSpec, space='sRGB', maxSignal=255, adapt='Bradford', ... ) }
\arguments{
\item{MunsellSpec}{ a numeric Nx3 matrix with HVC values in the rows, or a vector that can be converted to such a matrix, by row.\cr
\code{MunsellSpec} can also be a character vector with Munsell Notations,
which is converted to an Nx3 matrix using \code{\link{HVCfromMunsellName}()}. }

\item{space}{the name of an installed RGB space.
Spaces \code{'sRGB'} and \code{'AdobeRGB'} are pre-installed,
and others can be installed with \code{\link{installRGB}()}.
}

\item{maxSignal}{maximum of the non-linear signal RGB; Other popular values are 1, 1023, and 65535}

\item{adapt}{method for chromatic adaptation, see \code{\link[spacesXYZ:adaptation]{CAT}()} for valid values.
Also see \bold{Details}. }

\item{...}{other parameters passed to \code{\link{MunsellToxyY}()} }       
}

\value{
a \code{data.frame} with these columns
\item{xyY}{ an Nx3 matrix with xyY values in the rows, that are adapted to Illuminant C.
This is an intermediate result that is sometimes useful, e.g. it can be passed
to \code{\link{IsWithinMacAdamLimits}()}. }
\item{RGB}{an Nx3 matrix with non-linear RGB signal values in the rows. All values are clamped to the appropriate cube, e.g. \eqn{[0,255]^3} }
\item{OutOfGamut}{logical vector, TRUE means the result was out of gamut (the cube)
before clamping it}
In case of error, it returns \code{NULL}.
}


\details{
The conversion is done in these steps.
\itemize{
\item HVC  \rarrow  xyY using \code{\link{MunsellToxyY}()}. This xyY is for Illuminant C.
\item xyY is adapted from Illuminant C to the white-point of the RGB
\code{space} using \code{\link[spacesXYZ:adapt]{adaptxyY}()} and the given chromatic adaptation method
\item xyY  \rarrow  XYZ using \code{\link[spacesXYZ:toXYZ]{XYZfromxyY}()} 
\item XYZ  \rarrow  RGB using \code{\link[spacesRGB]{RGBfromXYZ}()} with the given \code{space}
and \code{maxSignal}
}
}

\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{https://www.munsellcolourscienceforpainters.com/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{installRGB}()},
\code{\link{MunsellToXYZ}()},
\code{\link[spacesRGB]{RGBfromXYZ}()},
\code{\link[spacesXYZ:toXYZ]{XYZfromxyY}()},
\code{\link[spacesXYZ:adaptation]{CAT}()},
\code{\link{IsWithinMacAdamLimits}()}
}


\examples{
MunsellToRGB('7.6P 8.9/2.2')
}

\keyword{RGB}

