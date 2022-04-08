\name{MunsellToRGB}
\alias{MunsellToRGB}

\title{Convert a Munsell specification to RGB coordinates}
\description{\code{MunsellToRGB} Converts a 
Munsell specification to RGB coordinates, by interpolating over the 
extrapolated Munsell renotation data
}

\usage{MunsellToRGB( MunsellSpec, space='sRGB', which='scene', maxSignal=255, 
    adapt='Bradford', ... ) }
\arguments{
\item{MunsellSpec}{ a numeric Nx3 matrix with HVC values in the rows, or a vector that can be converted to such a matrix, by row.\cr
\code{MunsellSpec} can also be a character vector with Munsell Notations,
which is converted to an Nx3 matrix using \code{\link{HVCfromMunsellName}()}. }

\item{space}{the name of an installed RGB space.
Spaces \code{'sRGB'} and \code{'AdobeRGB'} are pre-installed,
and others can be installed with \code{\link[spacesRGB:installRGB]{spacesRGB::installRGB}()}.
}

\item{which}{either \code{'scene'} or \code{'display'}; passed to \code{\link[spacesRGB:RGBfromXYZ]{spacesRGB::RGBfromXYZ}()} }

\item{maxSignal}{maximum of the non-linear signal RGB; Other popular values are 1, 100, 1023, and 65535}

\item{adapt}{method for chromatic adaptation, see \code{\link[spacesXYZ:CAT]{spacesXYZ::CAT}()} for valid values.
Also see \bold{Details}. }

\item{...}{other parameters passed to \code{\link{MunsellToxyY}()} }       
}

\value{
a \code{data.frame} with these columns:

\item{SAMPLE_NAME}{a character N-vector equal to the original \code{MunsellSpec} if that was a character vector.
Or the Munsell notation string converted from the input matrix HVC.}

\item{xyY}{ an Nx3 matrix with xyY values in the rows, that are adapted to Illuminant C.
This is an intermediate result that is sometimes useful, e.g. it can be passed
to \code{\link{IsWithinMacAdamLimits}()}.
The rownames of xyY are copied from those of HVC, unless the rownames of HVC are \code{NULL}
when they are set to \code{SAMPLE_NAME}. }

\item{RGB}{an Nx3 matrix with non-linear RGB signal values in the rows. 
All values are clamped to the appropriate cube, e.g. \eqn{[0,255]^3}.
The rownames of RGB are copied from those of xyY. }

\item{OutOfGamut}{logical vector, TRUE means the result was out of gamut (the cube) before clamping it}

In case of error, it returns \code{NULL}.
}


\details{
The conversion is done with these steps:
\itemize{
\item HVC  \rarrow  xyY using \code{\link{MunsellToxyY}()} with \code{...}. This xyY is for Illuminant C.

\item xyY  \rarrow  XYZ using \code{\link[spacesXYZ:toXYZ]{spacesXYZ::XYZfromxyY}()} 

\item XYZ is adapted from Illuminant C to the white-point (with \code{which}) of the RGB \code{space},
using \code{\link[spacesXYZ:adaptXYZ]{spacesXYZ::adaptXYZ}()}, with the given chromatic adaptation method \code{adapt}

\item XYZ  \rarrow  RGB using \code{\link[spacesRGB:RGBfromXYZ]{spacesRGB::RGBfromXYZ}()} with the given \code{space},
\code{which}, and \code{maxSignal}
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
\code{\link{IsWithinMacAdamLimits}()},
\code{\link[spacesRGB:installRGB]{spacesRGB::installRGB}()},
\code{\link[spacesRGB:RGBfromXYZ]{spacesRGB::RGBfromXYZ}()},
\code{\link[spacesXYZ:toXYZ]{spacesXYZ::XYZfromxyY}()},
\code{\link[spacesXYZ:CAT]{spacesXYZ::CAT}()}
}

\examples{
MunsellToRGB( c('7.6P 8.9/2.2', 'N 3/'), space='AdobeRGB' )
##     SAMPLE_NAME      xyY.x      xyY.y      xyY.Y     RGB.R     RGB.G     RGB.B OutOfGamut
##  1 7.6P 8.9/2.2  0.3109520  0.3068719 74.6134498 227.72419 220.18659 229.23297      FALSE
##  2         N 3/  0.3101000  0.3163000  6.3911778  73.01793  73.01793  73.01793      FALSE

}

\keyword{RGB}

