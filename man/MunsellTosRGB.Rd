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
\code{MunsellSpec} can also be a character N-vector with Munsell notations,
which is converted to an Nx3 matrix using \code{\link{HVCfromMunsellName}()}. }
  \item{maxSignal}{maximum of signal sRGB; Other popular values are 1, 100, 1023, and 65535}

  \item{...}{other parameters passed to \code{\link{MunsellToxyY}()} }
}

\value{
a \code{data.frame} with these columns:

\item{SAMPLE_NAME}{the original \code{MunsellSpec} if that was a character vector.
Or the Munsell notation string converted from the input matrix HVC.}

\item{xyY}{ an Nx3 matrix with xyY values in the rows, that are adapted to Illuminant C.
This is an intermediate result that is sometimes useful, e.g. it can be passed
to \code{\link{IsWithinMacAdamLimits}()}.
The rownames of xyY are copied from those of HVC, unless the rownames of HVC are \code{NULL}
when they are set to \code{SAMPLE_NAME}. }

\item{RGB}{an Nx3 matrix with non-linear signal sRGB values in the rows.
All values are clamped to the appropriate cube, e.g. \eqn{[0,255]^3}.
The rownames of RGB are copied from those of xyY.}

\item{OutOfGamut}{logical vector, TRUE means the RGB result was out of gamut (the cube) before clamping it}
}


\details{
The conversion is done in these steps.
\itemize{
\item HVC  \rarrow  xyY using \code{\link{MunsellToxyY}()} and the given \code{...}. This xyY is for Illuminant C.

\item xyY  \rarrow  XYZ  using \code{\link[spacesXYZ:XYZfromxyY]{spacesXYZ::XYZfromxyY}()}

\item XYZ is adapted from Illuminant C to Illuminant D65 (from the sRGB standard) using
\code{\link[spacesXYZ:adaptXYZ]{spacesXYZ::adaptXYZ}()} and the \emph{Bradford} chromatic adaptation method

\item XYZ  \rarrow  sRGB using \code{\link[spacesRGB:RGBfromXYZ]{spacesRGB::RGBfromXYZ}()} with the given \code{maxSignal}
}
}

\note{
The more general function \code{\link{MunsellToRGB}()} also performs this conversion.
The main reason to use \code{MunsellTosRGB()} is that it takes a little less time,
since the CAT (using the \emph{Bradford} method) is precomputed during
\code{\link[base:.onLoad]{base::.onLoad}()}.
}

\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{https://www.munsellcolourscienceforpainters.com/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{MunsellToXYZ}()},
\code{\link{MunsellToRGB}()},
\code{\link[spacesXYZ:XYZfromxyY]{spacesXYZ::XYZfromxyY}()},
\code{\link[spacesXYZ:adaptation]{spacesXYZ::CAT}()},
\code{\link[spacesXYZ:adaptXYZ]{spacesXYZ::adaptXYZ}()},
\code{\link[spacesRGB:RGBfromXYZ]{spacesRGB::RGBfromXYZ}()},
\code{\link{IsWithinMacAdamLimits}()}
}


\examples{
MunsellTosRGB( c('7.6P 8.9/2.2', 'N 3/') )
##     SAMPLE_NAME      xyY.x      xyY.y      xyY.Y     RGB.R     RGB.G     RGB.B OutOfGamut
##  1 7.6P 8.9/2.2  0.3109520  0.3068719 74.6134498 231.35746 221.14207 230.35011      FALSE
##  2         N 3/  0.3101000  0.3163000  6.3911778  71.50491  71.50491  71.50491      FALSE
}
