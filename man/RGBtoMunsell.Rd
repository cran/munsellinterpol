\name{RGBtoMunsell}
\alias{RGBtoMunsell}
\title{Convert RGB coordinates to Munsell HVC}

\description{\code{RGBtoMunsell} Converts RGB coordinates to Munsell HVC,
by interpolating over the extrapolated Munsell renotation data
}

\usage{RGBtoMunsell( RGB, space='sRGB', which='scene', maxSignal=255, adapt='Bradford', ... ) }

\arguments{
\item{RGB}{ a numeric Nx3 matrix with RGB coordinates in the rows,
or a vector that can be converted to such a matrix, by row.
These are non-linear display values, but they are not required to be integers.}

\item{space}{the name of an installed RGB space.
Spaces \code{'sRGB'} and \code{'AdobeRGB'} are pre-installed,
and others can be installed using
\code{\link[spacesRGB:installRGB]{spacesRGB::installRGB}()}.}

\item{which}{either \code{'scene'} or \code{'display'}; passed to \code{\link[spacesRGB:XYZfromRGB]{spacesRGB::XYZfromRGB}()} }

\item{maxSignal}{maximum value of signal RGB.
Other popular values are 1, 100, 1023, and 65535.
Even when 1, they are still taken to be non-linear signal values.}

\item{adapt}{method for chromatic adaptation,
see \code{\link[spacesXYZ:CAT]{spacesXYZ::CAT}()} for valid values.
Also see \bold{Details}.}
  
\item{...}{other parameters passed to \code{\link{XYZtoMunsell}()} }
}

\value{
a numeric Nx3 matrix with HVC coordinates in the rows.
The rownames are copied from input RGB to output HVC, unless the rownames are \code{NULL}
when they are set to the Munsell notations for HVC.\cr
In case of error, it returns \code{NULL}.
}


\details{
The conversion is done in these steps:
\itemize{
\item RGB  \rarrow  XYZ using \code{\link[spacesRGB:XYZfromRGB]{spacesRGB::XYZfromRGB}()} with the given
\code{space}, \code{which}, and \code{maxSignal}

\item XYZ is adapted from the white-point (with \code{which}) of \code{space} to Illuminant C,
using \code{\link[spacesXYZ:adaptXYZ]{spacesXYZ::adaptXYZ}()},
with the given chromatic adaptation method \code{adapt}

\item XYZ  \rarrow  HVC using \code{\link{XYZtoMunsell}()} with \code{...}
}
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
\code{\link{XYZtoMunsell}()},
\code{\link[spacesXYZ:CAT]{spacesXYZ::CAT}()},
\code{\link[spacesXYZ:adaptXYZ]{spacesXYZ::adaptXYZ}()},
\code{\link[spacesRGB:XYZfromRGB]{spacesRGB::XYZfromRGB}()}
}

\examples{
RGBtoMunsell( c(255,45,67) )
##                     H        V        C
##  5.4R 5.5/18 5.401135 5.477315 18.01984

RGBtoMunsell( c(255,45,67), space='Adobe' )
##                     H        V        C
##  5.9R 6.2/22 5.924749 6.214155 21.83907
}

\keyword{RGB}

