\name{RGBtoMunsell}
\alias{RGBtoMunsell}
\title{Convert sRGB coordinates to a Munsell specification}

\description{\code{RGBtoMunsell} Converts RGB coordinates to a Munsell 
specification, by interpolating over the extrapolated Munsell renotation data
}

\usage{RGBtoMunsell( RGB, space='sRGB', maxSignal=255, adapt='Bradford', ... ) }

\arguments{
\item{RGB}{ a numeric Nx3 matrix with RGB coordinates in the rows,
or a vector that can be converted to such a matrix, by row.
These are non-linear display values, but they are not required to be integers.}

\item{space}{the name of an installed RGB space.
Spaces \code{'sRGB'} and \code{'AdobeRGB'} are pre-installed,
and others can be installed with \code{\link{installRGB}()}.}

\item{maxSignal}{maximum value of RGB for display.
Other popular values are 1, 1023, and 65535.
Even when 1, they are still taken to be non-linear display values.}

\item{adapt}{method for chromatic adaptation,
see \code{\link[spacesXYZ:adaptation]{CAT}()} for valid values.
Also see \bold{Details}. }
  
\item{...}{other parameters passed to \code{\link{XYZtoMunsell}()}    }
}

\value{
a numeric Nx3 matrix with HVC coordinates in the rows.
The rownames are copied from input to output.\cr
In case of error, it returns \code{NULL}.
}


\details{
The conversion is done in 3 steps.
\itemize{
\item RGB  \rarrow  XYZ using \code{\link{XYZfromRGB}()} with the given
\code{space} and \code{maxSignal}
\item XYZ is adapted from the white-point of  \code{space} to Illuminant C
using the given chromatic adaptation method
\item XYZ  \rarrow  HVC using \code{\link{XYZtoMunsell}()}
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
\code{\link[spacesRGB]{XYZfromRGB}()},
\code{\link{XYZtoMunsell}()},
\code{\link[spacesXYZ:adaptation]{CAT}()}
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

