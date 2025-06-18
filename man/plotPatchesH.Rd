\name{plotPatchesH}
\alias{plotPatchesH}
\title{Plot Colored Patches for a fixed Munsell Hue}

\description{This plot simulates a page from the Munsell Book of Color.
The colors are best viewed on a display calibrated for
the RGB space given as the second argument.}

\usage{plotPatchesH( hue, space='sRGB', adapt='Bradford', background='gray50',
                        main="Hue \%s  (H=\%g)      [\%s   adapt=\%s]",
                        value=NULL, chroma=NULL, ... ) }

\arguments{

\item{hue}{a Munsell hue for which the plot is created.
It is automatically wrapped to the interval (0,100].
It does not have to be a multiple of 2.5.
\code{hue} can also be a numeric vector of such numbers,
and then a separate plot is made for each element of the vector.
\code{hue} can also be a character vector of Hue Names,
which is then converted to a numeric vector using
\code{\link{HueNumberFromString}()}.
}

\item{space}{the name of an installed RGB space.
Spaces \code{'sRGB'} and \code{'AdobeRGB'} are pre-installed,
and others can be installed using \code{\link[spacesRGB:installRGB]{spacesRGB::installRGB}()} 
}

\item{adapt}{method used to adapt xyY for Illuminant C to xyY for Illuminant D65.
It is passed to \code{\link{MunsellToRGB}()}.
}

\item{background}{background color for the plot.
It is passed to \code{par()} as argument \code{bg}.
}

\item{main}{a string used to set the main title of the plot.
The optional placeholder \code{'\%s'} is replaced by the Hue Name,
\code{'\%g'} is replaced by the Hue Number,
the next \code{'\%s'} is replaced by \code{space}, and
the last \code{'\%s'} is replaced by \code{adapt}.
}

\item{value}{a vector of Munsell Values to use for the plot.
Values must be in the interval [0,10], and in strictly increasing order.
Non-integers, such as 2.5 and 8.5, are valid.
If \code{value=NULL}, then it is set to \code{0:10}.
}

\item{chroma}{a vector of Munsell Chromas to use for the plot.
Chromas must be non-negative and in strictly increasing order.
Odd integers and non-integers are valid.
If \code{chroma=NULL}, then a sequence of even numbers is chosen,
starting at 0 and ending at the Chroma limit for the given Hue page, see \bold{Details}.
}

\item{...}{other arguments passed to the function \code{MunsellToRGB()}.
This includes \code{hcinterp}, \code{vinterp}, and \code{xyC}.
}
}

\details{
If \code{chroma=NULL},
for the closest discrete Hue in \code{real.dat},
the patches in \code{real.dat} are transformed to xyY using simple lookup with no interpolation.
These are then tested against the MacAdam Limits for Illuminant C
using \code{\link{IsWithinMacAdamLimits}()}.
The patches outside the limits are discarded,
and the maximum Chroma of the remaining patches,
which is always an even integer, determines \code{chroma}.
\cr
Patches inside the MacAdam Limits can still be outside the RGB cube;
those patches are drawn in outline only, and with the clamped RGB coordinates printed inside.
}

\value{
\code{TRUE} for success and \code{FALSE} for failure.
}


\author{Glenn Davis}

\seealso{
\code{\link{MunsellToRGB}()},
\code{\link{HueNumberFromString}()},
\code{\link{IsWithinMacAdamLimits}()},
\code{\link[spacesRGB:installRGB]{spacesRGB::installRGB}()}
}

\keyword{plot}