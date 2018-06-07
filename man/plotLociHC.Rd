\name{plotLociHC}
\alias{plotLociHC}
\title{Plot Curves of Constant Munsell Hue and Chroma}

\description{Plot Curves of Constant Munsell Hue and Chroma}

\usage{plotLociHC( value=5, hue=seq(2.5,100,by=2.5), chroma='auto', coords='xy',
                        main="Value \%g/", est=FALSE, ... ) }

\arguments{
\item{value}{a Munsell value for which the plot is created.
It must be in the interval (0,10].
\code{value} can also be a numeric vector of such numbers,
and then a separate plot is made for each element of the vector.}

\item{hue}{a numeric vector for which curves of contant Hue are plotted.
Each of these radial curves starts at Munsell Chroma = \code{min(chroma,1)}
and extends to \code{max(chroma)}.
\code{hue} can also be a character vector,
which is then converted to a numeric vector using
\code{\link{HueNumberFromString}()}.
}

\item{chroma}{a numeric vector for which ovoids of constant Chroma are plotted.
Each of these ovoids is closed; i.e. goes full circle from Hue=0 to Hue=100.
If \code{chroma='auto'} then a vector is chosen appropriate for the current \code{value}.
}

\item{coords}{either \code{'xy'} or \code{'ab'}.
If  \code{coords='xy'} then the plots are in the standard xy chromaticity plane.
If  \code{coords='ab'} then the plots are in the ab chrominance plane, from Lab.
Even when it is \code{'ab'}, the interpolation of loci takes place in \code{'xy'}
\emph{before} transformation to Lab.
}

\item{main}{a string used to set the main title of the plot.
The optional placeholder \code{'\%g'} is replaced by the current \code{value}.
}

\item{est}{if \code{TRUE}, initial estimates for the iteration used in
\code{\link{xyYtoMunsell}()} are plotted
}

\item{...}{other arguments passed to the function \code{MunsellToxyY()}.
This includes \code{hcinterp}, \code{vinterp}, and \code{xyC}.
However \code{warn=FALSE} is forced.
}
}

\details{
The plot limits (\code{xlim} and \code{ylim}) are set to include
all points where the Hue radials intersect the Chroma ovoids,
plus the white point.
\cr
If \code{value} is one of {0.2,0.4,0.6,0.8,1,2,3,4,5,6,7,8,9,10}
then published points from \code{real.dat}
are plotted with filled black points (real points),
and extrapolated points from \code{all.dat} are drawn with open circles (unreal points).
}

\note{
The option \code{hcinterp='bicubic'} makes fairly smooth ovoids,
and \code{hcinterp='bilinear'} makes 40-sided polygons (when \code{coords='xy'}).
Compare with the plots in  
\emph{Newhall et. al. (1943)},
\emph{Judd, et. al. (1956)},
and \emph{Judd, et. al. (1975)} p. 263.
}

\value{
\code{TRUE} for success and \code{FALSE} for failure.
}


\references{
Newhall, Sidney M.,
Dorothy Nickerson,
Deane B. Judd.
\bold{Final Report of the O.S.A. Subcommitte on the Spacing of the Munsell Colors}.
Journal of the Optical Society of America.
Vol. 33.  No. 7.  pp. 385-418.
July 1943.

Judd, Deane B. and Günther Wyszecki.
\bold{Extension of the Munsell Renotation System to Very Dark Colors}.
Journal of the Optical Society of America.
Vol. 46.  No. 4.  pp. 281-284.
April 1956.

Judd, Deane B. and Günther Wyszecki.
\bold{Color in Business, Science, and Industry}.
3rd edition.
John Wiley & Sons.
1975.
}

\author{Glenn Davis}

\seealso{
\code{\link{MunsellToxyY}()},
\code{\link{HueNumberFromString}()}
}

\keyword{plot}

