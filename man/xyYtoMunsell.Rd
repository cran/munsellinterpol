\name{xyYtoMunsell}
\alias{xyYtoMunsell}

\title{Convert xyY coordinates into a Munsell specification}

\description{\code{xyYtoMunsell} Converts xyY coordinates into a Munsell 
specification, by interpolating over the extrapolated Munsell renotation data
}

\usage{xyYtoMunsell( xyY, xyC='NBS', hcinterp='bicubic', vinterp='cubic',
                     VfromY='ASTM', warn=TRUE, perf=FALSE ) }

\arguments{
\item{xyY}{a numeric Nx3 matrix with CIE xyY coordinates in the rows,
or a vector that can be converted to such a matrix, by row.  
These are for viewing in an environment with Illuminant C, with Y=100.}

\item{xyC}{ a numeric 2-vector with xy chromaticity of Illuminant C.
It can also be one of the strings given in \code{\link{MunsellToxyY}()}. }

\item{hcinterp}{either \code{'bicubic'} or \code{'bilinear'} (partial matching enabled).
See \code{\link{MunsellToxyY}()} for details. }

\item{vinterp}{either \code{'cubic'} or \code{'linear'} (partial matching enabled).
See \code{\link{MunsellToxyY}()} for details. }

\item{VfromY}{passed as the parameter \code{which} to the function \code{VfromY()}.
See \code{\link{VfromY}()} for details.
Option \code{'MGO'} is not allowed because then Y>100 when V=10.
}

\item{warn}{if an xyY cannot be mapped
(usually because the root finder has wandered afar),
its H and V are set to \code{NA} in the returned \code{data.frame}.
Just before returning, if any rows have \code{NA}, and this argument is \code{TRUE},
then a warning is logged.
}

\item{perf}{if \code{perf} is \code{TRUE}, then extra performance related metrics
are appended to the returned \code{data.frame}, see \bold{Value}.
}
}

\value{
a \code{data.frame} with N rows and these columns:
\item{xyY}{The input xyY}
\item{HVC}{the computed HVC. H is automatically wrapped to (0,100].  In case of failure, H and C are \code{NA}.}
\item{SAMPLE_NAME}{a character vector - the Munsell notation for HVC}
If \code{perf} is \code{TRUE} then there are these additional columns:
\item{time.elapsed}{elapsed time in seconds. If available, the function \code{microbenchmark::get_nanotime()} is used.}
\item{iterations}{the number of iterations of \code{rootSolve::multiroot()}}
\item{evalations}{the number of forward (HVC  \rarrow  xyY) function evaluations}
\item{estim.precis}{the estimated precision from \code{rootSolve::multiroot()}.
This is in the HC plane for the Munsell Value easily computed from Y.}

If the rownames of xyY are not \code{NULL} and have no duplicates,
they are copied to the returned data frame.
}

\details{
See \code{\link{MunsellToxyY}()} and the \bold{User Guide - Appendix C}.
}

\source{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{https://www.munsellcolourscienceforpainters.com/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}\cr
\url{https://www.rit.edu/science/munsell-color-lab}\cr
\url{https://www.rit-mcsl.org/MunsellRenotation/all.dat}\cr
\url{https://www.rit-mcsl.org/MunsellRenotation/real.dat}
}

\section{Warning}{
Even when \code{vinterp='cubic'} the function xyY  \rarrow HVC is not \eqn{C^1}
on the plane V=1.
This is because of a change in Value spacing:
when V\eqn{\ge}1 the Value spacing is 1, but when V\eqn{\le}1 the Value spacing is 0.2.
}

\references{
Judd, Deane B.
\bold{The 1931 I.C.I. Standard Observer and Coordinate System for Colorimetry}.
Journal of the Optical Society of America.
Vol. 23. pp. 359-374.
October 1933.

Newhall, Sidney M.,
Dorothy Nickerson,
Deane B. Judd.
\bold{Final Report of the O.S.A. Subcommitte on the Spacing of the Munsell Colors}.
Journal of the Optical Society of America.
Vol. 33.  No. 7.  pp. 385-418.
July 1943.

Kelly, Kenneth L.
Kasson S. Gibson.
Dorothy Nickerson.
\bold{Tristimulus Specification of the Munsell Book of Color from Spectrophometric
Measurements}
National Bureau of Standards RP1549
Volume 31. August 1943.

Judd, Deane B. and Günther Wyszecki.
\bold{Extension of the Munsell Renotation System to Very Dark Colors}.
Journal of the Optical Society of America.
Vol. 46.  No. 4.  pp. 281-284.
April 1956.

Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{https://www.munsellcolourscienceforpainters.com/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}


\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{MunsellToxyY}()},
\code{\link[rootSolve]{multiroot}()},
\code{\link[microbenchmark]{get_nanotime}()}
}


\examples{
xyYtoMunsell(c(0.310897, 0.306510, 74.613450))
##       xyY.1     xyY.2     xyY.3     HVC.H     HVC.V     HVC.C    SAMPLE_NAME
## 1  0.310897  0.306510 74.613450 87.541720  8.900000  2.247428 7.5P 8.9/2.2
}

\keyword{datasets}
