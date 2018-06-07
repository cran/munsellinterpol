\name{MunsellToxyY}
\alias{MunsellToxyY}
\title{Convert a Munsell specification into xyY coordinates}
\description{\code{MunsellToxyY} Converts a 
Munsell specification into xyY coordinates, by interpolating over the 
extrapolated Munsell renotation data
}
\usage{MunsellToxyY( MunsellSpec, xyC='NBS', hcinterp='bicubic', vinterp='cubic',
                     YfromV='ASTM', warn=TRUE ) }

\arguments{
\item{MunsellSpec}{ a numeric Nx3 matrix or a vector that can be converted to
such a matrix.
Each row has Munsell HVC, where H is Hue Number,
and V and C are the standard Munsell Value and Chroma.
The Hue is automatically wrapped to the interval (0,100].
\code{MunsellSpec} can also be a character N-vector with standard Munsell notation; it is converted to an Nx3 matrix.}

\item{xyC}{ a numeric 2-vector with xy chromaticity of Illuminant C.
It can also be one of the strings in the first column of this table;
it is then replaced by the corresponding xy in the second column.
\tabular{lll}{
             \tab   xy white point \tab   reference \cr
\code{'NBS'} \tab  \code{c(0.3101,0.3163)} \tab Kelly, et. al. [RP1549] (1943). Rheinboldt et al. (1960) \cr
\code{'JOSA'} \tab  \code{c(0.31012,0.31631)} \tab Judd, Deane B. (1933) \cr
\code{'NTSC'} \tab  \code{c(0.310,0.316)}  \tab  NTSC (1953) \cr
\code{'CIE'}  \tab  \code{c(0.31006,0.31616)}  \tab  CIE:15 2004\cr
}
The default \code{'NBS'} is probably what is intended by \cite{Newhall et. al.}
although no xy for C appears in that paper.
This is the C used in the first computer program for conversion:
\emph{Rheinboldt et al. (1960)}.
The other options are provided so that a neutral Munsell chip has the xy that
the user expects.  Alternative values of \code{xyC} should not be \emph{too far}
from the above.
If \code{hcinterp} is \code{'bicubic'},
this parameter only affects chips with Chroma \eqn{<} 4 (except Chroma=2).
If \code{hcinterp} is \code{'bilinear'},
this parameter only affects chips with Chroma \eqn{<} 2.
}

\item{hcinterp}{either \code{'bicubic'} or \code{'bilinear'} (partial matching enabled).
In the bicubic case, for a general input point, the output value is interpolated
using a 4x4 subgrid of the lookup table, and the interpolation function is class \eqn{C^1}
(except at the neutrals).
In the bilinear case, the interpolation uses a 2x2 subgrid,
and the function is class \eqn{C^0}.
}

\item{vinterp}{either \code{'cubic'} or \code{'linear'} (partial matching enabled).
In the cubic case, for a general input point, the output value is interpolated
using 4 planes of constant Value, and the interpolation function is class \eqn{C^1}.
In the linear case, the interpolation uses 2 planes and the function is class \eqn{C^0}.
}

\item{YfromV}{passed as the parameter \code{which} to the function \code{YfromV()}.
See \code{\link{YfromV}()} for details.
Option \code{'MGO'} is not allowed because then Y>100 when V=10.
}

\item{warn}{if a chip cannot be mapped (usually because the Chroma is too large),
its x and y are set to \code{NA} in the returned
\code{data.frame}.
Just before returning, if any rows have \code{NA}, and this argument is \code{TRUE},
then a warning is logged.
}
}

\note{
When making plots in planes of constant Value,
option \code{hcinterp='bicubic'} makes fairly smooth ovals,
and \code{hcinterp='bilinear'} makes polygons.
The ovals are smooth even when \code{vinterp='linear'},
but the function is not class \eqn{C^1} at the planes of integer Value.
To get a fully \eqn{C^1} function (except at the neutrals and on the plane V=1),
\code{hcinterp} and \code{vinterp} must be set to the defaults.
}

\section{Warning}{
Even when \code{vinterp='cubic'} the function HVC \rarrow xyY is not \eqn{C^1}
on the plane V=1.
This is because of a change in Value spacing:
when V\eqn{\ge}1 the Value spacing is 1, but when V\eqn{\le}1 the Value spacing is 0.2.
}

\value{
a \code{data.frame} with these columns:
\item{SAMPLE_NAME}{the original \code{MunsellSpec} if that was a character vector.
Or the Munsell notation string converted from HVC.}
\item{HVC}{the input Nx3 matrix}
\item{xyY}{the computed output matrix, with CIE xyY coordinates of MunsellSpec illuminated by Illuminant C.  In case of error, x and y are set to \code{NA}.}
}

\details{
In case \code{hcinterp='bicubic'} or \code{vinterp='cubic'} a Catmull-Rom spline is used;
see the article \emph{Cubic Hermite spline}.
This spline has the nice property that it is \emph{local} and requires at most 4 points.
And if the knot spacing is uniform:
1) the resulting spline is \eqn{C^1},
2) if the knots are on a line, the interpolated points are on the line too.
}

\source{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{http://www.99main.com/~centore/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}\cr
\url{http://www.rit.edu/science/pocs/renotation}\cr
\url{http://www.rit-mcsl.org/MunsellRenotation/all.dat}\cr
\url{http://www.rit-mcsl.org/MunsellRenotation/real.dat}
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

National Television System Committee.
[Report and Reports of Panel No. 11, 11-A, 12-19,
with Some supplementary references cited in the Reports,
and the Petition for adoption of transmission standards
for color television before the Federal Communications Commission]
(1953)

Rheinboldt, Werner C. and 
John P. Menard.
\bold{Mechanized Conversion of Colorimetric Data to Munsell Renotations}.
Journal of the Optical Society of America.
Vol. 50, Issue 8, pp. 802-807.
August 1960.

Wikipedia.
Cubic Hermite spline.
\url{https://en.wikipedia.org/wiki/Cubic_Hermite_spline}

Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{http://www.99main.com/~centore/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{xyYtoMunsell}()}
}


\examples{
MunsellToxyY( '7.6P 8.9/2.2' )
}

\keyword{Munsell}

