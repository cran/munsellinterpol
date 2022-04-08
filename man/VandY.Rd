\name{VandY}
\alias{VfromY}
\alias{YfromV}
\title{Convert Munsell Value V to Luminance Factor Y, and back again}

\description{Convert non-linear Munsell Value V to linear Luminance Factor Y, and back again
}
\usage{
YfromV( V, which='ASTM' )
VfromY( Y, which='ASTM' )}

\arguments{
\item{V}{a numeric vector with elements in the interval [0,10]}

\item{Y}{a numeric vector with elements in the interval [0,100]}

\item{which}{ one of the strings in the first column of this table
\tabular{ll}{
              \tab    reference \cr
\code{'ASTM'} \tab    ASTM D-1535 (2008) \cr
\code{'OSA'} \tab     Newhall, et. al. (1943) \cr
\code{'MgO'} \tab     Newhall, et. al. (1943) \cr
\code{'Munsell'} \tab Munsell, et. al. (1933)\cr
\code{'Priest'} \tab Priest, et. al. (1920)\cr
}
Partial matching is enabled, and it is case insensitive.
}
}



\details{
\code{'Priest'} is the earliest (1920) transfer function in this package.
It is implemented as:
\deqn{V = sqrt(Y)  ~~~~~~and~~~~~~   Y = V^2 }
One readily checks that when V=10, Y=100, and vice-versa.
This transfer function has been implemented in colorimeters,
using analog electric circuits.
It is used in Hunter Lab - the precursor of CIE Lab.

\code{'Munsell'} is the next (1933) transfer function,
and was proposed by Munsell's son (Alexander Ector Orr Munsell) and co-workers.   
It is implemented as:
\deqn{V = sqrt( 1.474*Y - 0.00474*Y^2 )}
\deqn{Y = 50 * ( (1474 - sqrt(1474^2 - 4*4740*V^2)) / 474 )}
One readily checks that when V=10, Y=100, and vice-versa.
The luminance factor Y is \emph{absolute},
AKA \emph{relative to the perfect reflecting diffuser}.

\code{'Priest'} and \code{'Munsell'} are included in this package
for historical interest only.

The remaining three define Y as a quintic polynomial in V.

The next one historically - \code{'MgO'} -  is implemented as:
\deqn{ Y = (((((8404*V - 210090)*V  + 2395100)*V  - 2311100)*V + 10000000)*V ) / 10000000}
One readily checks that when V=10, Y=102.568.
This Y is larger than 100, because the authors decided to make Y relative
to a clean surface of MgO, instead of the perfect reflecting diffuser.
In their words:

\itemize{
\item It should be noted that the
reflectances indicated are not absolute but relative to magnesium oxide;
whereas the maximum at value 10/ was formerly 100 percent, it is now 102.57.
Use of this relation facilitates results and
also avoids the somewhat dubious conversion to
absolute scale, by permitting Y determinations
with a MgO standard to be converted directly to Munsell value.
}

Nowadays, the perfect reflecting diffuser is preferred over MgO.
For users who would like to modify this quintic as little as possible,
with the perfect reflecting diffuser in mind
(going back to \code{'MUNSELL'}),
I offer \code{'OSA'}, which is given by this quintic of my own design
\deqn{Y = (((((8404*V - 210090)*V  + 2395100)*V  - 2311100)*V + 10000000)*V ) / 10256800}

ASTM had a similar modification in mind,
but did it a little differently by scaling each coefficient.
\code{'ASTM'} is given by this quintic:
\deqn{Y = ( ((((81939*V - 2048400)*V  + 23352000)*V - 22533000)*V + 119140000)*V  ) / 1.e8 }
One readily checks that when V=10, Y=100 exactly (for both \code{'OSA'} and \code{'ASTM'}).

The inverses - from Y to V - of all 3 quintics are implemented as 
3 \code{splinefun()}s at a large number (about 300) of points.
These inverses are both fast and accurate.
The round-trip Y  \rarrow  V  \rarrow  Y is accurate to 7 digits after the decimal.
The round-trip V  \rarrow  Y  \rarrow  V is accurate to 8 digits after the decimal.
}

\value{
a numeric vector the same length as the input
}

\note{
The quintic functions \code{'ASTM'} and \code{'OSA'} are very close.
They agree at the endpoints 0 and 10 exactly, and the largest difference is near
V=6.767 where they differ by about 0.0007.
}


\references{
Priest, I. G. Gibson, K. S. and McNicholas, H. J.
\bold{An Examination of the Munsell Color System.
I. Spectral and and Total Reflection and the Munsell Scale of Value}.
Technologic Papers of the Bureau of Standards, No. 167.
pp. 1-33.
Washington D.C.
1920.

Munsell, A. E. O.,
L. L. Sloan,
and I. H. Godlove.
\bold{Neutral Value Scales. I. Munsell Neutral Value Scale}.
Journal of the Optical Society of America.
Vol. 23. Issue 11. pp. 394-411.
November 1933.

Newhall, Sidney M.,
Dorothy Nickerson,
Deane B. Judd.
\bold{Final Report of the O.S.A. Subcommitte on the Spacing of the Munsell Colors}.
Journal of the Optical Society of America.
Vol. 33.  No. 7.  pp. 385-418.
July 1943.

ASTM D 1535-08.
Standard Practice for Specifying Color by the Munsell System.
2008
}

\author{Glenn Davis}

\seealso{
\code{\link{MunsellToxyY}()},
\code{\link{xyYtoMunsell}()}
}

\examples{
VfromY( c(0,50,100) )
##  [1]  0.00000  7.53772 10.00000
}

\keyword{datasets}

