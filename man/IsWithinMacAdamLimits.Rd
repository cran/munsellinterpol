\name{IsWithinMacAdamLimits}
\alias{IsWithinMacAdamLimits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test xyY Coordinates for being Inside the MacAdam Limits
}
\description{
Test xyY Coordinates for being Inside the MacAdam Limits for Illuminants C and D65
}
\usage{
IsWithinMacAdamLimits( xyY, Illuminant='C' )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{xyY}{a numeric Nx3 matrix with CIE xyY coordinates in the rows,
or a vector that can be converted to such a matrix, by row. 
The reference white is asssumed to satisfy Y=100.}

\item{Illuminant}{either 'C' or 'D65'. Partial matching is enabled, and is
case-insensitive.
}
}

\details{
The \emph{MacAdam Limits} is the boundary of the
\emph{optimal color solid} (also called the \emph{Rösch Farbkörper}),
in XYZ coordinates.
The optimal color solid  is convex and depends on the illuminant.
Points on the boundary of the solid are called \emph{optimal colors}.
This color solid has a smooth boundary (in the physical sense), except at black and white.
It is symmetric about the midpoint of the segment joining black and white
(the 50\% gray point).
It can be \emph{approximated} as the convex hull of a large number
of optimal colors.
For Illuminant D65 about 200 points (in xyY coordinates)
were computed and published by Wyszecki&Stiles (1982)..
For Illuminant C a similar set of 994 points was computed using package \bold{colorSpec};
for details please see the header of file \code{OptimalColorsForIlluminantC.txt}.
Tesselations of the convex hulls are stored privately inside this package,
and ready to go.
}

\note{
The function eventually calls \code{geometry::tsearchn()}.
}

\value{
A logical vector of length N.
A value is \code{TRUE} iff. the corresponding row in \code{xyY}
is inside the optimal color solid for the illuminant.
}

\seealso{
\code{\link[geometry]{delaunayn}()},
\code{\link[geometry]{tsearchn}()}
}


\references{
Rösch, S.
\bold{Darstellung der Farbenlehre für die Zwecke des Mineralogen}.
Fortschr. Mineral. Krist. Petrogr.
Vol. 13 No. 143.
1929.

MacAdam, David L.
\bold{Maximum Visual Efficiency of Colored Materials}.
Journal of the Optical Society of America.
Vol 25, No. 11.
pp. 361-367.
November 1935.

Günther Wyszecki and W. S. Stiles.
\bold{Color Science: Concepts and Methods, Quantitative Data and Formulae, Second Edition}.
John Wiley & Sons, 1982.
Table II(3.7). p. 778-779.
}


\source{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{http://centore.isletech.net/~centore/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}


\author{
Jose Gama and Glenn Davis
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
IsWithinMacAdamLimits( c(0.6,0.3,10, 0.6,0.3,20, 0.6,0.3,30, 0.6,0.3,40  ), 'C' )
##  [1]  TRUE  TRUE FALSE FALSE
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MacAdam }
