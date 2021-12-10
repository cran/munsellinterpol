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
The reference white is asssumed to have Y=100.}

\item{Illuminant}{either \code{'C'} or \code{'D65'}. Partial matching is enabled, and is case-insensitive.}
}

\details{
The \emph{MacAdam Limit} is the boundary of the
\emph{optimal color solid} (also called the \emph{Rösch Farbkörper}),
in XYZ coordinates.
The optimal color solid is convex and depends on the illuminant.
Points on the boundary of the solid are called \emph{optimal colors}.
% This color solid has a smooth boundary (in the physical sense), except at black and white.
It is symmetric about the midpoint of the segment joining black and white
(the 50\% gray point).
It can be expressed as a \emph{zonohedron} - a convex polyhedron with a special form;
for details on zonohedra, see \cite{Centore}.

For each of the 2 illuminants, a zonohedron \eqn{Z} is pre-computed
(and stored in \code{sysdata.rda}).
The wavelengths used are 380 to 780 nm with 5nm step (81 wavelengths).
Each zonohedron has 81*80=6480 parallelogram faces, though some of them are coplanar.
\eqn{Z} is expressed as the intersection of 6480 halfspaces.
The plane equation of each parallelogram is pre-computed,
but redundant ones are not removed (in this version).

For testing a query point xyY, a pseudo-distance metric \eqn{\delta} is used.
Let the zonohedron \eqn{Z} be the intersection of the halfspaces
\langle \eqn{h_i,x} \rangle \eqn{\le b_i ~~ i=1,...,n}, 
where each \eqn{h_i} is a unit vector.
The point xyY is converted to XYZ, and \eqn{\delta}(XYZ) is computed as:
\eqn{\delta}(XYZ) := max( \langle \eqn{h_i},XYZ \rangle - \eqn{b_i} )
where the maximum is taken over all \eqn{i=1,...,n}.
This calculation can be optimized;
because the zonohedron is centrally symmetric, only half of the planes actually have to be stored,
and this cuts the memory and processing time in half.
It is clear that XYZ is within the zonohedron iff \eqn{\delta}(XYZ) \eqn{\le} 0,
and that XYZ is on the boundary iff \eqn{\delta}(XYZ)=0.
This pseudo-distance is part of the returned \code{data.frame}.

An interesting fact is that if \eqn{\delta}(XYZ)>0,
then \eqn{\delta}(XYZ) \eqn{\le} dist(XYZ,\eqn{Z}),
with equality iff the segment from \eqn{XYZ} to the
point \eqn{z} on the boundary of \eqn{Z} closest to XYZ is
normal to one of the faces of \eqn{Z} that contains \eqn{z}.
This is why we call \eqn{\delta} a \emph{pseudo-distance}.
Another interesting fact is that if \eqn{\delta}(XYZ) \eqn{\le} 0,
then \eqn{\delta}(XYZ) = -min( \eqn{\Psi_Z(u)} - \langle \eqn{u},XYZ \rangle ),
where the minimum is taken over all unit vectors \eqn{u} and
where \eqn{\Psi_Z} is the \emph{support function} of \eqn{Z}.
}


\value{
A \code{data.frame} with N rows and these columns:
\item{within}{a logical which is \code{TRUE} iff the corresponding row in \code{xyY}
is inside the optimal color solid for the illuminant.
If a point is exactly on the boundary (unlikely), \code{within=TRUE}.
Explicitly, \code{within = (delta<=0)}.}
\item{delta}{the pseudo-distance \eqn{\delta} discussed in \bold{Details}}

The row names of the output value are set equal to the row names of \code{xyY}.
}



\references{
Centore, Paul.
\emph{A zonohedral approach to optimal colours}.
\bold{Color Research & Application}.
Vol. 38.
No. 2.
pp. 110-119.
April 2013.

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

Wikipedia.
\bold{Support Function}.
\url{https://en.wikipedia.org/wiki/Support_function}
}


\author{
Glenn Davis and Jose Gama
}



\examples{
IsWithinMacAdamLimits( c(0.6,0.3,10, 0.6,0.3,20, 0.6,0.3,30, 0.6,0.3,40  ), 'C' )

##    within  delta
##  1   TRUE -1.941841
##  2   TRUE -1.332442
##  3  FALSE  3.513491
##  4  FALSE 12.826172
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MacAdam }
