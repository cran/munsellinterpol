\name{adaption}
\alias{adaption}
\alias{adaptxyY}
\alias{adaptXYZ}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Chromatic Adaption
}
\description{
Adapt XYZ or xyY from a source viewing enviroment with a given illuminant,
to a destination viewing environment with a different illuminant
}
\usage{
adaptXYZ( XYZ.src, white.src, white.dest, method="bradford" )
adaptxyY( xyY.src, white.src, white.dest, method="bradford" )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{XYZ.src}{
an Nx3 matrix, or a vector that can be converted to such a matrix, by rows.
Each row has an XYZ in the source viewing environment.}

\item{xyY.src}{
an Nx3 matrix, or a vector that can be converted to such a matrix, by rows.
Each row has an xyY in the source viewing environment.}

\item{white.src}{the XYZ or xyY of the illuminant in the source viewing environment}

\item{white.dest}{the XYZ or xyY of the illuminant in the destination viewing environment}

\item{method}{
the method used for the chromatic adaption. Available methods are:
\code{"Bradford"}, \code{"VonKries"}, \code{"MCAT02"}, and \code{"scaling"}.
See \bold{References}.
Partial matching is enabled, and matching is case-insensitive.
}
}



\details{
\code{adaptXYZ()} is the more fundamental of the two.
\code{adaptxyY()} simply calls \code{adaptXYZ()} using
\code{\link{xyY2XYZ}()} and \code{\link{XYZ2xyY}()} to do conversions.
However, \code{adaptxyY()} does do an addition check:
if the chromaticity of a row in \code{xyY.src} exactly matches that of \code{white.src},
then the returned chromaticity is set to be exactly the same as that in \code{white.dest}.  Because of numerical truncation, it might not be.
\cr
In \code{adaptXYZ()} the \code{white} coordinates are \code{XYZ}.
In \code{adaptxyY()} the \code{white} coordinates are \code{xyY}.
}

\value{
\code{adaptXYZ()} returns an Nx3 matrix with adapted XYZ in each row.
\code{adaptxyY()} returns an Nx3 matrix with adapted xyY in each row.
}

\note{
All these adaption methods are linear.
For \code{MCAT02} only the simplified linear variant is used.
Linear methods (transforms) have these desirable properties:
\itemize{
\item symmetry: If A and B are viewing enironments, then adapting XYZ from A to B
is the inverse of adapting XYZ from B to A.
\item  commutative triangle:
If there are 3 viewing environments, say A, B, and C,
then adapting from A to B, and then from B to C,
is the same as adapting from A to C.
Equivalently (assuming symmetry), a round-trip around the triangle is the identity.
The property reduces to properties of matrix multiplication.
The property makes it possible to use the intermediate
\emph{Profile Connection Space} (PCS) in v. 4 ICC color profiles.
In the above-mentioned diagram, B represents the PCS, and A and C represent
viewing environments for devices.
The profile tag \code{'chad'} is the \emph{chromatic adaptation matrix}.
}
According to \emph{Hunt} p. 591, the non-linear CAT97 adaption transforms
is not symmetric.
\cr
I do not know whether any non-linear adaption transforms violate the
second property, but guess that some do.
}

\section{Analogy}{
Chromatic adaption can be viewed as an Aristotelian analogy of proportions.
A general analogy of this type is usually written \emph{A:B = C:D} and read as
"\emph{A} is to \emph{B} as \emph{C} is to \emph{D}".
In our case the expression \emph{A:I} is read as
"the appearance of \emph{A} in a viewing environment with illuminant \emph{I}",
or more simply
"the appearance of \emph{A} under illuminant \emph{I}".
It is better to think of A not as an object color, but as a self-luminous color.
The analogy \eqn{A:I = B:J} can be read as
"the appearance of \eqn{A} under illuminant \eqn{I}
is the same as the appearance of \eqn{B} under illuminant \eqn{J}".
\cr
Solving this problem of chromatic adaption is solving
this analogy for \emph{X}:  \emph{A:I = X:J},
where \emph{I} is the the source illuminant
and \emph{J} is the destination illuminant.
}


\references{
Hunt, R. W. G.
\bold{The Reproduction of Colour}.  6th Edition.
John Wiley & Sons.
2004.

International Color Consortium.
ICC.1:2001-04.
File Format for Color Profiles.
2001.

Lindbloom, Bruce.
Chromatic Adaptation.
\url{http://brucelindbloom.com/Eqn_ChromAdapt.html}

Wikipedia.
CIECAM02.
\url{https://en.wikipedia.org/wiki/CIECAM02}
}

\author{
Glenn Davis
}



%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MunsellToRGB}()}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{adaption}% use one of  RShowDoc("KEYWORDS")

