\name{ColorBlockFromMunsell}
\alias{ColorBlockFromMunsell}

\title{Get ISCC-NBS Number and ISCC-NBS Name from Munsell Hue, Value, and Chroma}

\description{Get ISCC-NBS Number and ISCC-NBS Name from Munsell Hue, Value, and Chroma.
}

\usage{ColorBlockFromMunsell( MunsellSpec ) }

\arguments{
\item{MunsellSpec}{ a numeric Nx3 matrix or a vector that can be converted to
such a matrix.
Each row has Munsell HVC, where H is Hue Number,
and V and C are the standard Munsell Value and Chroma.
The Hue is automatically wrapped to the interval (0,100].
\code{MunsellSpec} can also be a character N-vector with standard Munsell notation; it is converted to an Nx3 matrix.}
}

\details{
The ISCC-NBS System is a partition of Munsell Color Solid into 267 color blocks.
Each block is a disjoint union of \emph{elementary blocks},
where an elementary block is defined
by its minimum and maximum limits in Hue, Value, and Chroma.
Some blocks are non-convex.
The \emph{peripheral blocks}, of which there are 120,
have arbitrary large chroma and are considered semi-infinite for this function;
there is no consideration of the MacAdam limits.
For each query vector HVC,
the function searches a private \code{data.frame} with 932 elementary blocks,
for the one elementary block that contains it.
}

\section{History}{
The \bold{Munsell Book of Color} was published in 1929.
The first ISCC-NBS partition, in 1939, had 319 blocks and names (including 5 neutrals).
There were no block numbers.
The aimpoints of the Munsell samples were thoroughly revised in 1943.
The ISCC-NBS partition was revised in 1955,
and this is the version used here.
}

\section{Future Work}{
It might be useful to compute the distance from the query point
to the boundary of the containing color block.
}

\value{
a \code{data.frame} with N rows and these columns:
\item{HVC}{the input Nx3 matrix, or such a matrix converted from Munsell notation}
\item{Number}{the corresponding ISCC-NBS color number - an integer from 1 to 267}
\item{Name}{the corresponding ISCC-NBS color name - a character string}
\item{Centroid}{the centroid of the block in Munsell Notation - a character string; see \code{\link{CentroidsISCCNBS}}}
The \code{rownames} are set to the input \code{MunsellSpec}.
}

\references{
Munsell Color Company, A.H. Munsell, and A.E.O. Munsell.
\bold{Munsell book of color: defining, explaining, and illustrating the fundamental characteristics of color}.
1929.

Judd, Deane B. and 
Kenneth L. Kelly.
\bold{Method of Designating Colors}.
Journal of Research of the National Bureau of Standards.
Research Paper 1239.
Volume 23  Issue 3.
pp. 355-385.
September 1939.

Newhall, Sidney M.,
Dorothy Nickerson,
Deane B. Judd.
\bold{Final Report of the O.S.A. Subcommitte on the Spacing of the Munsell Colors}.
Journal of the Optical Society of America.
Vol. 33.  No. 7.  pp. 385-418.
July 1943.

Kelly, Kenneth L. and 
Deane B. Judd
\bold{The ISCC-NBS Method of Designating Colors and a Dictionary of Color Names}.
National Bureau of Standards Circular 553.
Washington DC: US Government Printing Office.
November 1, 1955.
}

\examples{
ColorBlockFromMunsell( c( "3R 8/3", "7.4YR 3/4" ) )

##            HVC.H HVC.V HVC.C Number           Name      Centroid
##  3R 8/3      3.0   8.0   3.0      4     light pink  2.5R 8.6/5.2
##  7.4YR 3/4  17.4   3.0   4.0     58 moderate brown 5.5YR 3.5/3.9
}


\seealso{
\code{\link{CentroidsISCCNBS}}
}

\author{Glenn Davis}
\keyword{ISCC-NBS}

