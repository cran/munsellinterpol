\name{CentroidsISCCNBS}
\alias{CentroidsISCCNBS}
\docType{data}
\title{Centroid Notations for the Revised ISCC-NBS Color-Name Blocks}

% \usage{CentroidsISCCNBS}

\description{\code{CentroidsISCCNBS} is a table with the centroids
of the revised ISCC-NBS Color-Name Blocks, taken from \emph{Kelly (1958)}}

\format{
This \code{data.frame} has 267 rows and these columns:
\describe{
\item{Number}{ISCC-NBS number (an integer from 1 to 267)}
\item{Name}{ISCC-NBS name}
\item{MunsellSpec}{Munsell specification of the centroid of the block a (character string)}
}
}

\details{
The earliest paper I am aware of is by \emph{Nickerson, et. al.} in 1941.
After the big Munsell renotation in 1943, the name blocks were revised in 1955.
\cr
When the central colors were recomputed in
\emph{Kelly (1958)}, they were called the "Central Colors",
though the text makes it clear that most are truly centroids,
which were computed from the centroid of an "elementary shape",
which is a "sector of a right cylindrical annulus".
For the "peripheral blocks" of high chroma, the centroids were
"estimated graphically by plotting the MacAdam limits".
\cr
In \emph{Kelly (1965)} these were called "centroid colors",
and that is the name we will use here.
}



\references{
Nickerson, Dorothy and Sidney M. Newhall.
\bold{Central Notations for ISCC-NBS Color names}.
J Opt. Soc. Am.
Vol 31  Iss. 9.
pp. 597-591.
1941.

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

Kelly, Kenneth Low.
\bold{Central Notations for the Revised ISCC-NBS Color-Name Blocks}.
Journal of Research of the National Bureau of Standards.
Research Paper 2911.
Vol. 61  No. 5.
pp. 427-431.
November 1958.

Kelly, Kenneth Low.
\bold{A Universal Color Language}.
Color Engineering.
Vol. 3  No. 2.
pp. 2-7.
March-April, 1965.
}

\section{Contributor}{Glenn Davis}

\examples{
print( CentroidsISCCNBS[ 1:5, ] )

##    Number          Name   MunsellSpec
##  1      1    vivid pink     1.5R 7/13
##  2      2   strong pink  1.5R 7.5/9.1
##  3      3     deep pink 1.9R 6.0/11.1
##  4      4    light pink  2.5R 8.6/5.2
##  5      5 moderate pink  2.5R 7.2/5.2
}

\keyword{datasets}

