\name{roundHVC}
\alias{roundHVC}

\title{Round a Given Munsell HVC to the Closest Sample in One or More Books of Samples}

\description{
The functions in this package compute Munsell HVC coordinates with high precision.
When reporting the Munsell notation it is often desirable,
for historical consistency or other reasons,
to report a notation that is actually appears in one or more published books of
Munsell samples.
This rounding operation requires exact knowledge of such books.
We use a paper by \bold{Ferguson}, which gives the contents of books on these fields:
soil, rocks, beads, and plants.  
It also includes \emph{The New Munsell Student Color Set}, see \bold{References}.
}

\usage{roundHVC( HVC, books )}

\arguments{
\item{HVC}{a numeric Nx3 matrix with Munsell HVC values in the rows, or a vector that can be converted to such a matrix, by row.
\cr
\code{HVC} can also be a character N-vector with Munsell Notations,
which is converted to an Nx3 matrix using \code{\link{HVCfromMunsellName}()}. 
}

\item{books}{a comma-delimted string that specifies one or more books of samples.
The words in the string can be 
\code{'soil'}, \code{'rock'}, \code{'bead'}, \code{'plant'}, or \code{'newstudent'}.
Matching is partial (initial substring only) and case-insensitive.
See \bold{References} and \bold{Examples}.
}
}

\value{
A \code{data.frame} with a N rows and these columns:

\item{HVC}{the original Nx3 input matrix HVC.}

\item{ISCC-NBS Name}{the ISCC-NBS Name (the color name) of the given input sample. }

\item{MunsellRounded}{the rounded Munsell notation that is rounded from HVC.
It is guaranteed to be in one of the books given by the argument \code{books}.}

\item{Ferguson Name}{the color name of \code{MunsellRounded}, as assigned by \bold{Ferguson}}

The row names are set to the row names of HVC.
But if these are \code{NULL} they are set to the Munsell notation of the input HVC.
Finally, if these row names have duplicates, they are set to \code{1:N}.
}


\details{
First, the set of samples over all selected books is formed.
Then the closest sample in this set to the given notation is found.
For the color distance, the function
\code{\link{NickersonColorDifference}()} is used.
Unfortunately, the particular coefficients in this difference equation can lead to large regions
of Munsell color space with exact ties for the two closest samples.
These exact ties make the final result unpredictable and sometime unintuitive.
As a tie-breaker, a very small multiple (\code{1.e-6}) of plain
Eucliden distance is added to the Nickerson difference, which has proven effective.
}


\source{
The data on books of samples is provide in
\code{.../munsellinterpol/inst/extdata/Supplement1_3.6.2024.csv}.
This was converted from the original file \code{yjfa_a_11710256_sm0001.docx}
which is \bold{Online Supplement #1} from \bold{Ferguson}.

The sample counts in the 5 books are:

\tabular{lc}{
\bold{Book} \tab \bold{Sample Count} \cr
Soil   \tab 437  \cr
Rock    \tab 115 \cr
Bead  \tab   176 \cr
Plant   \tab  332 \cr
New Munsell Student \tab 246
}
The union of all 5 books has 869 samples.
}

\section{Acknowledgements}{
The conversion of the book data to 
\code{Supplement1_3.6.2024.csv} was done by Willie Ondricek,
who provided all the motivation for this function.
He also kindly provided photographs of the soil book, and much other material.
}


\references{
\bold{Ferguson, Jonathan}.
\emph{Munsell notations and color names: Recommendations for archaeological practice}.
\bold{Journal of Field Archaeology}. 39:4. 327-335. (2014).
Online Supplement #1.

\bold{Long, J. T.}
\emph{The New Munsell Student Color Set}.
New York: Fairchild Books. 2011.

\bold{Munsell Color}. 
\emph{Munsell Soil-Color Charts}.
Grand Rapids, MI:  2009.

\bold{Munsell Color}. 
\emph{Geological Rock-Color Chart}.
GrandRapids, MI: 2009.

\bold{Munsell Color}.
\emph{Munsell Bead Color Book}.
Grand Rapids, MI: 2012.

\bold{Munsell Color}.
\emph{Munsell Color Charts for Plant Tissues}.
New Windsor, NY:  1977.
}

\author{Glenn Davis}



\examples{
## search the soil book for the best match for 3 samples
roundHVC( c("7.7YR 3.4/6.1", "2.6PB 6.1/4.5", "N 6.6/" ), books='soil' )

##               HVC.H HVC.V HVC.C ISCC-NBS Name MunsellRounded Ferguson Name
## 7.7YR 3.4/6.1  17.7   3.4   6.1  strong brown      7.5YR 4/6  Strong brown
## 2.6PB 6.1/4.5  72.6   6.1   4.5     pale blue        5PB 6/1   Bluish gray
## N 6.6/          0.0   6.6   0.0    light gray           N 7/    Light gray

## The middle sample has a poor match, since "purplish-blue" with Chroma=4.5
## is not well-matched by soil color samples.
## Search again, but this time in all 5 books.

roundHVC( c("7.7YR 3.4/6.1", "2.6PB 6.1/4.5", "N 6.6/" ), books='s,b,r,p,n' )

##               HVC.H HVC.V HVC.C ISCC-NBS Name MunsellRounded Ferguson Name
## 7.7YR 3.4/6.1  17.7   3.4   6.1  strong brown      7.5YR 4/6  Strong brown
## 2.6PB 6.1/4.5  72.6   6.1   4.5     pale blue        5PB 6/4     Pale blue
## N 6.6/          0.0   6.6   0.0    light gray           N 7/    Light gray

## This time, for the middle sample, a better match was found.
## The other 2 sample are unchanged.
}


