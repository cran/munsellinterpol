\name{roundHVC}
\alias{roundHVC}

\title{Round a Given Munsell HVC to the Closest Chip in One or More Books of Chips}

\description{
The functions in this package compute Munsell HVC coordinates with high precision.
But when reporting the Munsell notation it is often desirable,
for historical consistency or other reasons,
to report a notation that is actually appears in one or more published books of
Munsell chips.
This rounding operation requires exact knowledge of such books,
which is taken from the data frame \code{\link{MunsellBooks}}.
}

\usage{roundHVC( HVC, books )}

\arguments{
\item{HVC}{a numeric Nx3 matrix with Munsell HVC values in the rows, or a vector that can be converted to such a matrix, by row.
\cr
\code{HVC} can also be a character N-vector with Munsell Notations,
which is converted to an Nx3 matrix using \code{\link{HVCfromMunsellName}()}. 
}

\item{books}{a comma-delimited string that specifies one or more books of chips.
The words in the string can be 
\code{'soil'}, \code{'rock'}, \code{'bead'}, \code{'plant'}, \code{'newstudent'},
\code{'glossy'}, and \code{'matte'}.
Matching is partial (initial substring only) and case-insensitive.
See \bold{References} and \bold{Examples}.
}
}

\value{
A \code{data.frame} with a N rows and these columns:

\item{HVC}{the given Nx3 input matrix HVC.}

\item{ISCC-NBS Name}{the ISCC-NBS Name (the color name) of the given input chip, see \code{\link{ColorBlockFromMunsell}()}. }

\item{MunsellRounded}{the Munsell notation that is rounded from HVC.
It is guaranteed to be in one of the books given by the argument \code{books}.}

\item{FergusonName}{the color name of \code{MunsellRounded}, as assigned by \bold{Ferguson}.
If the chip does not appear in one of the 5 specialized books in \code{\link{MunsellBooks}},
it is set to \code{NA}.}

The row names are set to the row names of HVC.
But if these are \code{NULL} they are set to the Munsell notation of the input HVC.
Finally, if these row names have duplicates, they are set to \code{1:N}.
}


\details{
First, the set of chips over all selected books is formed - a sub-frame of \code{\link{MunsellBooks}}.
Then the closest chip in this set to the given HVC is found.
For the color distance, the function
\code{\link{NickersonColorDifference}()} is used.
Unfortunately, the particular coefficients in this difference equation can lead to large regions
of Munsell color space with exact ties for the two closest chips.
These exact ties make the final result unpredictable and sometime unintuitive.
As a tie-breaker, a very small multiple (\code{1.e-6}) of plain
Eucliden distance is added to the Nickerson difference, which has proven effective.
}


\references{
\bold{Ferguson, Jonathan}.
\emph{Munsell notations and color names: Recommendations for archaeological practice}.
\bold{Journal of Field Archaeology}. 39:4. 327-335. (2014).
Online Supplement #1.
}

\author{Glenn Davis}


\seealso{
\code{\link{MunsellBooks}},
\code{\link{ColorBlockFromMunsell}()}
}

\examples{
## search the soil book for the best match for 3 chips
roundHVC( c("7.7YR 3.4/6.1", "2.6PB 6.1/4.5", "N 6.6/" ), books='soil' )

##               HVC.H HVC.V HVC.C ISCC-NBS Name MunsellRounded Ferguson Name
## 7.7YR 3.4/6.1  17.7   3.4   6.1  strong brown      7.5YR 4/6  Strong brown
## 2.6PB 6.1/4.5  72.6   6.1   4.5     pale blue        5PB 6/1   Bluish gray
## N 6.6/          0.0   6.6   0.0    light gray           N 7/    Light gray

## The middle chip has a poor match, since "purplish-blue" with Chroma=4.5
## is not well-matched by soil color chips.
## Search again, but this time in 5 books.

roundHVC( c("7.7YR 3.4/6.1", "2.6PB 6.1/4.5", "N 6.6/" ), books='s,b,r,p,n' )

##               HVC.H HVC.V HVC.C ISCC-NBS Name MunsellRounded Ferguson Name
## 7.7YR 3.4/6.1  17.7   3.4   6.1  strong brown      7.5YR 4/6  Strong brown
## 2.6PB 6.1/4.5  72.6   6.1   4.5     pale blue        5PB 6/4     Pale blue
## N 6.6/          0.0   6.6   0.0    light gray           N 7/    Light gray

## This time, for the middle chip, a better match was found.
## The other 2 chip are unchanged.
}


