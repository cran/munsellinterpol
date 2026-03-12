\name{MunsellBooks}

\alias{MunsellBooks}
\docType{data}
\title{The Munsell Chips present in Specialized and Full Munsell Books}

\description{
In addition to the standard glossy and matte books,
Munsell publishes books for specialized fields.
For 5 of these specialized books, we use a paper by \bold{Ferguson}, which gives the
chip contents of books for these fields:
\bold{soil}, \bold{rocks}, \bold{beads}, and \bold{plants}.
It also includes \bold{newstudent} \emph{Color Set}, see \bold{References}.
\cr
For the full \bold{glossy} and \bold{matte} books, see \bold{Source}.
}


\section{Acknowledgements}{
The conversion of the 5 specialized books data to 
\code{./inst/extdata/Supplement1_3.6.2024.csv} was done by Willie Ondricek,
who provided all the motivation for this function.
He also kindly provided photographs of the soil book, and much other material.
}


\format{
A data.frame with 1753 rows and these columns:

\tabular{ll}{
{\code{HVC}} \tab 
a matrix with 3 columns: Hue, Value, and Chroma.
For the chromatic chips, Hue is a positive multiple of 2.5,
except for 3 chips in the bead book which are odd multiples of 1.25.
For the achromatic neutral chips, Hue is really undefined, but here it is set to 0.
Value is always a positive multiple of 0.25, and between 0.5 and 9.5.
Chroma is a non-negative integer between 0 and 16, except 13 and 15 are omitted.
\cr
\cr
{\code{FergusonName}} \tab
the color name of the chip, as assigned by \bold{Ferguson}.
If a chip does not appear in one of the 5 specialized books, it is set to \code{NA}.
\cr
\cr
{\code{bead}}         \tab
a logical value. If \code{TRUE} then the chip appears in the \bold{bead} book
\cr
\cr
{\code{newstudent}}         \tab
a logical value. If \code{TRUE} then the chip appears in the \bold{newstudent} book
\cr
\cr
{\code{plant}}         \tab
a logical value. If \code{TRUE} then the chip appears in the \bold{plant} book
\cr
\cr
{\code{rock}}         \tab
a logical value. If \code{TRUE} then the chip appears in the \bold{rock} book
\cr
\cr
{\code{soil}}         \tab
a logical value. If \code{TRUE} then the chip appears in the \bold{soil} book
\cr
\cr
{\code{glossy}}         \tab
a logical value. If \code{TRUE} then the chip appears in the \bold{glossy} book
\cr
\cr
{\code{matte}}         \tab
a logical value. If \code{TRUE} then the chip appears in the \bold{matte} book
}
The rows are sorted first by Hue, then Value, and then Chroma, all in increasing order.
Thus the achromatic neutral chips come first, followed by the chromatics.
}


\details{
The chip counts in the 7 books are:

\tabular{lccc}{
\bold{Book}         \tab  \bold{| Neutral Chips |} \tab  \bold{Chromatic Chips |} \tab \bold{Total Chips |} \cr
Soil                \tab           10          \tab          427            \tab         437  \cr
Rock                \tab            9          \tab          106            \tab         115  \cr
Bead                \tab            5          \tab          171            \tab         176  \cr
Plant               \tab            0          \tab          332            \tab         332  \cr
New Munsell Student \tab            8          \tab          238            \tab         246  \cr
Glossy              \tab           37          \tab         1485            \tab        1522  \cr
Matte               \tab           32          \tab         1269            \tab        1301  \cr
}
The union of all 7 books has 1753 chips.
There are only 3 inclusions among the books:
Rock \eqn{\subset} Glossy,  Plant \eqn{\subset} Glossy,  and New_Munsell_Student \eqn{\subset} Matte.
}


\source{
The data on the first 5 books of chips in the following table was converted from the
file
\cr
\code{yjfa_a_11710256_sm0001.docx},
which is \bold{Online Supplement #1} from \bold{Ferguson}.

For the Matte book, the 1269 chromatic chips
are taken from:
\cr
\code{https://sites.uef.fi/spectral/databases-software/munsell-colors-matt-spectrofotometer-measured/}.
\cr
For the Glossy book, the 1485 chromatic chips
are taken from:
\cr
\code{https://www.munsellcolorscienceforpainters.com/MunsellResources/MunsellResources.html}
\cr
\code{/SpectralReflectancesOf2007MunsellBookOfColorGlossy.txt}.

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


\seealso{
\code{\link{roundHVC}()}
}


\examples{
str(MunsellBooks)

## 'data.frame':   1753 obs. of  9 variables:
## $ HVC         : num [1:1753, 1:3] 0 0 0 0 0 0 0 0 0 0 ...
##  ..- attr(*, "dimnames")=List of 2
##  .. ..$ : NULL
##  .. ..$ : chr [1:3] "H" "V" "C"
## $ FergusonName: chr  NA NA "Black" NA ...
## $ bead        : logi  FALSE FALSE TRUE FALSE FALSE FALSE ...
## $ newstudent  : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
## $ plant       : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
## $ rock        : logi  FALSE FALSE TRUE FALSE FALSE FALSE ...
## $ soil        : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
## $ glossy      : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
## $ matte       : logi  FALSE FALSE FALSE FALSE FALSE TRUE ...
}

\keyword{datasets}
