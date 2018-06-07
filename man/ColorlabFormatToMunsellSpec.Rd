\name{ColorlabFormatToMunsellSpec}
\alias{ColorlabFormatToMunsellSpec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert Colorlab Munsell Format to Munsell HVC}
\description{Convert Colorlab Munsell Format to Munsell HVC}

\usage{
ColorlabFormatToMunsellSpec( HVCH )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HVCH}{
a numeric Nx4 matrix, or a vector that can be converted to such a matrix, by row.
Each row of the matrix contains an HVCH vector.\cr
}
}
 
\value{
an Nx3 matrix, with each row an HVC vector.
Value and Chroma are simply copied unchanged.
The complex part is conversion of Colorlab Hue Step and Hue Index
to Hue Number.
For neutrals, both Hue Step and Hue Index are ignored.
Invalid input values, such as a Hue Index that is not
an integer from 0 to 10 (except for neutrals), are converted to \code{NA}s.
The rownames of the input are copied to the output,
but if these are \code{NULL}, the rownames are set to the Munsell notations.
}

\details{
Colorlab Munsell format uses 4 numbers.
\enumerate{
\item Hue Step, in the interval (0,10], or 0 for neutrals
\item Munsell Value, in the interval (0,10]
\item Munsell Chroma, non-negative
\item Hue Index, an integer from 1 to 10, or 0 for neutrals
}
}


\references{
Color Processing Toolbox.      
Colorlab 1.0.
\url{http://www.uv.es/vista/vistavalencia/software/colorlab.html}
}
\author{
Jose Gama and Glenn Davis
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MunsellSpecToColorlabFormat}()}
}

\examples{
ColorlabFormatToMunsellSpec( c( 3.2,3,2,1,  2,5.1,0,0, 2,5.1,0.1,0 ) )
##                     H   V  C
##  3.20B 3.00/2.00 63.2 3.0  2
##  N 5.10/          0.0 5.1  0
##  <NA>              NA  NA NA
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Colorlab }% use one of  RShowDoc("KEYWORDS")

