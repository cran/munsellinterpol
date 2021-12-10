\name{MunsellSpecToColorlabFormat}
\alias{MunsellSpecToColorlabFormat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert Munsell Specification to Colorlab Format
}
\description{
Convert Munsell Specification to Colorlab Format
}
\usage{
MunsellSpecToColorlabFormat( MunsellSpec )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MunsellSpec}{
a numeric Nx3 matrix, or a vector that can be converted to such a matrix, by row.
Each row of the matrix contains an HVC vector.
H is automatically wrapped to the interval (0,100].\cr
\code{MunsellSpec} can also be a character vector with Munsell Notation;
which is converted to an Nx3 matrix using \code{\link{HVCfromMunsellName}()}.
}
}

\details{
Colorlab Munsell format uses 4 numbers.
\enumerate{
\item Hue Step, in the interval (0,10], or 0 for neutrals.
In Colorlab documentation it is called the \emph{hue shade}.
It is also the Hue Number H mod 10 (unless H is an exact multiple of 10).
\item Munsell Value, in the interval [0,10]
\item Munsell Chroma, non-negative
\item Hue Index, an integer from 1 to 10, or 0 for neutrals.
This index defines the \emph{principal hue}, see \bold{Details}.
}

}

\value{
an Nx4 matrix, with rows as described in \bold{Details}.
Value and Chroma are simply copied unchanged.
The complex part is conversion of Hue Number to Colorlab Hue Step and Hue Index.
If Chroma is 0, both the Hue Step and Hue Index are set to 0.
Invalid input values are converted to \code{NA}s.\cr
If the input is a character vector,
the rownames of the returned matrix are set to that vector.
}

\references{
Color Processing Toolbox.      
Colorlab 1.0.
\url{https://www.uv.es/vista/vistavalencia/software/colorlab.html}

Cleland, T. M.
\bold{A Practical description of the Munsell Color System with Suggestions for its Use}.
(1921)
}

\author{
Jose Gama and Glenn Davis
}

\note{
The Colorlab format is closer to the
\emph{Munsell Book of Color (1929)} than HVC.
In the book the hue circle is divided into 10
\emph{principal hues} - 5 simple and 5 compound.
The 10 hue labels are: \bold{R, YR, Y, GY G, BG, B, PB, P}
(simple are 1 letter and compound are 2 letters).
In Colorlab these labels are replaced by the Hue Index.
WARNING: In the Munsell System, see \emph{Cleland}, there is a different Hue Index -
\bold{R} is 1, \bold{YR} is 2, ..., \bold{P} is 10.
The Colorlab index has a different origin, and goes around the circle in a different
direction !
\cr
Each one of these major hues corresponds to an arc on the circle,
with a 10-point hue scale.
The midpoint of the arc has hue step = 5.
Fortunately this 10-point hue scale is exactly the same as the Colorlab Hue Step.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~




\seealso{
\code{\link{HVCfromMunsellName}()},
\code{\link{ColorlabFormatToMunsellSpec}()}
}

\examples{
MunsellSpecToColorlabFormat( c(100,5,5, 10,3,4, 90,4,3, 77,1,2, 66,2,0, 0,1,2 ) )
##                    HN V C HI
##  10.00RP 5.00/5.00 10 5 5  8
##  10.00R 3.00/4.00  10 3 4  7
##  10.00P 4.00/3.00  10 4 3  9
##  7.00PB 1.00/2.00   7 1 2 10
##  N 2.00/            0 2 0  0
##  10.00RP 1.00/2.00 10 1 2  8
}
 

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Colorlab }% use one of  RShowDoc("KEYWORDS")

