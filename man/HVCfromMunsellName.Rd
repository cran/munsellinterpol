\name{HVCfromMunsellName}
\alias{HVCfromMunsellName}
\alias{MunsellHVC}
\alias{HueNumberFromString}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert Munsell Notation to numerical Munsell HVC}

\description{Convert Munsell Notation to numerical Munsell HVC}

\usage{
HVCfromMunsellName( MunsellName )
MunsellHVC( MunsellName )
HueNumberFromString( HueString )
}

\arguments{
\item{MunsellName}{
a character vector of length N > 0,
where each string should be a valid Munsell notation, e.g.
\code{'2.3P 5/2.3',  '9.2YR 3/6',  'N 2.3/'}, and \code{'N 4/0'}. 
Whitespace is optional and ignored.
It is OK for a neutral to end in either \code{'/'} or \code{'/0'}.}

\item{HueString}{
a character vector of length N > 0,
where each string should be the initial hue part of a Munsell notation,
e.g. '4.5GY', '2.5R', '10.3B', etc.
Whitespace is optional and ignored.
Neutrals, denoted by \code{'N'}, are invalid because the hue is undefined.}
}


\value{
\code{HVCfromMunsellName()}
returns a numeric Nx3 matrix with HVC in the rows.
For neutral colors, both H and C are set to 0.
If a string cannot be parsed, the entire row is set to \code{NA}s.
The rownames are set to \code{MunsellName}.
\cr
\code{MunsellHVC()}
returns a character Nx3 matrix with HVC in the rows,
and is there for backward compatibility with older versions of the package.
For neutral colors, H is set to \code{'N'} and C is set to \code{'0'}.
\cr
\code{HueNumberFromString()} returns the hue number H (in (0,100]).
If the string cannot be parsed, or the color is neutral, the output is set to \code{NA}.
\cr
For all functions the Hue Number is wrapped to (0,100].
}

\note{
Ever since the \emph{Munsell Book of Color (1929)},
the Munsell hue circle has been divided into 10 \emph{principal hues} or arcs.
And each principal hue has been assigned a 10-point scale,
with 5 at the midpoint of the arc.
Moreover, the hue "origin" has been at '10RP'.
So a 100-point scale (with no letters) for the entire hue circle is obvious and trivial
to construct,
but I have been unable to determine the first explicit mention of such a scale.
The earliest I have have found is from \emph{Nimeroff (1968)},
Figures 20 and 21 on page 27.

\if{html}{\out{<div style="text-align: center">}}

\if{html}{\figure{hue-circle.jpg}{options: alt="Hue Circle"}}

\if{html}{\out{</div>}}

\if{latex}{\figure{hue-circle.jpg}}

There is a reference to \emph{ASTM D 1535} in the \bold{References} of \emph{Nimeroff},
but it is not dated, and the 2 figures are not attributed to it.
There was an \emph{ASTM D 1535} in 1968 but I have not been able to locate it;
it is possible that the 100-point scale first appeared in \emph{ASTM D-1535 (1968)},
or even earlier in \emph{ASTM D 1525-58T (1958)}.
\cr
Interestingly, in the \emph{Atlas of the Munsell Color System (1915)} there were only
5 principal hues, and each arc was assigned a 10-point scale.
If the entire hue circle of 1915 were assigned a scale, it would have been a 50-point scale.


}



\references{
Nimeroff, I.
\bold{Colorimetry}.
National Bureau of Standards Monograph 104.
January 1968. 35 cents.

ASTM D 1535-80.
Standard Practice for Specifying Color by the Munsell System.
1980.

Munsell Book of Color: defining, explaining, and illustrating the fundamental characteristics of color.
Munsell Color Co.
1929.

Atlas of the Munsell Color System.
Malden, Mass., Wadsworth, Howland & Co., inc., Printers.
1915.
}

\author{
Glenn Davis
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MunsellNameFromHVC}()},
\code{\link{HueStringFromNumber}()}
}

\examples{
HVCfromMunsellName( c( "4.2P 2.9/3.8", "N 2.3/", "N 8.9/0" ) )
##                  H   V   C
##  4.2P 2.9/3.8 84.2 2.9 3.8
##  N 2.3/        0.0 2.3 0.0
##  N 8.9/0       0.0 8.9 0.0

HueNumberFromString( c('4B','4.6GY','10RP','N') )
##  [1]  64.0  34.6 100.0    NA
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MunsellName}% __ONLY ONE__ keyword per line
