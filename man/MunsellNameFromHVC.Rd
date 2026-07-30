\name{MunsellNameFromHVC}
\alias{MunsellNameFromHVC}
\alias{HueStringFromNumber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert Munsell Numeric to Munsell String Notation
}
\description{
Convert Munsell Numeric to Munsell String Notation
}

\usage{
MunsellNameFromHVC( HVC, format='g', digits=2, ctol=0 )

HueStringFromNumber( Hue, format='g', digits=2 )
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{HVC}{a numeric Nx3 matrix or a vector that can be converted to
such a matrix.
Each row has an HVC vector, where H is Hue Number,
and V and C are the standard Munsell Value and Chroma.
The Hue is automatically wrapped to the interval (0,100].}

\item{Hue}{a numeric vector of Hue Numbers,
which are automatically wrapped to the interval (0,100].}

\item{format}{determines the meaning of the argument \code{digits},
and usually equal to \code{'g'} or \code{'f'} though other options are available,
see \code{\link[base]{formatC}()} for details.
When \code{format='f'} trailing 0s might be displayed.
\cr
For \code{MunsellNameFromHVC()}, \code{format} may have length 2,
and in that case the first one is applied to Hue,
and the second one is applied to Value and Chroma.
It may also have length 3, and in that case the 3 values are applied
to Hue, Value, and Chroma in that order.
}

\item{digits}{
when \code{format='g'} the number of significant digits,
and when \code{format='f'} the number of digits displayed after the decimal point.
Both \code{format} and \code{digits} are passed as arguments to
\code{\link[base]{formatC}()}.
\cr
For \code{MunsellNameFromHVC()}, \code{digits} may have length 2,
and in that case the first one is applied to Hue,
and the second one is applied to Value and Chroma.
It may also have length 3, and in that case the 3 values are applied
to Hue, Value, and Chroma in that order.
}

\item{ctol}{
a small tolerance that is used to "snap" a very small Chroma to 0,
which is then converted to an exact Munsell Neutral, see \bold{Examples}.
The default value, \code{ctol=0} has no effect on the Chroma. 
}
}


\value{
Both functions return  a character vector of length N.
\code{MunsellNameFromHVC()} returns the full notation.
\code{HueStringFromNumber()} returns just initial the hue part;
which is useful for labeling plots.
}



\section{Condition Signaling}{
For an overview see:
\href{../doc/munsellinterpol-guide.html#appendix-e---logging-and-condition-signaling}{\bold{Appendix E - Logging and Condition Signaling}} in the User Guide.

In case of an \code{ERROR} event, an error condition is signaled with the additional
S3 class \code{"invalid_argument"}.
The invalid argument is added to the condition object.

In case of an \code{WARN} event, a warning condition is signaled with the additional
S3 class \code{"incomplete"}.
The condition object also has the submatrix of the invalid rows of \code{HVC},
or the subvector of invalid elements of \code{Hue},
as item \code{invalid}, and their indexes as item \code{indexes}.
}

\note{
If \code{format='f'} for Chroma, then Chroma is first \emph{rounded} to the given \code{digits}.
Chromas close to 0 may then become 0 and be converted to an exact Munsell Neutral,
see \bold{Examples}.
\cr
The \code{width} argument of \code{\link[base]{formatC}()} is always set to 1,
to suppress leading spaces.
}

\references{
ASTM D 1535-97.
Standard Practice for Specifying Color by the Munsell System.
1997
}

\author{Glenn Davis and Hugo Rodrigues}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[base]{formatC}()},
\code{\link{HVCfromMunsellName}()},
\code{\link{HueNumberFromString}()}
}

\examples{
white.D65 = c(95.047,100.000,108.883)

HVC = XYZtoMunsell( c( 0.01*white.D65, 0.1*white.D65), white=white.D65 )
HVC
##                            H         V            C
##  6.6GY 0.86/2.8e-15 36.62033 0.8634072 2.753467e-15
##  N 3.7/              0.00000 3.7206664 0.000000e+00

# note that the Chroma in the 2nd row is 0 as it should be,
# but the Chroma in the first row is a tiny bit off.
# This is annoying when converted to Munsell notation, using defaults.

MunsellNameFromHVC( HVC )
##  [1] "6.6GY 0.86/2.8e-15" "N 3.7/"   

# Convert again with a positive ctol.

MunsellNameFromHVC( HVC, ctol=1.e-4 )
##  [1] "N 0.86/" "N 3.7/"      # now *both* are converted to Munsell Neutrals


# setting format='f' can achieve a similar goal

MunsellNameFromHVC( c(39,5.1,7.3,  14.1234,5.1234,0.0003 )  )
##  [1] "9GY 5.1/7.3"      "4.1YR 5.1/0.0003"

MunsellNameFromHVC( c(39,5.1,7.34,  14.1234,5.1234,0.003 ), format='f' )
##  [1] "9.00GY 5.10/7.34" "N 5.12/"     


#  On the other hand, one may *want* to see more precision for V and C, but not H
MunsellNameFromHVC( c(39,5.1,7.3,  14.1234,5.1234,0.0003 ), digits=c(2,6)  )
##  [1] "9GY 5.1/7.3"         "4.1YR 5.1234/0.0003"


# make nice labels for a plot with all the Hues
HueStringFromNumber( seq( 2.5, 100, by=2.5 ) )  
##   [1] "2.5R"  "5R"    "7.5R"  "10R"   "2.5YR" "5YR"   "7.5YR" "10YR"  "2.5Y" 
##  [10] "5Y"    "7.5Y"  "10Y"   "2.5GY" "5GY"   "7.5GY" "10GY"  "2.5G"  "5G"   
##  [19] "7.5G"  "10G"   "2.5BG" "5BG"   "7.5BG" "10BG"  "2.5B"  "5B"    "7.5B" 
##  [28] "10B"   "2.5PB" "5PB"   "7.5PB" "10PB"  "2.5P"  "5P"    "7.5P"  "10P"  
##  [37] "2.5RP" "5RP"   "7.5RP" "10RP" 
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MunsellName }% use one of  RShowDoc("KEYWORDS")

