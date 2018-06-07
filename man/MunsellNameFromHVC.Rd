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
MunsellNameFromHVC( HVC, format='g', digits=2 )

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
When \code{format='f'} trailing 0s might be displayed.}

\item{digits}{
when \code{format='g'} the number of significant digits,
and when \code{format='f'} the number of digits displayed after the decimal point.
Both \code{format} and \code{digits} are passed as arguments to
\code{\link[base]{formatC}()}.}
}


\value{
Both functions return  a character vector of length N.
\code{MunsellNameFromHVC()} returns the full notation.
\code{HueStringFromNumber()} returns just initial the hue part;
which is useful for labeling plots.
}

\note{
If \code{format='f'}, then Chroma is first \emph{rounded} to to the given \code{digits}.
Chromas close to 0 may then become 0 and be displayed as a neutral, see \bold{Examples}.
\cr
The \code{width} argument of \code{\link[base]{formatC}()} is always set to 1,
to suppress leading spaces.
}

\references{
ASTM D 1535-97.
Standard Practice for Specifying Color by the Munsell System.
1997
}

\author{Glenn Davis}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[base]{formatC}()},
\code{\link{HVCfromMunsellName}()},
\code{\link{HueNumberFromString}()}
}

\examples{
MunsellNameFromHVC( c(39,5.1,7.3,  0,5.1234,0.003 )  )
##  [1] "9GY 5.1/7.3"    "10RP 5.1/0.003"

MunsellNameFromHVC( c(39,5.1,7.34,  0,5.1234,0.003 ), format='f' )
##  [1] "9.00GY 5.10/7.34" "N 5.10/"

HueStringFromNumber( seq( 2.5, 100, by=2.5 ) )   # make nice labels for a plot
##   [1] "2.5R"  "5R"    "7.5R"  "10R"   "2.5YR" "5YR"   "7.5YR" "10YR"  "2.5Y" 
##  [10] "5Y"    "7.5Y"  "10Y"   "2.5GY" "5GY"   "7.5GY" "10GY"  "2.5G"  "5G"   
##  [19] "7.5G"  "10G"   "2.5BG" "5BG"   "7.5BG" "10BG"  "2.5B"  "5B"    "7.5B" 
##  [28] "10B"   "2.5PB" "5PB"   "7.5PB" "10PB"  "2.5P"  "5P"    "7.5P"  "10P"  
##  [37] "2.5RP" "5RP"   "7.5RP" "10RP" 
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MunsellName }% use one of  RShowDoc("KEYWORDS")

