\name{LuvToMunsell}
\alias{LuvToMunsell}
\alias{LuvtoMunsell}
\title{Convert CIE Luv coordinates to Munsell HVC}

\description{\code{LuvToMunsell} Converts CIE Luv coordinates to Munsell HVC,
by interpolating over the extrapolated Munsell renotation data
}

\usage{LuvToMunsell( Luv, white='D65', adapt='Bradford', ... ) }
\arguments{
  \item{Luv}{CIE Luv coordinates An Nx3 matrix, or a vector that can be converted to such a matrix. Each row of the matrix has Luv.}
  
  \item{white}{XYZ for the source white - a numeric 3-vector with scaling irrelevant.
    \code{white} can also be a numeric 2-vector; which is interpreted as xy chromaticity and converted to XYZ with Y=100.
    \code{white} can also be the name of any standard illuminant recognized by
    \code{\link[spacesXYZ:standardXYZ]{spacesXYZ::standardXYZ}()} which is scaled so that Y=100.
    The default \code{'D65'} is converted to \code{c(95.047,100,108.883)}. }
  
  \item{adapt}{method for chromatic adaptation, see \bold{Details}}
  
  \item{...}{other parameters passed to \code{\link{XYZtoMunsell}()} }  
}

\value{
An Nx3 matrix with the Munsell HVC coordinates in each row.
The rownames are set to those of Luv.
}

\details{
The conversion is done in these steps:
\itemize{
\item Luv \rarrow XYZ using \code{\link[spacesXYZ:toXYZ]{spacesXYZ::XYZfromLuv}()} with the given \code{white}.
\item XYZ is then adapted from the given \code{white} to Illuminant C using
the given chromatic adaptation method, see \code{\link[spacesXYZ:adaptation]{spacesXYZ::CAT}()}.
\item XYZ \rarrow HVC using \code{\link{XYZtoMunsell}()}
}
}

\note{
The case of the letter 't' in the function name was recently changed from lower to upper.
The function \code{LuvtoMunsell()} is equivalent but deprecated and provided for a limited time.
}


\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{https://www.munsellcolourscienceforpainters.com/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\seealso{
\code{\link{XYZtoMunsell}()},
\code{\link[spacesXYZ:toXYZ]{spacesXYZ::XYZfromLuv}()},
\code{\link[spacesXYZ:standardXYZ]{spacesXYZ::standardXYZ}()},
\code{\link[spacesXYZ:adaptation]{spacesXYZ::CAT}()}
}


\author{Jose Gama and Glenn Davis}
\examples{
LuvToMunsell( c( 74.613450, -5.3108, 10.6,  55, 0, 0 ) )
##                      H        V        C
##  8.7GY 7.4/1.2 38.6599 7.383948 1.229376
##  N 5.4/         0.0000 5.395003 0.000000
}

