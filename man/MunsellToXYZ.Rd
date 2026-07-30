\name{MunsellToXYZ}
\alias{MunsellToXYZ}

\title{Convert a Munsell specification to CIE XYZ coordinates}

\description{\code{MunsellToXYZ} Converts a 
Munsell specification to XYZ coordinates, by interpolating over the 
extrapolated Munsell renotation data.

Note that \emph{by default} the returned \code{XYZ} is relative to the white point of Illuminant C.
Since Illuminant C is obsolete, it is almost certain that if this \code{XYZ} is used in subsequent
calculations, it must be chromatically adapted to a different white point.

But if the white point for subsequent calculations is supplied,
then the chromatic adaptation is done internally by
\code{\link[spacesXYZ:adaptXYZ]{spacesXYZ::adaptXYZ()}}.
}

\usage{MunsellToXYZ( MunsellSpec, white=NULL, adapt='Bradford', ... ) }
\arguments{
\item{MunsellSpec}{ a numeric Nx3 matrix with HVC values in the rows, or a vector that can be converted to such a matrix, by row.\cr
\code{MunsellSpec} can also be a character vector with Munsell Notations,
which is converted to an Nx3 matrix using \code{\link{HVCfromMunsellName}()}. }


\item{white}{XYZ for the output white point.  
  If \code{NULL} then the output XYZ are for viewing in an environment with Illuminant C,
  and the output \code{XYZ} is not chromatically adapted.
  If \code{white} is not \code{NULL}, then the interpolated \code{XYZ} is chromatically adapted from
  the white point of Illuminant C to the given white point.  
  \code{white} can be a numeric 3-vector; by convention \code{Y=100} for Munsell work.
  \code{white} can also be a numeric 2-vector; which is interpreted as xy chromaticity and converted to XYZ with Y=100.
  \code{white} can also be the name of any standard illuminant recognized by
  \code{\link[spacesXYZ:standardXYZ]{spacesXYZ::standardXYZ}()} which is then scaled so that Y=100.
  }

\item{adapt}{method for chromatic adaptation; for valid values see \code{\link[spacesXYZ:adaptation]{spacesXYZ::CAT}()}.
  Also see \bold{Details}.
  If \code{white=NULL} then \code{adapt} is ignored.
  }

\item{...}{other parameters passed to \code{\link{MunsellToxyY}()}        
}
}

\value{
an Nx3 matrix with XYZ values in the rows.
The rownames of XYZ are copied from the input HVC matrix,
unless the rownames are \code{NULL} and
then the output rownames are the Munsell notations for HVC.
}


\section{Condition Signaling}{
For an overview see:
\href{../doc/munsellinterpol-guide.html#appendix-e---logging-and-condition-signaling}{\bold{Appendix E - Logging and Condition Signaling}} in the User Guide.

In case of an \code{ERROR} event, a condition is signaled with one of these additional
S3 classes: \code{"invalid_argument"} or \code{"package_unavailable"},
whichever is appropriate.
The invalid argument, or the name of the unavailable package, is added to the condition object.
}




\details{
The conversion is done in these steps:
\itemize{
\item HVC  \rarrow  xyY using \code{\link{MunsellToxyY}()} and \code{...}. This xyY is for Illuminant C.
\item xyY  \rarrow  XYZ using \code{\link[spacesXYZ:toXYZ]{spacesXYZ::XYZfromxyY}()} 
\item If \code{white} is not \code{NULL}, XYZ is chromatically adapted from the white point of Illuminant C to the given \code{white} using the given \code{adapt} method.
}
}

\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\MKMTB
}

\author{Jose Gama and Glenn Davis}
\examples{
## make expected XYZ for a Perfect Reflecting Diffuser, under Illuminant D65
white.D65 = 100 * as.numeric( spacesXYZ::standardXYZ( 'D65' ) )
white.D65
##  [1]  95.047 100.000 108.883

## convert Munsell 'N 10/' - the Perfect Reflecting Difuser - to XYZ
MunsellToXYZ( "N 10/" )
##               X   Y        Z
##  N 10/ 98.03984 100 118.1157
##  this is the white point of Illuminant C, and not equal to white.D65

## Redo and supply the true white point.
MunsellToXYZ( "N 10/", white=white.D65 )
##             X   Y       Z
##  N 10/ 95.047 100 108.883
##
##  now we have a match


MunsellToXYZ('7.6P 8.9/2.2')
##                     X        Y        Z
##  7.6P 8.9/2.2 75.6055 74.61345 92.92308
}

\seealso{
\code{\link{MunsellToxyY}()},
\code{\link[spacesXYZ:toXYZ]{spacesXYZ::XYZfromxyY}()},
\code{\link[spacesXYZ:adaptXYZ]{spacesXYZ::adaptXYZ()}}
}


\keyword{XYZ}

