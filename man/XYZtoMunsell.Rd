\name{XYZtoMunsell}
\alias{XYZtoMunsell}

\title{Convert XYZ coordinates to Munsell HVC}

\description{\code{XYZtoMunsell()} converts XYZ coordinates to Munsell HVC,
by interpolating over the  Munsell renotation data.

Note that \emph{by default} the input \code{XYZ} is relative to the white point of Illuminant C.
Since Illuminant C is obsolete, it is almost certain that an \code{XYZ} appearing
in modern calculations \emph{must} be chromatically adapted
to the white point of Illuminant C before interpolating.
This chromatic adaptation is done internally if the white point of XYZ is supplied as argument \code{white}.
The adaptation is performed by \code{\link[spacesXYZ:adaptXYZ]{spacesXYZ::adaptXYZ()}}.

If the input XYZ is computed from a reflectance spectrum,
then ideally the white point should be set to the XYZ computed when that reflectance spectrum
has constant value 1 - the so-called \emph{Perfect Reflecting Diffuser}.
The white points in published standards, e.g. D65 white = (95.047, 100, 108.883),
can differ from spectrally-derived white points
because of different wavelength step sizes and ranges.
For an example, see the vignette \href{../doc/soil.html}{\bold{Soil Colors}}.
}

\usage{XYZtoMunsell( XYZ, white=NULL, adapt='Bradford', ... ) }
\arguments{
\item{XYZ}{a numeric Nx3 matrix with CIE XYZ coordinates in the rows,
or a vector that can be converted to such a matrix, by row.
Ideally, XYZ should be scaled so that Y=100 for the perfect reflecting diffuser.
}

\item{white}{XYZ for the source white point.  
  If \code{NULL} then the input XYZ are for viewing in an environment with Illuminant C, and the input \code{XYZ} is not chromatically adapted.
  If \code{white} is not \code{NULL}, then the input \code{XYZ} is chromatically adapted from
  the given \code{white} to the white point of Illuminant C before conversion to HVC.  
  \code{white} can be a numeric 3-vector; by convention \code{Y=100} for Munsell work.
  \code{white} can also be a numeric 2-vector; which is interpreted as xy chromaticity and converted to XYZ with Y=100.
  \code{white} can also be the name of any standard illuminant recognized by
  \code{\link[spacesXYZ:standardXYZ]{spacesXYZ::standardXYZ}()} which is then scaled so that Y=100.
  
  If the input XYZ is computed from a reflectance spectrum,
  then ideally \code{white} should be set to the XYZ computed
  when that reflectance spectrum has constant value 1 (see \bold{Description}).
  }

\item{adapt}{method for chromatic adaptation, see \code{\link[spacesXYZ:adaptation]{spacesXYZ::CAT}()} for valid values.
  Also see \bold{Details}.
  If \code{white=NULL} then \code{adapt} is ignored.
  }

\item{...}{other parameters passed to \code{\link{xyYtoMunsell}()} }
}

\value{
an Nx3 matrix with Munsell HVC in the rows.
The rownames are copied from input to output, unless the input rownames are \code{NULL}
when the ouput rownames are set to the Munsell notations for HVC.\cr
In case of error, it returns \code{NULL}.
}


\section{Condition Signaling}{
For an overview see:
\href{../doc/munsellinterpol-guide.html#appendix-e---logging-and-condition-signaling}{\bold{Appendix E - Logging and Condition Signaling}} in the User Guide.

In case of an \code{ERROR} or \code{FATAL} event, a condition is signaled with one of these additional
S3 classes: \code{"invalid_argument"} or \code{"package_unavailable"},
whichever is appropriate.
The invalid argument, or the name of the unavailable package, is added to the condition object.
}


\details{
The conversion is done in these steps:
\itemize{
\item If \code{white} is not \code{NULL}, XYZ is chromatically adapted from the given \code{white} to the white point of Illuminant C using \code{\link[spacesXYZ:adaptXYZ]{spacesXYZ::adaptXYZ()}} and the given \code{adapt} method. 
\item XYZ  \rarrow  xyY using \code{\link[spacesXYZ:xyYfromXYZ]{spacesXYZ::xyYfromXYZ}()} 
\item xyY  \rarrow  HVC using \code{\link{xyYtoMunsell}()} and \code{...}.  This function does the interpolation over the  Munsell renotation data.
}
}

\references{
Paul Centore 2014.
The Munsell and Kubelka-Munk Toolbox.
\MKMTB
}

\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{MunsellToXYZ}()},
\code{\link{xyYtoMunsell}()},
\code{\link[spacesXYZ:xyYfromXYZ]{spacesXYZ::xyYfromXYZ}()},
\code{\link[spacesXYZ:adaptXYZ]{spacesXYZ::adaptXYZ()}}
}

\examples{
## make XYZ for an 18% gray card, under Illuminant D65
white.D65 = 100 * as.numeric( spacesXYZ::standardXYZ( 'D65' ) )
gray18 = 0.18 * white.D65
XYZtoMunsell( gray18 )
##                       H        V         C
##  7.6GY 4.9/0.65 37.6456 4.851927 0.6466169

## This is not neutral gray, as expected.  The hue is green-yellow.
## Chroma = 0.6466 is too large; it should be 0, up to numerical truncation.
## The problem is that the white points are not equal.
## Redo and supply the true white point.

XYZtoMunsell( gray18, white=white.D65 )
##         H        V C
##  N 4.9/ 0 4.851927 0        #  much better; an exact neutral

##  Although the hue H=0 here, it is actually undefined because C=0
##  (as in polar coordinates, theta is undefined when r=0).


## When XYZ=(0,0,0) then Value should be exactly 0.
## And when XYZ is equal to the white point, then Value should be exactly 10.
XYZtoMunsell( c( 0,0,0, white.D65 ), white=white.D65 )
##        H  V C
##  N 0/  0  0 0
##  N 10/ 0 10 0


##  By convention, Y=100 for the reference white; but this is not strictly required.
gray18
##              X  Y        Z
##  [1,] 17.10846 18 19.59894

XYZtoMunsell( gray18, white=gray18 )
##        H  V C
##  N 10/ 0 10 0
##  The Value is still exactly 10.
}

\keyword{XYZ}

