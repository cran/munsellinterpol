\name{MunsellToLuv}
\alias{MunsellToLuv}

\title{Convert a Munsell specification to CIE Luv coordinates}

\description{\code{MunsellToLuv} Converts a 
Munsell specification to CIE Luv coordinates, by interpolating over the 
extrapolated Munsell renotation data
}
\usage{MunsellToLuv( MunsellSpec,  white=c(95.047,100,108.883), adapt='Bradford', ... ) }

\arguments{
  \item{MunsellSpec}{ a numeric Nx3 matrix with HVC values in the rows, or a vector that can be converted to such a matrix, by row.\cr
\code{MunsellSpec} can also be a character vector with Munsell Notations,
which is converted to an Nx3 matrix using \code{\link{HVCfromMunsellName}()}}
  \item{white}{ XYZ for the destination white.  The default is Illuminant=D65, 2 observer }
  \item{adapt}{ method for chromatic adaptation, see \code{\link[spacesXYZ:adaptation]{CAT}()} for valid values.
  Also see \bold{Details}. }
  \item{...}{other parameters passed to \code{\link{MunsellToXYZ}()}  }
}

\value{
An Nx3 matrix with the Luv coordinates in each row.
The rownames are copied from input to output.
}

\details{
The conversion is done in 3 steps.
\itemize{
\item HVC  \rarrow  XYZ using \code{\link{MunsellToXYZ}()}
\item XYZ is adapted from Illuminant C to the given \code{white} using
the given chromatic adaptation method
\item XYZ  \rarrow  Luv using \code{\link[spacesXYZ:fromXYZ]{LuvfromXYZ}()} with the given \code{white}
}
}


\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{http://centore.isletech.net/~centore/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\author{Jose Gama and Glenn Davis}
\examples{
MunsellToLuv('7.6P 8.9/2.2')
}

\seealso{
\code{\link{LuvtoMunsell}()},
\code{\link{MunsellToXYZ}()},
\code{\link[spacesXYZ:adaptation]{CAT}()},
\code{\link[spacesXYZ:fromXYZ]{LuvfromXYZ}()}
}

\keyword{datasets}

