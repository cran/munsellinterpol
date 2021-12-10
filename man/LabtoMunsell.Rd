\name{LabtoMunsell}
\alias{LabtoMunsell}
\title{Convert CIE Lab coordinates into a Munsell specification}
\description{\code{LabtoMunsell} Converts CIE Lab coordinates to a Munsell 
specification, by interpolating over the extrapolated Munsell renotation data
}
\usage{LabtoMunsell( Lab, white=c(95.047,100,108.883), adapt='Bradford', ... ) }
\arguments{
  \item{Lab}{CIE Lab coordinates An Nx3 matrix, or a vector that can be converted to such a matrix. Each row of the matrix has Lab.}
  \item{white}{XYZ for the source white.  The default is Illuminant=D65, 2 observer}
  \item{adapt}{method for chromatic adaptation, see \code{\link[spacesXYZ:adaptation]{CAT}()} for valid values.
  Also see \bold{Details}.}
  \item{...}{other parameters passed to \code{\link{XYZtoMunsell}()} }
}

\value{
An Nx3 matrix with the Munsell HVC coordinates in each row.
The rownames are set to those of Lab.
}

\details{
The conversion is done in 3 steps.
\itemize{
\item Lab  \rarrow  XYZ using \code{\link[spacesXYZ:toXYZ]{XYZfromLab}()} with the given \code{white}.
\item XYZ is then adapted from the given \code{white} to Illuminant C using the given \code{adapt} method.
\item XYZ \rarrow  HVC using \code{\link{XYZtoMunsell}()}.
}
}


\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{https://www.munsellcolourscienceforpainters.com/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\seealso{
\code{\link[spacesXYZ:adaptation]{CAT}()},
\code{\link[spacesXYZ:toXYZ]{XYZfromLab}()},
\code{\link{XYZtoMunsell}()}
}

\author{Jose Gama and Glenn Davis}
\examples{
LabtoMunsell( c(74.613450, -20.4, 10.1) )
}

\keyword{datasets}

