\name{LuvtoMunsell}
\alias{LuvtoMunsell}
\title{Convert CIE Luv coordinates into a Munsell specification}
\description{\code{LuvtoMunsell} Converts CIE Luv coordinates into a Munsell 
specification, by interpolating over the extrapolated Munsell renotation data
}

\usage{LuvtoMunsell( Luv, white=c(95.047,100,108.883), adapt='Bradford', ... ) }
\arguments{
  \item{Luv}{CIE Luv coordinates An Nx3 matrix, or a vector that can be converted to such a matrix. Each row of the matrix has Luv.}
  \item{white}{XYZ for the reference white.  The default is Illuminant=D65, 2 observer}
  \item{adapt}{method for chromatic adaptation, see \bold{Details}}
  \item{...}{other parameters passed to \code{\link{XYZtoMunsell}()} }  
}

\value{
An Nx3 matrix with the Munsell HVC coordinates in each row.
The rownames are set to those of Luv.
}

\details{
The conversion is done in 3 steps.
\itemize{
\item Luv \rarrow XYZ using \code{\link[spacesXYZ:toXYZ]{XYZfromLuv}()} with the given \code{white}.
\item XYZ is then adapted from the given \code{white} to Illuminant C using
the given chromatic adaptation method, see \code{\link[spacesXYZ:adaptation]{CAT}()}.
\item XYZ \rarrow HVC using \code{\link{XYZtoMunsell}()}
}
}


\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{https://www.munsellcolourscienceforpainters.com/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\seealso{
\code{\link{XYZtoMunsell}()},
\code{\link[spacesXYZ:toXYZ]{XYZfromLuv}()},
\code{\link[spacesXYZ:adaptation]{CAT}()}
}


\author{Jose Gama and Glenn Davis}
\examples{
LuvtoMunsell( c( 74.613450, -5.3108, 10.6 ) )
}
\keyword{datasets}

