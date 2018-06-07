\name{labtoMunsell}
\alias{labtoMunsell}
\title{Convert CIE Lab coordinates into a Munsell specification}
\description{\code{labtoMunsell} Converts CIE Lab coordinates to a Munsell 
specification, by interpolating over the extrapolated Munsell renotation data
}
\usage{labtoMunsell( Lab, white=c(95.047,100,108.883), adaption='bradford', ... ) }
\arguments{
  \item{Lab}{CIE Lab coordinates An Nx3 matrix, or a vector that can be converted to such a matrix. Each row of the matrix has Lab.}
  \item{white}{XYZ for the source white.  The default is Illuminant=D65, 2 observer}
  \item{adaption}{method for chromatic adaption, see \code{\link{adaption}} for valid values.
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
\item Lab  \rarrow  XYZ using \code{\link{lab2xyz}()} with the given \code{white}.
\item XYZ is then adapted from the given \code{white} to Illuminant C using
the given \code{\link{adaption}} method.
\item XYZ \rarrow  HVC using \code{\link{XYZtoMunsell}()}.
}
}


\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{http://www.99main.com/~centore/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\author{Jose Gama and Glenn Davis}
\examples{
labtoMunsell( c(74.613450, -20.4, 10.1) )
}

\keyword{datasets}

