\name{sRGBtoMunsell}
\alias{sRGBtoMunsell}
\title{Convert sRGB coordinates to a Munsell specification}

\description{\code{sRGBtoMunsell} Converts sRGB coordinates to a Munsell 
specification, by interpolating over the extrapolated Munsell renotation data
}

\usage{sRGBtoMunsell( sRGB, maxValue=255, adaption='bradford', ... ) }

\arguments{
  \item{sRGB}{ a numeric Nx3 matrix with sRGB coordinates in the rows, or a vector that can be converted to such a matrix, by row.
These are non-linear display values,
but they are not required to be integers.}

\item{maxValue}{maximum value of sRGB for display.
Other popular values are 1, 1023, and 65535.
Even when 1, they are still taken to be non-linear display values.}

\item{adaption}{method for chromatic adaption, see \code{\link{adaption}} for valid values.
Also see \bold{Details}. }
  
\item{...}{other parameters passed to \code{\link{XYZtoMunsell}()}    }
}

\value{
a numeric Nx3 matrix with HVC coordinates in the rows.
The rownames are copied from input to output.
}


\details{
The conversion is done in 3 steps.
\itemize{
\item sRGB  \rarrow  XYZ using \code{\link{srgb2xyz}()} with the given \code{maxValue}
\item XYZ is adapted from Illuminant D65 (from the sRGB standard) to Illuminant C
using the given \code{\link{adaption}} method
\item XYZ  \rarrow  HVC using \code{\link{XYZtoMunsell}()}
}
}


\references{
Wikipedia. \bold{sRGB}.
\url{https://en.wikipedia.org/wiki/SRGB}.

Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{http://www.99main.com/~centore/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{srgb2xyz}()},
\code{\link{adaption}},
\code{\link{XYZtoMunsell}()}
}

\examples{
sRGBtoMunsell( c(255,45,67) )
sRGBtoMunsell( c(1,0,1), maxValue=1 )   
}

\keyword{sRGB}

