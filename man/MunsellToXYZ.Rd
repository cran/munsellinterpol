\name{MunsellToXYZ}
\alias{MunsellToXYZ}
\title{Convert a Munsell specification to CIE XYZ coordinates}
\description{\code{MunsellToXYZ} Converts a 
Munsell specification to XYZ coordinates, by interpolating over the 
extrapolated Munsell renotation data
}

\usage{MunsellToXYZ( MunsellSpec, ... ) }
\arguments{
\item{MunsellSpec}{ a numeric Nx3 matrix with HVC values in the rows, or a vector that can be converted to such a matrix, by row.\cr
\code{MunsellSpec} can also be a character vector with Munsell Notations,
which is converted to an Nx3 matrix using \code{\link{HVCfromMunsellName}()}. }
\item{...}{other parameters passed to \code{\link{MunsellToxyY}()}        
}
}

\value{
an Nx3 matrix with XYZ values in the rows.
The rownames are copied from the input HVC matrix.
Exception: If the input matrix rownames are \code{NULL},
then the output rownames are the Munsell notation.
Note that these XYZ values are for viewing under Illuminant C, with Y=100.
}

\details{
This function calls  \code{\link{MunsellToxyY}()} 
and \code{\link{xyY2XYZ}()}.
}

\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{http://www.99main.com/~centore/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\author{Jose Gama and Glenn Davis}
\examples{
MunsellToXYZ('7.6P 8.9/2.2')
}

\seealso{
\code{\link{MunsellToxyY}()},
\code{\link{xyY2XYZ}()}
}


\keyword{XYZ}

