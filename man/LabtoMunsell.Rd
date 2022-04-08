\name{LabToMunsell}
\alias{LabToMunsell}
\alias{LabtoMunsell}

\title{Convert CIE Lab coordinates to Munsell HVC}

\description{\code{LabToMunsell} Converts CIE Lab coordinates to Munsell HVC,
by interpolating over the extrapolated Munsell renotation data
}

\usage{LabToMunsell( Lab, white='D65', adapt='Bradford', ... )}
\arguments{
  \item{Lab}{CIE Lab coordinates.  An Nx3 matrix, or a vector that can be converted to such a matrix. Each row of the matrix has Lab.}
  
  \item{white}{XYZ for the source white - a numeric 3-vector with scaling irrelevant.
    \code{white} can also be a numeric 2-vector; which is interpreted as xy chromaticity and converted to XYZ with Y=100.
    \code{white} can also be the name of any standard illuminant recognized by
    \code{\link[spacesXYZ:standardXYZ]{spacesXYZ::standardXYZ}()} which is scaled so that Y=100.
    The default \code{'D65'} is converted to \code{c(95.047,100,108.883)}. }
  
  \item{adapt}{method for chromatic adaptation, see \code{\link[spacesXYZ:adaptation]{spacesXYZ::CAT}()} for valid values.
    Also see \bold{Details}.}
  
  \item{...}{other parameters passed to \code{\link{XYZtoMunsell}()} }
}

\value{
An Nx3 matrix with the Munsell HVC coordinates in each row.
The rownames are set to those of Lab.
}

\details{
The conversion is done in these steps.
\itemize{
\item Lab  \rarrow  XYZ using \code{\link[spacesXYZ:toXYZ]{spacesXYZ::XYZfromLab}()} with the given \code{white}.
\item XYZ is then adapted from the given \code{white} to Illuminant C using the given \code{adapt} method.
\item XYZ \rarrow  HVC using \code{\link{XYZtoMunsell}()}.
}
}

\note{
The case of the letter 't' in the function name was recently changed from lower to upper.
The function \code{LabtoMunsell()} is equivalent but deprecated and provided for a limited time.
}

\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{https://www.munsellcolourscienceforpainters.com/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\seealso{
\code{\link{MunsellToLab}()},
\code{\link{XYZtoMunsell}()},
\code{\link[spacesXYZ:standardXYZ]{spacesXYZ::standardXYZ}()},
\code{\link[spacesXYZ:adaptation]{spacesXYZ::CAT}()},
\code{\link[spacesXYZ:toXYZ]{spacesXYZ::XYZfromLab}()}
}

\author{Jose Gama and Glenn Davis}
\examples{
LabToMunsell( c(74.613450, -20.4, 10.1,    80, 0, 0) )
##                      H        V       C
##  3.1G 7.4/3.6 43.13641 7.379685 3.62976
##  N 7.9/        0.00000 7.945314 0.00000
}

\keyword{datasets}

