\name{MunsellToLab}
\alias{MunsellToLab}
\title{Convert a Munsell specification to CIE Lab coordinates}
\description{\code{MunsellToLab} Converts a 
Munsell specification to CIE Lab coordinates, by interpolating over the 
extrapolated Munsell renotation data
}

\usage{MunsellToLab( MunsellSpec, white='D65', adapt='Bradford', ... ) }
\arguments{
  \item{MunsellSpec}{ a numeric Nx3 matrix with HVC values in the rows, or a vector that can be converted to such a matrix, by row.\cr
\code{MunsellSpec} can also be a character vector with Munsell Notations,
which is converted to an Nx3 matrix HVC using \code{\link{HVCfromMunsellName}()} }

  \item{white}{XYZ for the source white - a numeric 3-vector with scaling irrelevant.
    \code{white} can also be a numeric 2-vector; which is interpreted as xy chromaticity and converted to XYZ with Y=100.
    \code{white} can also be the name of any standard illuminant recognized by
    \code{\link[spacesXYZ:standardXYZ]{spacesXYZ::standardXYZ}()} which is scaled so that Y=100.
    The default \code{'D65'} is converted to \code{c(95.047,100,108.883)}. }
  
  \item{adapt}{ method for chromatic adaptation, see \code{\link[spacesXYZ:adaptation]{spacesXYZ::CAT}()} for valid values.
  Also see \bold{Details}. }
  
  \item{...}{other parameters passed to \code{\link{MunsellToXYZ}()} }
}

\value{
An Nx3 matrix with the Lab coordinates in each row.
The rownames of Lab are copied from the input HVC matrix,
unless the rownames are \code{NULL} and
then the output rownames are the Munsell notations for HVC.
}

\details{
The conversion is done in these steps:
\itemize{
\item HVC \rarrow XYZ using \code{\link{MunsellToXYZ}()}
\item XYZ is adapted from Illuminant C to the given \code{white} using \code{\link[spacesXYZ:adapt]{spacesXYZ::adaptXYZ}()}
and the given chromatic adaptation method
\item XYZ \rarrow Lab using \code{\link[spacesXYZ:fromXYZ]{spacesXYZ::LabfromXYZ}()} with the given \code{white}
}
}

\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{https://www.munsellcolourscienceforpainters.com/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\author{Jose Gama and Glenn Davis}

\seealso{
\code{\link{LabToMunsell}()},
\code{\link{MunsellToXYZ}()},
\code{\link[spacesXYZ:fromXYZ]{spacesXYZ::LabfromXYZ}()},
\code{\link[spacesXYZ:standardXYZ]{spacesXYZ::standardXYZ}()},
\code{\link[spacesXYZ:adapt]{spacesXYZ::adaptXYZ}()},
\code{\link[spacesXYZ:adaptation]{spacesXYZ::CAT}()}
}


\examples{
MunsellToLab( c('7.6P 8.9/2.2', 'N 5/' ) )
##                      L       a         b
##  7.6P 8.9/2.2 89.19097 5.09879 -3.250468
##  N 5/         51.00375 0.00000  0.000000
}



