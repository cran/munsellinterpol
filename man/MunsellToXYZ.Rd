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
The rownames of XYZ are copied from the input HVC matrix,
unless the rownames are \code{NULL} and
then the output rownames are the Munsell notations for HVC.
Note that these XYZ values are for viewing under Illuminant C, with Y=100.
There is no chromatic adaptation.
}

\details{
The conversion is done in these steps:
\itemize{
\item HVC  \rarrow  xyY using \code{\link{MunsellToxyY}()} and \code{...}. This xyY is for Illuminant C.
\item xyY  \rarrow  XYZ using \code{\link[spacesXYZ:toXYZ]{spacesXYZ::XYZfromxyY}()} 
}
}

\references{
Paul Centore 2014
The Munsell and Kubelka-Munk Toolbox
\url{https://www.munsellcolourscienceforpainters.com/MunsellAndKubelkaMunkToolbox/MunsellAndKubelkaMunkToolbox.html}
}

\author{Jose Gama and Glenn Davis}
\examples{
MunsellToXYZ('7.6P 8.9/2.2')
##                     X        Y        Z
##  7.6P 8.9/2.2 75.6055 74.61345 92.92308
}

\seealso{
\code{\link{MunsellToxyY}()},
\code{\link[spacesXYZ:toXYZ]{spacesXYZ::XYZfromxyY}()} 
}


\keyword{XYZ}

