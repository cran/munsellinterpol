\name{NickersonColorDifference}
\alias{NickersonColorDifference}

\title{Calculate the Nickerson Color Difference between two Colors}

\description{Calculate the Nickerson Color Difference between two colors, given in Munsell HVC;
see \cite{Nickerson}.
}

\usage{NickersonColorDifference( HVC0, HVC1, symmetric=TRUE, coeffs=c(0.4,6,3) ) }

\arguments{
\item{HVC0}{a numeric Nx3 matrix with HVC values in the rows, or a vector that can be converted to such a matrix, by row.
\code{HVC0} can also be a numeric 3-vector with a single HVC, and it is then replicated
to match the size of \code{HVC1}.
\cr
\code{HVC0} can also be a character N-vector with Munsell Notations,
which is converted to an Nx3 matrix using \code{\link{HVCfromMunsellName}()}. 
}

\item{HVC1}{a numeric Nx3 matrix with HVC values in the rows, or a vector that can be converted to such a matrix, by row.
\code{HVC1} can also be a numeric 3-vector with a single HVC, and it is then replicated
to match the size of \code{HVC0}.
\cr
\code{HVC1} can also be a character N-vector with Munsell Notations,
which is converted to an Nx3 matrix using \code{\link{HVCfromMunsellName}()}.
}

\item{symmetric}{if \code{FALSE} then use the original Nickerson difference formula,
and if \code{TRUE} then use a symmetrized version;
see \bold{Details}.
}

\item{coeffs}{the coefficients for Hue, Value, and Chroma in the difference equation,
see \bold{Details}.  The defaults are the ones used by Nickerson, but tweaking them is possible
in special circumstances.  The argument must have length 3 and all coefficents positive.
}
}

\value{
A numeric N-vector with the pairwise differences, 
i.e. between row i of \code{HVC0} and  row i of \code{HVC1}.
}

\details{
If \code{HVC0}=\eqn{H_0,V_0,C_0} and If \code{HVC1}=\eqn{H_1,V_1,C_1} then the original Nickerson formula is:
\deqn{NCD(HVC0,HVC1) = 0.4 C_0 \Delta H ~+~ 6 \Delta V  ~+~  3 \Delta C}
where \eqn{ \Delta H = |H_0 - H_1| } (on the circle),
\eqn{ \Delta V = |V_0 - V_1| } and
\eqn{ \Delta C = |C_0 - C_1| }.
Unfortunately, if \code{HVC0} and \code{HVC1} are swapped, the color difference is different.
The first color is considered to be the \emph{reference color} and the second one is the \emph{test color}.
The difference is not symmetric.
\cr
Another problem is that the difference is not continuous when the second color is a neutral gray,
for rectangular coordinates on a plane of constant V.

Both of these problems are fixed with a slightly modified formula:
\deqn{NCD(HVC0,HVC1) = 0.4 \min(C_0,C_1) \Delta H ~+~ 6 \Delta V  ~+~  3 \Delta C}

For the first formula set \code{symmetric=FALSE} and for the second formula set \code{symmetric=TRUE}.
}

\references{
Nickerson, Dorothy.
The Specification of Color Tolerances.
\bold{Textile Research}.
Vol 6. pp. 505-514.
1936.
}

\author{Jose Gama and Glenn Davis}



\examples{
NickersonColorDifference( '7.6P 8.9/2.2', '8P 8.2/3'  )
##  [1] 6.952
}
\keyword{datasets}

