\name{Munsell2xy}
\alias{Munsell2xy}
\docType{data}
\title{The Munsell HVC to xy 3D Lookup Table}

\description{
This is the discrete data for the Munsell Renotation System,
which is often considered to be the most perceptually uniform
color atlas.
It was created by the NBS and OSA from "3,000,000 color judgments"
by 40 observers.
}

% \usage{<none>}

\format{
A data frame with 4995 observations of the following 6 variables.
\itemize{
\item{\code{H}} the Munsell Hue. Each \code{H} is a multiple of 2.5 and in the interval (0,100].
\item{\code{V}} the Munsell Value.  Each \code{V} is an integer from 1 to 10, or one of 0.2, 0.4, 0.6, 0.8
\item{\code{C}} the Munsell Chroma.  Each \code{C} is a positive even integer.
\item{\code{x}} the x chromaticity coordinate, for Illuminant C.
\item{\code{y}} the y chromaticity coordinate, for Illuminant C.
\item{\code{real}} a logical value. If \code{TRUE} then x,y were published, otherwise they have been extrapolated.
}
Note that the luminance factor Y is *not* here,
since Y is a simple function of V,
see \code{\link{YfromV}()}.
}

\details{
All the (x,y) data here comes from the file \code{all.dat} downloaded from
Rochester Institute of Technology, see \bold{Source}.
The file \code{real.dat} is a subset, and contains the (x,y) published
in \cite{Newhall, et. al. (1943)}.
These rows have \code{real=TRUE} and are only for Value \eqn{\ge} 1.
There are 2734 of these.
\cr
Similarly, for Value<1 (very dark colors), (x,y) data from
the paper \cite{Judd et. al. (1956)} also have \code{real=TRUE}.
There are 355 of these.
\cr
So \code{all.dat} has 4995 colors,
of which 2734+355=3089 are "real" colors,
and the remaining 1906 are extrapolated.
I am confident that the extrapolation was done by
\emph{Schleter et. al. (1958)} at the NBS,
and put online by the Rochester Institute of Technology.
For more details, and the abstract of the 1958 article,
see the \bold{munsellinterpol User Guide}.
}

\note{
For the purpose of this package, I have found that the extrapolated (x,y)
for V\eqn{\ge}1 work well.
But for V<1 they did not work so well, and I was able to get better results with
my own extrapolation.
Moreover, to get reliable results in this package for high Chroma,
it was necessary to extrapolate past the data in \code{all.dat}.
}

\author{
Glenn Davis
}

\source{
Rochester Institute of Technology.
Program of Color Science.
Munsell Renotation Data.
\url{http://www.rit.edu/science/pocs/renotation}
}

\references{
Newhall, Sidney M.,
Dorothy Nickerson,
Deane B. Judd.
\bold{Final Report of the O.S.A. Subcommitte on the Spacing of the Munsell Colors}.
Journal of the Optical Society of America.
Vol. 33.  No. 7.  pp. 385-418.
July 1943.

Judd, Deane B. and Gunter Wyszecki.
\bold{Extension of the Munsell Renotation System to Very Dark Colors}.
Journal of the Optical Society of America.
Vol. 46.  No. 4.  pp. 281-284.
April 1956.


Schleter, J. C, D. B. Judd, D. B., H. J. Keegan.
\bold{Extension of the Munsell Renotation System (Abstract)}.
J. Opt. Soc. Am.
Vol 48. Num. 11.
pp. 863-864.
presented at the Forty-Third Annual Meeting of the Optical Society of America.
Statler Hilton Hotel, Detroit, Michigan.  October 9, 10, and 11, 1958.
}

\seealso{
\code{\link{YfromV}()}
}

\examples{
str(Munsell2xy)

##  'data.frame':	4995 obs. of  6 variables:
##   $ H   : num  32.5 35 37.5 37.5 40 40 42.5 42.5 45 45 ...
##   $ V   : num  0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 ...
##   $ C   : num  2 2 2 4 2 4 2 4 2 4 ...
##   $ x   : num  0.713 0.449 0.262 -0.078 0.185 -0.257 0.144 -0.235 0.117 -0.209 ...
##   $ y   : num  1.414 1.145 0.837 2.16 0.676 ...
##   $ real: logi  FALSE FALSE FALSE FALSE TRUE FALSE ...
}
\keyword{datasets}
