\name{bg.rma}
\alias{bg.rma}
\title{MBCB - Robust Multi-Array Average Background Correction for Illumina Beadarray}
\description{
  This function provides the means of using only the Robust Multi-Array Average background correction method for the Illumina platform.\cr
}
\usage{
bg.rma(pm, n.pts = 2^14)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pm}{ The data to be background corrected. }
  \item{n.pts}{ Pertains to the specificity or accuracy of the method. 2^14 is the default. }
}
\value{
  A list with two values:
 \describe{
    \item{ex.rma:}{A vector containing the background-corrected values.}
    \item{para:}{The statistical summary of the computation.}
 }
}
\note{
  This method does not make use of negative control beads.
}
\seealso{
  \code{\link{mbcb.main}}
}
\examples{
  data(MBCBExpressionData)
  
  bg.rma(expressionSignal[,2]); 
}
\author{Yang Xie \email{ Yang.Xie@UTSouthwestern.edu}, Min Chen \email{ min.chen@phd.mccombs.utexas.edu}, Jeff Allen \email{ Jeffrey.Allen@UTSouthwestern.edu} }
\keyword{ models }
