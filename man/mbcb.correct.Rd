\name{mbcb.correct}
\alias{mbcb.correct}
\title{MBCB - Model-Based Background Correction for Illumina Beadarray}
\description{
  This function is used to background-correct the provided data using the selected correction methods. Normalization is \emph{not} applied.\cr
}
\usage{
mbcb.correct (g, 
  control, 
  npBool=TRUE, 
  rmaBool=FALSE, 
  mleBool=FALSE, 
  bayesBool=FALSE, 
  gmleBool=FALSE, 
  iter=500, 
  burn=200, 
  isRawBead=FALSE)
}
\arguments{
  \item{g}{ The data representing the signal file. }
  \item{control}{ The data representing the control file. }
  \item{npBool}{ A boolean value representing the desire to compute the non-parametric background correction values. Set to true if you wish to compute Non-Parametric background correction.}
  \item{rmaBool}{ A boolean value representing the desire to compute the RMA background correction values. Set to true if you wish to compute RMA background correction.}
  \item{mleBool}{ A boolean value representing the desire to compute the MLE background correction values. Set to true if you wish to compute MLE background correction.}
  \item{bayesBool}{ A boolean value representing the desire to compute the Bayes background correction values. Set to true if you wish to compute Bayes background correction.}
  \item{gmleBool}{ A boolean value representing the desire to compute the GMLE background correction values. Set to true if you wish to compute GMLE background correction.}
  \item{iter}{ The iteration count; only used in Bayesian correction. }
  \item{burn}{ The number of iterations which will be burned; only used in Bayesian correction. }
  \item{isRawBead}{A boolean value representing whether the input files are bead-level or bead-type. If the input is bead-level, set this value to True so that the raw bead-level values can be summarized to bead-type data.}
}
\value{
  This function returns a complex list which can be grouped into two categories:\cr
    Background-corrected Values:\cr
    \describe{
        \item{NP:}{The background corrected values of the Non-Parametric method (or an empty data.frame if this method was not used).}
        \item{RMA:}{The background corrected values of the RMA method (or an empty data.frame if this method was not used).}
        \item{MLE:}{The background corrected values of the MLE method (or an empty data.frame if this method was not used).}
        \item{Bayes:}{The background corrected values of the Bayesian method (or an empty data.frame if this method was not used).}
    }
    Average Values\cr
    \describe{
        \item{AvgNP:}{A data.frame of statistics pertaining to the average of the NP normalization method (or an empty data.frame if this method was not used).}
        \item{AvgRMA:}{A data.frame of statistics pertaining to the average of the RMA method (or an empty data.frame if this method was not used).}
        \item{AvgMLE:}{A data.frame of statistics pertaining to the average of the MLE method (or an empty data.frame if this method was not used).}   
        \item{AvgBayes:}{A data.frame of statistics pertaining to the average of the Bayesian method (or an empty data.frame if this method was not used).} 
    }
    These values have not been normalized or log2 transformed. See \code{\link{mbcb.main}} for such functionality.
}
\examples{
  data(MBCBExpressionData)
  mbcb.correct(expressionSignal, negativeControl);                            
}
\seealso{
  \code{\link{mbcb.main}}
}
\author{Yang Xie \email{ Yang.Xie@UTSouthwestern.edu}, Min Chen \email{ min.chen@phd.mccombs.utexas.edu}, Jeff Allen \email{ Jeffrey.Allen@UTSouthwestern.edu} }
\keyword{ models }
