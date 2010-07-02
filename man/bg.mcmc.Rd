\name{bg.mcmc}
\alias{bg.mcmc}                                                                                                                                                                                                 
\title{MBCB - Bayesian Background Correction for Illumina Beadarray}
\description{
  This function provides the means of using only the MCMC (Bayesian) background correction method for the Illumina platform.\cr
}
\usage{
bg.mcmc(iter=500, burn=200)
}
\arguments{
  \item{iter}{ The iteration count for the Baysian correction. }
  \item{burn}{ The number of iterations to burn for the Bayesian correction. }
}
\value{
  This function returns an array of alpha, mu, and sigma values representing the values computed during the mcmc trial.
}
\note{
  This function makes use of two global variables. It will expect \code{obsbead} and \code{obsnc} are both established prior to calling this function. Obviously, this is not ideal, but R's pass-by-value functionality hindered the ability to pass these matrices as parameters. Using global variables increases performance substantially.
}
\seealso{
  \code{\link{mbcb.main}}
}
\examples{
  data(MBCBExpressionData)
  # Use of global variables is obviously not ideal, but with R's pass-by-value 
  #   setup, we quickly run out of memory without using them on such large 
  #   arrays

  #all of the signals from sample #2
  obsbead <<- expressionSignal[,2]
  
  #the negative control values for this sample
  obsnc <<- negativeControl[,2]         
  
  #compute the alpha, mu, and sigma values
  bg.mcmc();                            
}
\author{Yang Xie \email{ Yang.Xie@UTSouthwestern.edu}, Min Chen \email{ min.chen@phd.mccombs.utexas.edu}, Jeff Allen \email{ Jeffrey.Allen@UTSouthwestern.edu} }
\keyword{ models }
