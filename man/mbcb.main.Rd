\name{mbcb.main}
\alias{mbcb.main}
\title{MBCB - Model-Based Background Correction for Illumina Beadarray}
\description{
  This is the main function which incorporates all the others. This should be the most straightforward and autonomous function in the MBCB package.\cr
}
\usage{
mbcb.main (signal, 
  control, 
  npBool=TRUE, 
  rmaBool=FALSE, 
  mleBool=FALSE, 
  bayesBool=FALSE, 
  gmleBool=FALSE, 
  paramEstFile="param-est", 
  bgCorrectedFile="bgCorrected", 
  iter=500, 
  burn=200, 
  normMethod="none", 
  isRawBead=FALSE)
}
\arguments{
  \item{signal}{ The data representing the signal file. }
  \item{control}{ The data representing the control file. }
  \item{npBool}{ A boolean value representing the desire to compute the non-parametric background correction values. }
  \item{rmaBool}{ A boolean value representing the desire to compute the RMA background correction values. }
  \item{mleBool}{ A boolean value representing the desire to compute the MLE background correction values. }
  \item{bayesBool}{ A boolean value representing the desire to compute the Bayes background correction values. }
  \item{gmleBool}{ A boolean value representing the desire to compute the GMLE background correction values. }
  \item{paramEstFile}{The \emph{base} file name to which suffixes and a file extension will be appended (i.e. 'C:/output'). These files will store the parameter estimates of each background correction method selected. }
  \item{bgCorrectedFile}{The \emph{base} file name to which suffixes and a file extension will be appended (i.e. 'C:/output'). These files will store the background corrected intensities }
  \item{iter}{ The iteration count; only used in Bayesian correction. }
  \item{burn}{ The number of iterations which will be burned; only used in Bayesian correction. }
  \item{normMethod}{ The normalization method to be used. By default, none will be applied. The choices are:\cr
    \var{none} - no normalization will be applied.\cr
    \var{quant} - Quantile-Quantile normalization will be applied (requires the affy and affyio packages be present).\cr
    \var{median} - Median or Global normalization will be applied.}
  \item{isRawBead}{A boolean value representing whether the input files are bead-level or bead-type. If the input is bead-level, set this value to True so that the raw bead-level values can be summarized to bead-type data.}
}
\value{
  The function will compute and output log2-tranformed values for the desired background correction methods. Rather than returning this data as R objects, the output is written to files based on the names given in \var{paramEstFile} and \var{bgCorrectedFile}. 
}
\note{
  You can use \code{\link{mbcb.parseFile}} to create the signal and control matrices from the given files.
}
\examples{
  data(MBCBExpressionData)
  #Use of global variables is obviously not ideal, but with R's pass-by-value 
  #   setup, we quickly run out of memory without using them on such large 
  #   arrays

  mbcb.main(expressionSignal, negativeControl);
}
\seealso{
  \code{\link{mbcb.correct}} \code{\link{mbcb.parseFile}}
}
\author{Yang Xie \email{ Yang.Xie@UTSouthwestern.edu}, Min Chen \email{ min.chen@phd.mccombs.utexas.edu}, Jeff Allen \email{ Jeffrey.Allen@UTSouthwestern.edu} }
\keyword{ models }
