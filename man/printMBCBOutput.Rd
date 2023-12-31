\name{printMBCBOutput}
\alias{printMBCBOutput}
\title{MBCB - Model-Based Background Correction for Illumina Beadarray}
\description{
  This function is used to neatly output the values created by the other methods in the MBCB package.\cr
}
\usage{
  printMBCBOutput(sig, 
    average, 
    rmaBool, 
    npBool, 
    mleBool, 
    bayesBool, 
    gmleBool, 
    avgOutputFile, 
    detailOutputFile)
}
\arguments{
  \item{sig}{ A list corresponding to the backgrdoun-corrected signal values (as is generated by the mbcb.correct function). The list should have elements for all desired background correction methods named 'NP', 'RMA', 'Bayes', and/or 'MLE'. }
  \item{average}{ A list containing the average values (as is generated by the mbcb.correct function). The list should have elements for all desired background correction methods named 'NP', 'RMA', 'Bayes', and/or 'MLE'. }
  \item{npBool}{ A boolean value representing the desire to compute the non-parametric background correction values. }
  \item{rmaBool}{ A boolean value representing the desire to compute the RMA background correction values. }
  \item{mleBool}{ A boolean value representing the desire to compute the MLE background correction values. }
  \item{bayesBool}{ A boolean value representing the desire to compute the Bayes background correction values. }
  \item{gmleBool}{ A boolean value representing the desire to compute the GMLE background correction values. }
  \item{avgOutputFile}{ The \emph{base} file name to which suffixes and a file extension will be appended (i.e. 'C:/output'). These files will store the average values of each background correction method selected. }
  \item{detailOutputFile}{  The \emph{base} file name to which suffixes and a file extension will be appended (i.e. 'C:/output'). These files will store the background corrected intensities }
}
\value{
  This function prints corresponding CSV files based on which background-correction methods were selected.
}
\seealso{
  \code{\link{mbcb.main}}
}
\author{Yang Xie \email{ Yang.Xie@UTSouthwestern.edu}, Min Chen \email{ min.chen@phd.mccombs.utexas.edu}, Jeff Allen \email{ Jeffrey.Allen@UTSouthwestern.edu} }
\keyword{ models }
