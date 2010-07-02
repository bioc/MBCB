\name{mbcb.parseFile}
\alias{mbcb.parseFile}
\title{MBCB - Model-Based Background Correction for Illumina Beadarray}
\description{
  This function is used to read the files provided into a format which will be usable by the MBCB package (data.frame).\cr
}
\usage{
mbcb.parseFile(sigFile, conFile, isRawBead = FALSE)
}
\arguments{
  \item{sigFile}{ The file-name (character string) representing the signal file. }
  \item{conFile}{ The file-name (character string) representing the control file. }
  \item{isRawBead}{A boolean value representing whether the input files are bead-level or bead-type. If the input is bead-level, set this value to True so that the raw bead-level values can be summarized to bead-type data.}
}
\value{
  This function will return a list containing two data frames: \var{sig} and \var{con}.
  \describe{
    \item{sig:}{Represents the signal file.}
    \item{con:}{Represents the negative control data.frame.}
  }
}
\note{
  The input files should be tab-delimited files in the following form:\cr\cr
	\var{signal.txt}\cr
\preformatted{				g.1	g.2	g.3	g.4
10181072_239_rc-S	160.3	776.4	135.8	407.7
10181072_290-S		138.1	219.8	122.1	142.4
...			...	...	...	...
}
	\var{neg-con.txt}\cr
\preformatted{id	g.1	g.2	g.3	g.4
50133	127	213.5	82	103
50315	232	295	143.5	156
...	...	...	...	...
}
Note that both of these examples come from a test involving four trials. Your data could use whatever number you'd like.
}
\examples{
  data(MBCBExpressionData)

  # Create files from the data provided in this package in the current 
  #     directory
  # Obviously, this is the opposite of what the function does, but we need to 
  #     write sample files for the sake of the demonstration
  # The signal and negative control files can be used by calling the above 
  #     command without the writing and reading of the data.
  write.table(expressionSignal, 'signal.txt', sep="\t");
  write.table(negativeControl, 'negative.control.txt', sep="\t");

  #read in those files just created.
  data <- mbcb.parseFile('signal.txt', 'negative.control.txt');
  signal <- data$sig;
  negCon <- data$con;       
}
\seealso{
  \code{\link{mbcb.correct}}
}
\author{Yang Xie \email{ Yang.Xie@UTSouthwestern.edu}, Min Chen \email{ min.chen@phd.mccombs.utexas.edu}, Jeff Allen \email{ Jeffrey.Allen@UTSouthwestern.edu} }
\keyword{ models }
