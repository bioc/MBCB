\name{mbcb.gui}
\alias{mbcb.gui}
\title{MBCB - Model-Based Background Correction for Illumina Beadarray}
\description{
  This file provides the Graphical-User-Interface for the
  MBCB package.
}
\usage{
mbcb.gui()
}
\value{
  This GUI will allow you to easily input data files and a negative-control file and will output a log2-transformed background-corrected file.
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
\seealso{
  \code{\link{mbcb.main}}
}
\examples{
  data(MBCBExpressionData)

  #create files from the data provided in this package in the current directory
  write.table(expressionSignal, 'signal.txt', sep="\t");
  write.table(negativeControl, 'negative.control.txt', sep="\t");

  #open the GUI; you can use the files just created as input for the signal and
  #   negative control files.
  mbcb.gui();
}
\author{Yang Xie \email{ Yang.Xie@UTSouthwestern.edu}, Min Chen \email{ min.chen@phd.mccombs.utexas.edu}, Jeff Allen \email{ Jeffrey.Allen@UTSouthwestern.edu} }
\keyword{ models }
