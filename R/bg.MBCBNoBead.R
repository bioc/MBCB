#bg.MBCB.R
#  This file contains the actual implementation of the algorithms

# MBCB package
# Algorithm design by Yang Xie et. al.
# Package implementation and maintenence by Jeff Allen, (Bo Yao, the following maintainer)
# 2008 - University of Texas Southwestern Medical Center


mbcb.parseFile <- function(sigFile, conFile, isRawBead = FALSE){
  #if possible, update the progress bar. (only works if the user is using the GUI)
  TryNoError(setTkProgressBar(mbcbPB, title = "Reading...", label = "Reading Data.\nThis may take a while.", 0))
  
  #preliminary check
  if (is.null(conFile) && isRawBead){
    stop ("Error: In order to convert the raw, bead-level data into bead-type data you must have the signal file.")
    return()
  }
    
  #count the width
  fWidth <- length(scan(sigFile, skip=1,nlines=1, what="raw"))

  #create the list that will define the type of data I'm reading in.
  types <<- list("character")
  for (i in 2:fWidth) {
    types <<- c(types, "numeric")
  }
  
  if (is.null(conFile)){
    #no control file has been provided, don't import
    con <- NULL 
  }
  else{
      if (!isRawBead){

        TryNoError(setTkProgressBar(mbcbPB, title = "Reading...", label = "Reading Signal File.\nThis may take a while.", 0))
        sig <<- read.table(sigFile,sep="\t",header=TRUE, colClasses=types, row.names=1)    
        TryNoError(setTkProgressBar(mbcbPB, title = "Reading...", label = "Reading Control File.\nThis may take a while.", 0))
        con <<- read.table(conFile,sep="\t",header=TRUE, colClasses=types, row.names=1)  
        
        
      }else{
        warning ("WARNING: Bead-type input is still in development and will likely not function correctly.")
        TryNoError(setTkProgressBar(mbcbPB, title = "Reading...", label = "Reading Control File.\nThis may take a while.", 0))

        #try-catch here. If there's an error reading in via CSV, alert that they didn't click Bead-level data.
        con <<- read.csv(conFile, header=TRUE)                                                                           



        TryNoError(setTkProgressBar(mbcbPB, title = "Reading...", label = "Reading Signal File.\nThis may take a while.", 0))
        sig <- read.table(sigFile, header=TRUE, sep="\t")
                
                
        
        TryNoError(setTkProgressBar(mbcbPB, title = "Reading...", label = "Converting to bead-type data.", 0))
        
        #must convert to bead-type data
        sig.n <- sig[sig$Code%in%con[,2],]

   #     print (head(sig.n))

#        sigS1 <- read.table("C:/Documents and Settings/jeffreya/My Documents/School/UTS Research/Illumina Processing (08-07)/illumina/mcmc/1671586275_A_1.txt", header=T, sep="\t")
#        sigS2 <- read.table("C:/Documents and Settings/jeffreya/My Documents/School/UTS Research/Illumina Processing (08-07)/illumina/mcmc/1671586275_A_2.txt", header=T, sep="\t")

#        s1.n<- sigS1[sigS1$Code%in%con[,2],]
#        s2.n<- sigS2[sigS2$Code%in%con[,2],]
 

#        sigB <- rbind(s1.n, s2.n)
#        sig <- aggregate(sigB$Grn, list(id=sigB$Code), median)


        sig2 <<- aggregate(sig.n$Grn, list(id=sig.n$Code), median)

        sig <<- cbind(sig2[,2])
         head(sig)
        
        row.names(sig) <<- sig2[,1]

      }
 	}

 
  
  #multi-argument returns are deprecated, must return a named list.
  return (list(sig=sig, con=con))
}





mbcb.correct <- function(g, control, npBool=TRUE, rmaBool=FALSE, mleBool=FALSE, bayesBool=FALSE, gmleBool=FALSE, iter=500, burn=200, isRawBead=FALSE){

  #by default, there is no error  
	error <- FALSE
  #simple argument checking
	if (rmaBool != 1 && rmaBool != 0)
		error <- TRUE
	if (npBool != 1 && npBool != 0)
		error <- TRUE
	if (mleBool != 1 && mleBool != 0)
		error <- TRUE
	if (bayesBool != 1 && bayesBool != 0)
		error <- TRUE
	if (gmleBool != 1 && gmleBool != 0) 
        error <- TRUE
	if (error == TRUE){
		stop ("Invalid parameters passed to the background function. All boolean values must be 1 or 0.")
		return()
	}

  trialNames <<- colnames(g)
    
	#the number of arrays in the data file
	n.array <<- (dim(g)[2])
                                
	truesignal <<- NULL
	truenoise <<- NULL
	obsbead <<- NULL
	obsnc <<- NULL #global because the optimization takes a function as an argument and can't pass more arguments

	nbead <<- dim(g)[1]  ## sample size of beadarray
	SignalEst.rma <- SignalEst.np <- SignalEst.mle <- SignalEst.bayes <- SignalEst.gam <- matrix(0,dim(g)[1],n.array)

	#build the average table. Should be g[2]-1 wide (the number of trials in the dataset) and 3 tall
  avgVarsNP  <- avgVarsRMA <- avgVarsMLE <- avgVarsBayes <- avgVarsGMLE <- data.frame(matrix(0,3,dim(g)[2]-1))
  
	#may be possible to optimize the memory by putting signalEst into this table as well... Shouldn't really be too memory-intensive either way...
	rownames(avgVarsNP) <- c("muhat", "sigsqhat", "lenghat")
	rownames(avgVarsMLE) <- c("muhat", "sigsqhat", "lenghat")
	rownames(avgVarsRMA) <- c("muhat", "sigsqhat", "lenghat")
	rownames(avgVarsBayes) <- c("muhat", "sigsqhat", "lenghat")
  rownames(avgVarsGMLE) <- c("muhat", "sigsqhat", "lenghat")
  

  if (is.null(control)){ 
    #strictly do RMA background correction, nothing else.
    rmaBool <- TRUE 
    if (npBool || mleBool || mleBool){
      warning ("You cannot use MLE, NP, or Bayesian normalization without providing a signal file.")
      warning ("Using RMA Normalization")
      npBool <- mleBool <- bayesBool <- FALSE
    }
  
    for (i in 1:n.array) {
  	  TryNoError(setTkProgressBar(mbcbPB, title = "Processing...", label = paste(round((i-1)*90/n.array), "% Complete", sep=""), (i-1)/n.array))
  	  obsbead <<- g[,i]
   
		  ##2. RMA
		  rma <- bg.rma(pm=obsbead)
		  muhat.rma <- rma$para[2]
		  sigsqhat.rma <- rma$para[3]^2
		  lenghat.rma <- rma$para[1]
		  SignalEst.rma[,i] <- truesig(muhat.rma,sigsqhat.rma,lenghat.rma)
	    avgVarsRMA[,i] <- c(muhat.rma, sigsqhat.rma, lenghat.rma)
 	    
	  }
  	outRMA <- outNP <- outMLE <- outBayes <- data.frame()
  	
		outRMA <- (data.frame(g[,0], SignalEst.rma))
  
   # print(head(outNP))
    return (list(RMA = outRMA, NP = outNP, MLE = outMLE, Bayes = outBayes,
        AvgNP = NULL, AvgRMA = avgVarsRMA, AvgMLE = NULL, AvgBayes = NULL))
	  
  } 
  
  nnc <<- dim(control)[1]   ## sample size of negative control
  

  
	for (i in 1:n.array) {

    TryNoError(setTkProgressBar(mbcbPB, title = "Processing...", label = paste(round((i-1)*90/n.array), "% Complete", sep=""), (i-1)/n.array))
	  rkbead <- rank(g[,i],ties.method ="first")

    obsbead <<- g[,i]
    sortedObsbead <<- sort(g[,i])

    obsnc <<- control[,i]
    



 	  ##1. non-parametric
	  muhat.np <- mean(obsnc);
	  sigsqhat.np <- var(obsnc);
	  lenghat.np <- mean(obsbead)-muhat.np
	  avgVarsNP[,i] <- c(muhat.np, sigsqhat.np, lenghat.np)
	  SignalEst.np[,i] <- truesig(muhat.np,sigsqhat.np,lenghat.np)

	  if (rmaBool == 1){

		  ##2. RMA
		  rma <- bg.rma(pm=obsbead)
		  muhat.rma <- rma$para[2]
		  sigsqhat.rma <- rma$para[3]^2
		  lenghat.rma <- rma$para[1]
		  SignalEst.rma[,i] <- truesig(muhat.rma,sigsqhat.rma,lenghat.rma)
	    avgVarsRMA[,i] <- c(muhat.rma, sigsqhat.rma, lenghat.rma)
 	    
	  }


	  if (mleBool == 1){
		  sd.alpha <- sqrt((lenghat.np^2+sigsqhat.np)/nbead+sigsqhat.np/nnc)
		  mle.1stp <- optimize(liksig, lower=lenghat.np-2*sd.alpha, upper=lenghat.np+2*sd.alpha, maximum=TRUE,
				       mu=muhat.np, sig=sqrt(sigsqhat.np))
		  lenghat.1stp <- mle.1stp$maximum

	  

		  ##5. MLE
		  mle <- optim(c(lenghat.1stp,muhat.np,sqrt(sigsqhat.np)), liksig.all)
		  ## in case it does not converge
		  if (mle$convergence>0) {warning("MLE convergence error: ","i=",i,"\n");}
		  muhat.mle <- mle$par[2]
		  sigsqhat.mle <- mle$par[3]^2
		  lenghat.mle <- mle$par[1]
		  SignalEst.mle[,i] <- truesig(muhat.mle, sigsqhat.mle, lenghat.mle)
 	    avgVarsMLE[,i] <- c(muhat.mle, sigsqhat.mle, lenghat.mle) 
 	    
 	    
	}

	 if (bayesBool == 1){
		  ##6. bayes
		  bayes <- bg.mcmc(iter=iter, burn=burn)
		  muhat.bayes <- bayes[2]
		  sigsqhat.bayes <- bayes[3]^2
		  lenghat.bayes <- bayes[1]
		  SignalEst.bayes[,i] <- truesig(muhat.bayes, sigsqhat.bayes, lenghat.bayes)	 
	  	avgVarsBayes[,i] <- c(muhat.bayes, sigsqhat.bayes, lenghat.bayes)
	 }
	 if (gmleBool == 1) {
            if(mleBool != 1) {
              sd.alpha <- sqrt((lenghat.np^2 + sigsqhat.np)/nbead + 
                               sigsqhat.np/nnc)
              mle.1stp <- optimize(liksig, lower = lenghat.np - 
                                   2 * sd.alpha, upper = lenghat.np + 2 * sd.alpha, 
                                   maximum = TRUE, mu = muhat.np, sig = sqrt(sigsqhat.np))
              lenghat.1stp <- mle.1stp$maximum
              mle <- optim(c(lenghat.1stp, muhat.np, sqrt(sigsqhat.np)), 
                           liksig.all)
              muhat.mle <- mle$par[2]
              sigsqhat.mle <- mle$par[3]^2
              lenghat.mle <- mle$par[1]
            }
            initval <- c(lenghat.mle, muhat.mle^2/sigsqhat.mle, sigsqhat.mle/muhat.mle)
            ##initval <- c(lenghat.np, muhat.np^2/sigsqhat.np, sigsqhat.np/muhat.np) 
            mle <- optim(initval, liksigall.gamnoise, method = c("L-BFGS-B"), lower=0.01)

            gamLeng.mle <- mle$par[1]
            gamAlf.mle  <- mle$par[2]
            gamBet.mle  <- mle$par[3]
            SignalEst.gam[,i] <- estsig.gamnoise(gamAlf.mle,gamBet.mle,gamLeng.mle, sortedObsbead)[rkbead]
            avgVarsGMLE[, i] <- c(gamAlf.mle, gamBet.mle, gamLeng.mle)
        }
        
        message(paste(i, "/", n.array, " complete."), sep="");

  }
   
	outRMA <- outNP <- outMLE <- outGMLE<- outBayes <- data.frame()
	
  if (rmaBool ==TRUE)
		outRMA <- (data.frame(g[,0], SignalEst.rma))
	if (npBool == TRUE)
		outNP <- (data.frame(g[,0], SignalEst.np))	
	if (mleBool ==TRUE)
		outMLE <- (data.frame(g[,0], SignalEst.mle))		
	if (bayesBool == TRUE)
		outBayes <- (data.frame(g[,0], SignalEst.bayes))	
	if (gmleBool == TRUE) 
    outGMLE <- (data.frame(g[, 0], SignalEst.gam))
  
    return (list(RMA = outRMA, NP = outNP, MLE = outMLE, Bayes = outBayes, GMLE= outGMLE,
    AvgNP = avgVarsNP, AvgRMA = avgVarsRMA, AvgMLE = avgVarsMLE, AvgBayes = avgVarsBayes, AvgGMLE = avgVarsGMLE))
  
}



mbcb.main <- function(signal, control, npBool=TRUE, rmaBool=FALSE, mleBool=FALSE, bayesBool=FALSE, gmleBool=FALSE, paramEstFile="param-est", bgCorrectedFile="bgCorrected", iter=500, burn=200, normMethod="none", isRawBead=FALSE){
  
  if (is.null(control)){ 
    #strictly do RMA background correction, nothing else.
    rmaBool <- TRUE

    if (npBool || mleBool || mleBool){
      warning ("You cannot use MLE, NP, or Bayesian normalization without providing a signal file.")
      warning ("Using RMA Normalization")
      npBool <- mleBool <- bayesBool <- FALSE
    }
  }

  dat <- mbcb.correct(signal, control, npBool=npBool, rmaBool=rmaBool, mleBool=mleBool, bayesBool=bayesBool, gmleBool = gmleBool, iter=iter, burn=burn)


  
  if(normMethod=="none" || normMethod=="n" || normMethod==FALSE){
    #no normalization
    if (npBool){
        #take the log2 transformation
        dat$NP<- log2(dat$NP)
      }

      if (rmaBool){
        #take the log2 transformation
        dat$RMA<- log2(dat$RMA)
      }
  
      if (mleBool){
        #take the log2 transformation
        dat$MLE<- log2(dat$MLE)
      }

      if (bayesBool){
        #take the log2 transformation
        dat$Bayes<- log2(dat$Bayes)
      }
      if (gmleBool) {
        dat$GMLE <- log2(dat$GMLE)
      }

  }
  else if(normMethod=="quant" || normMethod=="q"){              ### Quantile-Quantile Normalization ###
    
    if (npBool){
        #take the log2 transformation
        lognew<- log2(dat$NP)
        
        #preserve row and column names. New version of normalize.quantiles loses them.
        rn <- row.names(lognew);
        cn <- colnames(lognew);
        q <- normalize.quantiles(as.matrix(lognew))
        colnames(q) <- cn;
        row.names(q) <- rn;
        dat$NP <- q;
        
      }

      if (rmaBool){
        #take the log2 transformation
        lognew<- log2(dat$RMA)
        
        #preserve row and column names. New version of normalize.quantiles loses them.
        rn <- row.names(lognew);
        cn <- colnames(lognew);
        q <- normalize.quantiles(as.matrix(lognew))
        colnames(q) <- cn;
        row.names(q) <- rn;
        dat$RMA <- q;
      }
  
      if (mleBool){
        #take the log2 transformation
        lognew<- log2(dat$MLE)
        
        
        #preserve row and column names. New version of normalize.quantiles loses them.
        rn <- row.names(lognew);
        cn <- colnames(lognew);
        q <- normalize.quantiles(as.matrix(lognew))
        colnames(q) <- cn;
        row.names(q) <- rn;
        dat$MLE <- q;
      }

      if (bayesBool){
        #take the log2 transformation
        lognew<- log2(dat$Bayes)
        
        
        #preserve row and column names. New version of normalize.quantiles loses them.
        rn <- row.names(lognew);
        cn <- colnames(lognew);
        q <- normalize.quantiles(as.matrix(lognew))
        colnames(q) <- cn;
        row.names(q) <- rn;
        dat$Bayes <- q;
      }
      if (gmleBool) {
        lognew <- log2(dat$GMLE)
        
        
        #preserve row and column names. New version of normalize.quantiles loses them.
        rn <- row.names(lognew);
        cn <- colnames(lognew);
        q <- normalize.quantiles(as.matrix(lognew))
        colnames(q) <- cn;
        row.names(q) <- rn;
        dat$GMLE <- q;
      }
         
    
      }
  else if(normMethod=="median" || normMethod=="m"){              ### Global/Median Normalization ###
    if (npBool){

      ### First take log2 transformation
      logsig<- log2(dat$NP)
      
  
      ## Get median of each array     
      med.logsig<- apply(logsig,2,median,na.rm=TRUE)
      G<- length(logsig[,1])
      
      ## Get mean of median of array    
      diff.d<- med.logsig-mean(med.logsig)
      
      ## Normalize to make the median of each array the same 
      normSig<- logsig-matrix(rep(diff.d,G),ncol=n.array,byrow=TRUE)
      dat$NP <- normSig
    }
    
    if (rmaBool){     
      #log2 transformation
      logsig<- log2(dat$RMA)
      
  
      ## Get median of each array     
      med.logsig<- apply(logsig,2,median,na.rm=TRUE)
      G<- length(logsig[,1])
      
      ## Get mean of median of array    
      diff.d<- med.logsig-mean(med.logsig)
      
      ## Normalize to make the median of each array the same 
      normSig<- logsig-matrix(rep(diff.d,G),ncol=n.array,byrow=TRUE)  ####################################################WIDTH???
      dat$RMA <- normSig
    }

    if (mleBool){
      #log2 transformation
      logsig<- log2(dat$MLE)
       
      ## Get median of each array     
      med.logsig<- apply(logsig,2,median,na.rm=TRUE)
      G<- length(logsig[,1])
      
      ## Get mean of median of array    
      diff.d<- med.logsig-mean(med.logsig)
      
      ## Normalize to make the median of each array the same 
      normSig<- logsig-matrix(rep(diff.d,G),ncol=n.array,byrow=TRUE)
      dat$MLE <- normSig
    }
    if (bayesBool){
      #log2 transformation
      logsig<- log2(dat$Bayes)
      
  
      ## Get median of each array     
      med.logsig<- apply(logsig,2,median,na.rm=TRUE)
      G<- length(logsig[,1])
      
      ## Get mean of median of array    
      diff.d<- med.logsig-mean(med.logsig)
      
      ## Normalize to make the median of each array the same 
      normSig<- logsig-matrix(rep(diff.d,G),ncol=n.array,byrow=TRUE)
      dat$Bayes <- normSig
    }
    if (gmleBool) {
            logsig <- log2(dat$GMLE)
            med.logsig <- apply(logsig, 2, median, na.rm = TRUE)
            G <- length(logsig[, 1])
            diff.d <- med.logsig - mean(med.logsig)
            normSig <- logsig - matrix(rep(diff.d, G), ncol = n.array, 
                byrow = TRUE)
            dat$GMLE <- normSig
    }
  }
  else{
    stop("Invalid normalization; must be \"none\", \"quant\", or \"median\".")
  }   

  avg <- list (NP = dat$AvgNP, RMA = dat$AvgRMA, MLE = dat$AvgMLE, Bayes=dat$AvgBayes, GMLE = dat$AvgGMLE)       # MUST UPDATE AVG after norm
  bgCorrSig <<- list (NP = dat$NP, RMA = dat$RMA, MLE = dat$MLE, Bayes = dat$Bayes,GMLE=dat$GMLE)
  printMBCBOutput (bgCorrSig, avg, rmaBool, npBool, mleBool, bayesBool, gmleBool, paramEstFile, bgCorrectedFile)
}



TryNoError <- function(expr)
{
    if (data.class(result<-try(expr,TRUE))=="try-error"){}
    else{
        return (result)
    }
}



printMBCBOutput <- function(sig, average, rmaBool, npBool, mleBool, bayesBool, gmleBool, avgOutputFile, detailOutputFile){
	#sig and average should be lists containing the data.frames with the information and should be named RMA, MLE, NP, and Bayes  so that they can be indexed as average$RMA, etc.
	#see the help file for more info
	
	dete <<- detailOutputFile;
  	
  #Try to update the progressBar (if using the GUI)
	TryNoError(setTkProgressBar(mbcbPB, title = "Writing...", label = "Preparing results...", .90))
  
	error <- FALSE
	
	#if the boolean value = T, then the value is present
	if (rmaBool != 1 && rmaBool != 0)
		error <- TRUE
	if (npBool != 1 && npBool != 0)
		error <- TRUE
	if (mleBool != 1 && mleBool != 0)
		error <- TRUE
	if (bayesBool != 1 && bayesBool != 0)
		error <- TRUE
	if (gmleBool != 1 && gmleBool != 0) 
    error <- TRUE
        
	if (error == TRUE){
		stop ("Invalid parameters passed to the print function. All boolean values must be TRUE or FALSE.")
		return()
	}
	
	
	

	trialCount <- length(trialNames)


	# a buffer to which we'll append the column names
	OutColnamesBuffNP <- list()
	OutColnamesBuffRMA <- list()
	OutColnamesBuffMLE <- list()
	OutColnamesBuffBayes <- list()
  OutColnamesBuffGMLE <- list()

	#count up the number of columns in the average file based on the presence of the boolean values
	if (rmaBool + npBool + mleBool + bayesBool + gmleBool <= 0){
		stop ("No methods selected to output!")
		return()
	}
	thisTrialName <- ""

	#start at 1, because that is the first column's title
	for (trial in 1:trialCount){
		thisTrialName <- trialNames[trial]

		if (rmaBool == TRUE){
			OutColnamesBuffRMA <- c(OutColnamesBuffRMA,paste(thisTrialName,"rma", sep="-"))
		}
		if (npBool == TRUE){
			OutColnamesBuffNP <- c(OutColnamesBuffNP,paste(thisTrialName,"np", sep="-"))
		}
		if (mleBool == TRUE){
			OutColnamesBuffMLE <- c(OutColnamesBuffMLE,paste(thisTrialName,"mle", sep="-"))
		}
		if (bayesBool == TRUE){
			OutColnamesBuffBayes <- c(OutColnamesBuffBayes,paste(thisTrialName,"bayes", sep="-"))

		}
		if (gmleBool == TRUE) {
         OutColnamesBuffGMLE <- c(OutColnamesBuffGMLE, paste(thisTrialName, 
         "gmle", sep = "-"))
    }
	}

	if (rmaBool == 1){
		TryNoError(setTkProgressBar(mbcbPB, title = "Writing...", label = "Writing RMA Results...", .92))
    sig$RMA <- format(sig$RMA, digits=3, nsmall=2)
    colnames(sig$RMA) <- OutColnamesBuffRMA
		write.csv(sig$RMA, paste(detailOutputFile, "-RMA.csv", sep=""), row.names=TRUE)
	}
	if (npBool == 1){
		TryNoError(setTkProgressBar(mbcbPB, title = "Writing...", label = "Writing NP Results...", .94))
		sig$NP <- format(sig$NP, digits=3, nsmall=2)
		colnames(sig$NP) <- OutColnamesBuffNP
		write.csv(sig$NP, paste(detailOutputFile, "-NP.csv",sep=""), row.names=TRUE)	
	}
	if (mleBool == 1){
		TryNoError(setTkProgressBar(mbcbPB, title = "Writing...", label = "Writing MLE Results...", .96))
		sig$MLE <- format(sig$MLE, digits=3, nsmall=2)
		colnames(sig$MLE) <- OutColnamesBuffMLE
		write.csv(sig$MLE,  paste(detailOutputFile, "-MLE.csv", sep=""), row.names=TRUE)
	}
	if (bayesBool == 1){
		TryNoError(setTkProgressBar(mbcbPB, title = "Writing...", label = "Writing Bayes Results...", .98))
		sig$Bayes <- format(sig$Bayes, digits=3, nsmall=2)
		colnames(sig$Bayes) <- OutColnamesBuffBayes
		write.csv(sig$Bayes, paste(detailOutputFile, "-Bayes.csv", sep=""), row.names=TRUE)
	}
	if (gmleBool == 1) {
        sig$GMLE <- format(sig$GMLE, digits = 3, nsmall = 2)
#        print(dim(sig$GMLE))
        colnames(sig$GMLE) <- OutColnamesBuffGMLE
        write.csv(sig$GMLE, paste(detailOutputFile, "-GMLE.csv", 
        sep = ""), row.names = TRUE)
  }

	TryNoError(setTkProgressBar(mbcbPB, title = "Writing...", label = "Writing Averages...", .99))
#  print ("check1")
	if (rmaBool == 1){
		colnames(average$RMA) <- OutColnamesBuffRMA
		write.csv(average$RMA, file=paste(avgOutputFile, "-RMA.csv", sep=""))
	}
	if (npBool == 1){		
		colnames(average$NP) <- OutColnamesBuffNP
		write.csv(average$NP, file=paste(avgOutputFile, "-NP.csv", sep=""))
	}
	if (mleBool == 1){
		colnames(average$MLE) <- OutColnamesBuffMLE
		write.csv(average$MLE, file=paste(avgOutputFile, "-MLE.csv", sep=""))
	}
	if (bayesBool == 1){
		colnames(average$Bayes) <- OutColnamesBuffBayes
		write.csv(average$Bayes, file=paste(avgOutputFile, "-Bayes.csv", sep=""))
	}
	if (gmleBool == 1) {
#      print(dim(average$GMLE))
        colnames(average$GMLE) <- OutColnamesBuffGMLE
        write.csv(average$GMLE, file = paste(avgOutputFile, "-GMLE.csv", 
            sep = ""))
   }

#     print ("check2")
	b <- TryNoError(close(mbcbPB))
}

bg.rma<- function (pm, n.pts = 2^14){
  max.density <- function(x, n.pts) {
    aux <- density(x, kernel = "epanechnikov", n = n.pts,
                   na.rm = TRUE)
    aux$x[order(-aux$y)[1]]
  }
  pmbg <- max.density(pm, n.pts)
  bg.data <- pm[pm < pmbg]
  pmbg <- max.density(bg.data, n.pts)
  bg.data <- pm[pm < pmbg]
  bg.data <- bg.data - pmbg
  bgsd <- sqrt(sum(bg.data^2)/(length(bg.data) - 1)) * sqrt(2)
  sig.data <- pm[pm > pmbg]
  sig.data <- sig.data - pmbg
  expmean <- max.density(sig.data, n.pts)
  alpha <- 1/expmean
  mubg <- pmbg

  a.rma<- pm-mubg-bgsd^2*alpha
  b.rma<- bgsd
  ### note: 3.992521e-22+1-1=0, 3.992521e-22+(1-1)=3.992521e-22
  Ex.rma<- a.rma+b.rma*(dnorm(a.rma/b.rma)-dnorm((pm-a.rma)/b.rma))/(pnorm(a.rma/b.rma)+(pnorm((pm-a.rma)/b.rma)-1))

  ### Note: in Y. Xie's notation, alpha is rate; but Sherry's is length.
  para<- cbind(alpha=1/alpha,mubg,bgsd)

  return(list(ex.rma=Ex.rma, para=para))
}

bg.mcmc<- function(iter=500, burn=200){
	nbead <- length(obsbead);
	nnc <- length(obsnc);
	
  if (iter <= burn){
		#out of order
		stop ("Iteration Count must be greater than the burn count.")
		return()
	}
	alpha <- rep(0,iter)
	mu <- rep(0,iter)
	sigma <- rep(0,iter)

	B <- 3
	G <- length(nbead)
	N <- length(nnc)

	accept <- matrix(0,B,iter)
	acc.rate <- matrix(0,B,iter)
	sample.sd <- matrix(0,B,iter)
	 
	sample.sd[,1] <- c(0.1,0.1,0.1)
	accept[,1] <- rep(0,B)

	acc.r <- rep(0,B)

	x.sample <- matrix(0,G,iter)

	mu[1] <- mean(obsnc)
	sigma[1] <- sd(obsnc)
	alpha[1] <- mean(obsbead) -mu[1]

	for(t in 2:iter)
	{
		
		alpha.new <- exp(rnorm(1,log(alpha[t-1]),sample.sd[1,t-1]))
		acc.u <- runif(1,0,1)
		acc.r <- exp( -liksig.all(c(alpha.new,mu[t-1],sigma[t-1])) + liksig.all(c(alpha[t-1],mu[t-1],sigma[t-1])) )
		accept[1,t] <- acc.u < acc.r
		acc.rate[1,t] <- ifelse(t > 20, sum(accept[1,(t-20):(t-1)])/20,0.5)
		sample.sd[1,t] <- sample.sd[1,t-1] * exp(0.2*(acc.rate[1,t]-0.2))
		alpha[t] <- ifelse(accept[1,t] == 1,alpha.new, alpha[t-1])

		mu.new <- rnorm(1,mu[t-1],sample.sd[2,t-1])
		acc.u <- runif(1,0,1)
		acc.r <- exp( -liksig.all(c(alpha[t],mu.new,sigma[t-1])) + liksig.all(c(alpha[t],mu[t-1],sigma[t-1])) )
		accept[2,t] <- (acc.u < acc.r)
		acc.rate[2,t] <- ifelse(t > 20, sum(accept[2,(t-20):(t-1)])/20,0.5)
		sample.sd[2,t] <- sample.sd[2,t-1] * exp(0.2*(acc.rate[2,t]-0.2))
		mu[t] <- ifelse(accept[2,t] == 1,mu.new, mu[t-1])

		sigma.new <- exp(rnorm(1,log(sigma[t-1]),sample.sd[3,t-1]))
		acc.u <- runif(1,0,1)
		acc.r <- exp( -liksig.all(c(alpha[t],mu[t],sigma.new)) + liksig.all(c(alpha[t],mu[t],sigma[t-1])) )
		accept[3,t] <- (acc.u < acc.r)
		acc.rate[3,t] <- ifelse(t > 20, sum(accept[3,(t-20):(t-1)])/20,0.5)
		sample.sd[3,t] <- sample.sd[3,t-1] * exp(0.2*(acc.rate[3,t]-0.2))
		sigma[t] <- ifelse(accept[3,t] == 1,sigma.new, sigma[t-1])

	}

	e.alpha<- mean(alpha[burn:iter])
	e.mu<- mean(mu[burn:iter])
	e.sigma<- mean(sigma[burn:iter])
	c(e.alpha,e.mu,e.sigma)
	return(c(e.alpha,e.mu,e.sigma))

}

####################################################
## log likelihood of bead data
liksig <- function(alf,mu,sig) { ## maximization
  length(obsbead)*(-log(alf)+sig*sig/2/alf/alf)+sum(-(obsbead-mu)/alf+
      pnorm((obsbead-mu-sig*sig/alf)/sig,log.p=TRUE))
}

## log likelihood of bead data AND negative control
## minimization b/c it returns (- loglik)
liksig.all <- function(para) { 

  alf <- abs(para[1])
  mu <- abs(para[2])
  sig <- abs(para[3])
  loglik <- NA

  if (sig > 0 && alf > 0 && mu > 0)
  {
    alfsq <- alf*alf
    sigsq <- sig*sig
 
    loglik <-  -(length(obsbead)*(-log(alf)+sigsq/2/alfsq)+
      sum(-(obsbead-mu)/alf+ pnorm((obsbead-mu-sigsq/alf)/sig,log.p=TRUE)) +
      (-length(obsnc)*log(sig))+sum(-(obsnc-mu)^2)/2/sigsq
      )
  }

  loglik
}

## estimate the true signal
esttruesig <- function(mu.est, sigsq.est, alf.est, sigest=NULL){
  if(is.null(sigest)) {
    a <- obsbead - mu.est - sigsq.est/alf.est
    b <- sqrt(sigsq.est)
    a.div.b <- a/b
    ##s.a.b <- mu.est/b+b/alf.est
    sigest <- a+b*(dnorm(a.div.b)/pnorm(a.div.b));
  }
  list(mse=mean((sigest-truesignal)^2), rsq=cor(sigest, truesignal)^2);
}
#========================== simulation ==========================

truesig <- function(mu.est, sigsq.est, alf.est){
    a <- obsbead - mu.est - sigsq.est/alf.est
    b <- sqrt(sigsq.est)
    a.div.b <- a/b
    ##s.a.b <- mu.est/b+b/alf.est
    a+b*(dnorm(a.div.b)/pnorm(a.div.b));
}

liksig.gamnoise <- function(theta,gamma.a,gamma.b) { ## maximization
  intfun <- function(y) {
    exp(  (gamma.a-1)*log(y) - y/newbet - lognmlconst  )
    ##exp( y/theta) * dgamma(y, gamma.a, scale=gamma.b)
  }

  newbet <- 1/(1/gamma.b - 1/theta)
  lognmlconst <- lgamma(gamma.a) + gamma.a*log(gamma.b)
  
  if (gamma.b<theta) {
    nbead*(-log(theta)+gamma.a*(log(theta)-log(theta-gamma.b))) - sum(sortedObsbead)/theta +
      sum(pgamma(sortedObsbead, gamma.a, scale=newbet, log.p=TRUE))
  }
  else {
    val <- integrate(intfun,0,sortedObsbead[1],stop.on.error=FALSE)$value    
    sumlog <- log(val)
    for (j in 2:nbead) {
       val <- val + integrate(intfun,sortedObsbead[j-1],sortedObsbead[j],stop.on.error=FALSE)$value
       sumlog <- sumlog+log(val)
    }
    nbead*(-log(theta) )  - sum(sortedObsbead)/theta + sumlog
  }
}


liksigall.gamnoise <- function(para) { ## minimization
  if (min(para)<0) {return(Inf)}
  else {
    theta <- para[1]
    gamma.a <- para[2]
    gamma.b <- para[3]
    -(liksig.gamnoise(theta,gamma.a,gamma.b)+ sum(dgamma(obsnc, gamma.a, scale=gamma.b, log=TRUE)) )
  }
}

estsig.gamnoise <- function(gamma.alf.est, gamma.bet.est, theta.est, obead){
  intfun <- function(y) {
    y*(gamma.alf.est-1) * exp(-y/newbet)  
  }
  intfun2 <- function(y) {
    y*(gamma.alf.est) * exp(-y/newbet)  
  }
  
  newbet <- theta.est*gamma.bet.est/(theta.est - gamma.bet.est)

  if (gamma.bet.est<theta.est) {
    sigest<- obead- newbet * gamma.alf.est *
      pgamma(obead, gamma.alf.est+1, scale=newbet)/pgamma(obead, gamma.alf.est, scale=newbet)
  }
  else {
    sigest<- obead - sapply(obead, function(intup) integrate(intfun2,0,intup)$value)/
      sapply(obead, function(intup) integrate(intfun,0,intup)$value)
  }

  return(sigest) 
}
