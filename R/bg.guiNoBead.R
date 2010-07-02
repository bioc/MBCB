require (tcltk)
require (tcltk2)
mbcb.gui <- function() {
negControlFile <- tclVar("")
sigFile <- tclVar("")
bLvl <- tclVar("type")

#importedData <- tclVar()
#tclObj(importedData) <- c("file1.txt", "file2.txt")
haveNeg <- FALSE
haveSig <- FALSE

import.sig <- function(){
	sigFile <<- as.character(tclvalue(tkgetOpenFile()))
	#don't import a blank if the user cancelled the box
	if (!(sigFile == "")){
		tkconfigure(importTab.f.signal, text=tclVar(sigFile))
		haveSig <<- TRUE
	}
}

import.negCon <- function(){
	negControlFile <<- tclvalue(tkgetOpenFile())
	#don't import a blank if the user cancelled the box
	if (!(negControlFile == "")){
		tkconfigure(importTab.neg.negTxt, text=tclVar(negControlFile))
		haveNeg <<- TRUE
	}	

}


actions.norm <- function(){

  #check for errors
	if (!haveSig){
		tkmessageBox(message="You have not imported a signal file.")
		return()
	}

################ BEGIN BG CORRECT TAB     ##################
            
  BGCTab <- tk2notetab(main.nb, "BG Correct")
  tk2notetab.select(main.nb, "BG Correct")

  iter <- tclVar("500")
	burn <- tclVar("200")
	nPoints <- tclVar("14")
	
	MLE_checked <- tclVar(FALSE)
	bayes_checked <- tclVar(FALSE)
	RMA_checked <- tclVar(FALSE)
	NP_checked <- tclVar(TRUE)
	GMLE_checked <- tclVar(FALSE)
 
  if (!haveNeg){
  		tkmessageBox(message="You have not imported a negative control file.\nYou will only be able to use RMA Background-Correction.")
  		RMA_checked <- tclVar(TRUE)
  		NP_checked <- tclVar(FALSE)
  		negControlFile <<- NULL
  }
  else{
    	RMA_checked <- tclVar(FALSE)
  		NP_checked <- tclVar(TRUE)

  }	

	#declare the objects

  MLE <- tkframe(BGCTab, relief="groove", borderwidth=3, padx=5, pady=5)
	MLE.chk <- tkcheckbutton(MLE, text="MLE (Maximum Likelihood Estimation)", variable=MLE_checked, command=function() updateChecks())
	
	bayes <- tkframe(BGCTab, relief="groove", borderwidth=3, padx=5, pady=5)
	bayes.chk <- tkcheckbutton(bayes, text="Bayes (MCMC)", variable=bayes_checked, command=function() updateChecks())
	bayes.iter <- tkentry(bayes, textvariable = iter, width=5, state="disabled") 
	bayes.burn <- tkentry(bayes, textvariable = burn, width=5, state="disabled")

	RMA <- tkframe(BGCTab, relief="groove", borderwidth=3, padx=5, pady=5)
	RMA.chk <- tkcheckbutton(RMA, text="RMA (Robust Multi-Array Average)", variable=RMA_checked, command=function() updateChecks())
#	RMA.points <- tkentry(RMA, textvariable = nPoints, width=5, state="disabled")
	
	NP <- tkframe(BGCTab, relief="groove", borderwidth=3, padx=5, pady=5)
	NP.chk <- tkcheckbutton(NP, text="NP (Non-Parametric)", variable=NP_checked, command=function() updateChecks())
	
  GMLE <- tkframe(BGCTab, relief="groove", borderwidth=3, padx=5, pady=5)
	GMLE.chk <- tkcheckbutton(GMLE, text="GMLE", variable=GMLE_checked, command=function() updateChecks())
	
	startBut <- tkbutton(BGCTab, text="Next >>", command=function() initiate(), pady=4, width=13)

  if (!haveNeg){
    tkconfigure(bayes.chk, state="disabled")
    tkconfigure(RMA.chk, state="disabled")
    tkconfigure(NP.chk, state="disabled")		
    tkconfigure(MLE.chk, state="disabled")
    tkconfigure(GMLE.chk, state="disabled")
        		
	}
	else{
	  tkconfigure(bayes.chk, state="normal")
    tkconfigure(RMA.chk, state="normal")
    tkconfigure(NP.chk, state="normal")		
    tkconfigure(MLE.chk, state="normal")
    tkconfigure(GMLE.chk, state="normal")    
  }



	tkgrid (tklabel(BGCTab, text="Select your desired background correction methods", padx=5, pady=3), columnspan=3)
	tkgrid (MLE, row=1, column=0, sticky="nsew", columnspan=3)

  tkgrid (MLE.chk)
	
	tkgrid (bayes, row=2, column=0, sticky="nsew", columnspan=3)
	tkgrid (bayes.chk)
	tkgrid (tklabel(bayes, text="Iteration Count:"), bayes.iter)	
	tkgrid (tklabel(bayes, text="Burn Value:"), bayes.burn)

	tkgrid (RMA, row=3, column=0, sticky="nsew", columnspan=3)
	tkgrid (RMA.chk, columnspan=4)
#	tkgrid (tklabel(RMA, text="Points:      2^"), RMA.points)
	
	tkgrid (NP, row=4, column=0, sticky="nsew", columnspan=3)
	tkgrid (NP.chk)
	
	tkgrid (GMLE, row=5, column=0, sticky="nsew", columnspan=3)
	tkgrid (GMLE.chk);

	tkgrid (startBut, row=8, column=2, pady=15)
	
	

	initiate <- function(){
		if(as.numeric(tclvalue(iter)) <= as.numeric(tclvalue(burn))){
			tkmessageBox(message="The iteration count must be bigger than the burn value.")
		}
		else{
			normalizeData(as.numeric(tclvalue(MLE_checked)), as.numeric(tclvalue(bayes_checked)), as.numeric(tclvalue(RMA_checked)), as.numeric(tclvalue(NP_checked)), as.numeric(tclvalue(GMLE_checked)), as.numeric(tclvalue(iter)), as.numeric(tclvalue(burn)))
		}
	}
	
	updateChecks <- function(){
		startEn <- FALSE
		if (1 == as.numeric(tclvalue(MLE_checked))){
			startEn <- TRUE
		}

		if (1 == as.numeric(tclvalue(bayes_checked))){
			tkconfigure (bayes.iter, state="normal")
			tkconfigure (bayes.burn, state="normal")
			
			startEn <- TRUE
		}
		else{
			tkconfigure (bayes.iter, state="disabled")
			tkconfigure (bayes.burn, state="disabled")
		}

		if (1 == as.numeric(tclvalue(RMA_checked))){
			#tkconfigure (RMA.points, state="normal")
			startEn <- TRUE
		}
		else{
			#tkconfigure (RMA.points, state="disabled")
		}

		if (1 == as.numeric(tclvalue(NP_checked))){
			startEn <- TRUE
		}
		
		if (1 == as.numeric(tclvalue(GMLE_checked))){
			startEn <- TRUE
		}
		
		if (startEn){
			#at least one is checked
			tkconfigure(startBut, state="normal")
		}
		else{
			tkconfigure(startBut, state="disabled")
		}

		
	}	

}

normalizeData <- function(MLEBOOL, BAYESBOOL, RMABOOL, NPBOOL, GMLEBOOL, ITER, BURN){#ADD GMLE
  TryNoError(NormTab <- tk2notetab(main.nb, "Normalize"))
  tk2notetab.select(main.nb, "Normalize")
  normVar <- tclVar ("none")
  f <- tkframe(NormTab, relief="groove", borderwidth=3, padx=5, pady=5)
	noNorm <- tkradiobutton(f, variable = normVar, val="none", text="No Normalization", justify="left", width=25)
	medNorm <- tkradiobutton(f, variable = normVar, val="median", text="Global (Median) Normalization", justify="left", width=25)
	quantNorm <- tkradiobutton(f, variable = normVar, val="quant", text="Quantile-Quantile Normalization", justify="left", width=25)
	
	
	nextBut <- tkbutton(f, text="Next >>", command=function() exportData(MLEBOOL, BAYESBOOL, RMABOOL, NPBOOL, GMLEBOOL, ITER, BURN, as.character(tclvalue(normVar))), pady=4, width=13)

	tkgrid (f, sticky="nsew", pady=5, padx=5, ipady=3, ipadx=3)
	tkgrid (tklabel(f, text="Select your preferred normalization method:", padx=5, pady=3), columnspan=3)
	tkgrid (noNorm, columnspan=2)
	tkgrid (medNorm, columnspan=2 )
	tkgrid (quantNorm, columnspan=2)
	
	tkgrid (nextBut, row=7, column=1, pady=45)  

 }






 exportData <- function(MLEBOOL, BAYESBOOL, RMABOOL, NPBOOL, GMLEBOOL, ITER, BURN, NORM){
   outputFile <- tclVar("")	
	 detailedOutputFile <- tclVar("")
   
   outputSelected <- FALSE
 	 detailedOutputSelected <- FALSE
 
 	 TryNoError(expTab <- tk2notetab(main.nb, "Export"))
	 tk2notetab.select(main.nb, "Export")
   
  outEntry <- tkentry (expTab, textvariable=outputFile, width=28, state="disabled")
	detailedOutEntry <- tkentry (expTab, textvariable=detailedOutputFile, width=28, state="disabled")
	outBrowse <- tkbutton (expTab, text="Browse", command=function() setOut())
	detailedOutBrowse <- tkbutton (expTab, text="Browse", command=function() setDetOut())
  submitBut <- tkbutton(expTab, text="Go", command=function() eval(MLEBOOL, BAYESBOOL, RMABOOL, NPBOOL, GMLEBOOL, ITER, BURN, NORM, outFile, detailedOutputFile), pady=4, width=13)
	setOut <- function(){
		outFile <<- as.character(tclvalue(tkgetSaveFile(initialfile="parameters" , defaultextension="")))
		if (!(outFile == "")){
			tkconfigure(outEntry, text=tclVar(outFile))
			outputSelected <<- TRUE
		}
	}
	setDetOut <- function(){
		detailedOutputFile <<- as.character(tclvalue(tkgetSaveFile(initialfile="Background_Corrected_Intensities", defaultextension="")))
		if (!(detailedOutputFile == "")){
			tkconfigure(detailedOutEntry, text=tclVar(detailedOutputFile))
			detailedOutputSelected <<- TRUE
		}
	}	

  tkgrid (tklabel(expTab, text="Parameters: "), pady=5)
  tkgrid (outEntry, outBrowse,pady=5)
	tkgrid (tklabel(expTab, text="Background-Corrected Values: "), pady=5)
  tkgrid (detailedOutEntry, detailedOutBrowse,pady=5)
	 	
	tkgrid (submitBut, row=7, column=1, pady=45, columnspan=2) 
  
  tkgrid (tklabel(expTab, text="Note: The progress bar may show \"0%\" for a few "     ), columnspan=2)
  tkgrid (tklabel(expTab, text="minutes depending the correction methods you selected."), columnspan=2)
	
	  eval <- function (MLEBOOL, BAYESBOOL, RMABOOL, NPBOOL, GMLEBOOL, ITER, BURN, norm, out, detOut){
  	  if (!outputSelected || !detailedOutputSelected){		
  			tkmessageBox(message="You must select where you would like to save the output for both files.")
  		}
  		else{
        mbcbPB <<- tkProgressBar("Reading...", "Reading Data.\nThis may take a few minutes.", 0, 1, 0)
  			tkconfigure (main, cursor="watch")                                                           
        #give the window a second to load
        Sys.sleep(1)			
  			#correct the data

       isBead <- FALSE
       if (tclvalue(bLvl)=="bead"){
          #the data is in the bead form
          isBead <- TRUE
        }
        else{
          #the data is in the bead-type form
          isBead <- FALSE
        }
   

        #read in the files to the dataRes variable.


        if (is.null(negControlFile)){    
          dataRes <- mbcb.parseFile(as.character(sigFile), conFile = NULL, isRawBead = isBead)
          mbcb.main(dataRes$sig, NULL, npBool =NPBOOL, rmaBool = RMABOOL, mleBool = MLEBOOL, bayesBool = BAYESBOOL, gmleBool=GMLEBOOL, paramEstFile=out, bgCorrectedFile= detOut, iter = ITER, burn = BURN, normMethod = norm, isRawBead=isBead)       
        }
        else{                                                                                 
           dataRes <- mbcb.parseFile(as.character(sigFile), as.character(negControlFile), isRawBead = isBead)
           mbcb.main(dataRes$sig, dataRes$con, npBool =NPBOOL, rmaBool = RMABOOL, mleBool = MLEBOOL, bayesBool = BAYESBOOL, gmleBool = GMLEBOOL, paramEstFile=out, bgCorrectedFile= detOut, iter = ITER, burn = BURN, normMethod = norm, isRawBead=isBead)
        }
        tkconfigure (main, cursor="arrow")
  			tkmessageBox(message="Background-correction complete!")
        tkdestroy(main)
  		}
 		}
}



help.about <- function(){
	tkmessageBox(message="The background-correction algorithm behind MCBC was designed by Yang Xie\nThe R Package and GUI was developed by Jeff Allen\n\n2008 - UT Southwestern")

}


removeItem <- function(){
	sel <-as.integer(tkcurselection(importTab.f.dataList));
	TryNoError(tkdelete(importTab.f.dataList, sel))
	
}

Try <- function(expr)
{
    if (data.class(result<-try(expr,TRUE))=="try-error")
        tkmessageBox(title="An error has occured!",message=as.character(result),icon="error",type="ok")
    else
        return (result)
}

TryNoError <- function(expr)
{
    if (data.class(result<-try(expr,TRUE))=="try-error"){}
#        tkmessageBox(title="An error has occured!",message=as.character(result),icon="error",type="ok")
    else{
        return (result)
    }
}



############### BUILD THE ROOT GUI #############3

minHeight = 225
minWidth = 200
#importedData <- array()
#fileList <- array()

main <- tktoplevel()                   #Top level tk window used to load files
tkwm.resizable(main, FALSE, FALSE)
# ====Build Window Components====
#build parent frame

tkfocus(main)

main.nb <<- tk2notebook(main, c("Import", "BG Correct", "Normalize", "Export"))
importTab <- tk2notetab(main.nb, "Import")


importTab.f <- tkframe (importTab,relief="groove", borderwidth=3, pady=5)
main.statusFrame <- tkframe (main)
importTab.neg <- tkframe(importTab, relief="groove", borderwidth=3)
importTab.stage <- tkframe(importTab, relief="raised", width=10, height=50)

statusTxt <- tclVar("")
main.statusFrame.stat <- tklabel (main.statusFrame, text="Welcome to MBCB.", anchor="w", width=28, textvariable=statusTxt)
importTab.neg.negLbl <- tklabel (importTab.neg ,text="Negative Control File:")
importTab.neg.negTxt <- tkentry (importTab.neg, textvariable=negControlFile, width=35, state="disabled")
importTab.neg.browseNeg <- tkbutton  (importTab.neg,text="Browse", width=6 , command=function() import.negCon())

importTab.f.dataLbl <- tklabel (importTab.f, text="Signal file:")
#importTab.f.bead <- tkradiobutton(importTab.f, variable=bLvl, value="bead", text="Bead-Level")
#importTab.f.beadType <- tkradiobutton(importTab.f, variable=bLvl, value="type", text="Bead-Type")
importTab.f.signal <- tkentry (importTab.f, width=35, state="disabled")
importTab.f.browseSig <- tkbutton (importTab.f, text="Browse", command=function() import.sig())
importTab.process <- tkbutton (importTab,text="Next >>" , command=function() actions.norm(), width="13", pady="3")



tkgrid (main.nb, sticky="nsew")


tkgrid (importTab.neg, column=0, row=0, columnspan =2, pady=10 , ipady=5, ipadx=5)
tkgrid (importTab.neg.negLbl , column=0 , row=0 , sticky="nsew")
tkgrid (importTab.neg.negTxt , column=0 , row=1 , sticky="nsew")
tkgrid (importTab.neg.browseNeg , column=1 , row=1 , sticky="nsew" , padx=5)



tkgrid (importTab.f , column=0 , columnspan=2, row=3 , sticky="nsew", ipady=5, ipadx=5)
tkgrid (importTab.f.dataLbl , columnspan=4 , sticky="nsew", pady=0)
#tkgrid (importTab.f.bead, ttkseparator(importTab.f, orient="vertical"), importTab.f.beadType,  row=1, ipady=7)
tkgrid (importTab.f.browseSig, column=4, row =3, sticky="nsew", padx=5)
tkgrid (importTab.f.signal, column=0, columnspan=4, row=3, sticky="nsew")



tkgrid (importTab.process , column=1, columnspan=1, row=5 , sticky="nsew", pady=15 )   

tkgrid (importTab.stage , column=2 , row=0 , sticky="nswe" , rowspan=4)

tkgrid (main.statusFrame , column=0 , row=4 , sticky="swe")
tkgrid (main.statusFrame.stat , column=0 , row=0 , sticky="w" , padx=3)


################    END IMPORT TAB #################


}
