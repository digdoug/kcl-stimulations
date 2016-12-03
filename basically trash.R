#####basically trash


agg <- function(x,primerName,exptName,c){

   concs <- function(x,primerName,exptName,c){
  
   			x <- x[which(x$primer==primerName),]
            x <- x[x$Value != 999,]
            
           ##until the expt parse is run this will have to stay in this format

   	 	   exptName$exptname <- sapply(exptName$names,function(x){strsplit(x,split='_')[[1]][1]})
   	       exptName$concentrate <- sapply(exptName$names,function(x){strsplit(x,split='_')[[1]][2]})
  	       exptName$timepoints <- sapply(exptName$names,function(x){strsplit(x,split='_')[[1]][3]})
           exptName$day <- sapply(exptName$names,function(x){strsplit(x,split='_')[[1]][5]})
           

  	   	   x <- x[which(!is.na(match(x$timepoints,exptName$timepoints))),]
           x <- x[which(!is.na(match(x$exptname,exptName$exptname))),]
           x <- x[which(!is.na(match(x$concentrate,exptName$concentrate))),]

   		   use <- all(!is.na(x$timepoints == exptName[2]))
 		   concs <- tapply(x$realConc[use],factor(x$timepoints[use]),print,na.rm=T)

		   concs <- concs[order(as.numeric(names(concs)))] 
		   concs
		   }  
		   
	    test <- concs(x[[1]],primerName,exptName[[1]],c[[1]])
	    test2 <- concs(x[[2]],primerName,exptName[[2]],c[[2]])
	    test3 <- concs(x[[3]],primerName,exptName[[3]],c[[3]])
		
		listc1 <- normalize_conc(test,c[[1]])
		listc2 <- normalize_conc(test2,c[[2]])
		listc3 <- normalize_conc(test3,c[[3]])
		
        
		newlist <- mapply(c,listc1,listc2,listc3)
		means <- unlist(lapply(newlist,mean))
	    stderrs <- unlist(lapply(newlist,function(x){sd(x,na.rm=TRUE)/sqrt(length(x))}))
	    means <- means / means[1]
	    stderrs <- stderrs / means[1]
	    df <- data.frame(means=means,stderrs=stderrs,timepoints=c(0,20,40,60,120,240,360))
	    
	       ggplot(df,aes(x=timepoints,means)) + 
           geom_point() + 
           geom_line() +
           scale_x_continuous(name='Time (min)') +
           scale_y_continuous(name='Fold Induction') +
           geom_errorbar(aes(ymin=means-stderrs,ymax=means+stderrs),width=2) +
           opts(title = primerName)    
   }
  


   concs <- function(x,primerName,exptName,c){
  
   			x <- x[which(x$primer==primerName),]
            x <- x[x$Value != 999,]
            
           ##until the expt parse is run this will have to stay in this format

   	 	   exptName$exptname <- sapply(exptName$names,function(x){strsplit(x,split='_')[[1]][1]})
   	       exptName$concentrate <- sapply(exptName$names,function(x){strsplit(x,split='_')[[1]][2]})
  	       exptName$timepoints <- sapply(exptName$names,function(x){strsplit(x,split='_')[[1]][3]})
           exptName$day <- sapply(exptName$names,function(x){strsplit(x,split='_')[[1]][5]})
           

  	   	   x <- x[which(!is.na(match(x$timepoints,exptName$timepoints))),]
           x <- x[which(!is.na(match(x$exptname,exptName$exptname))),]
           x <- x[which(!is.na(match(x$concentrate,exptName$concentrate))),]

   		   use <- all(!is.na(x$timepoints == exptName[2]))
 		   concs <- tapply(x$realConc[use],factor(x$timepoints[use]),print,na.rm=T)

		   concs <- concs[order(as.numeric(names(concs)))] 
		   concs
		   }  
		   
	    test <- concs(c1,'Fosb',fiveminpulse_c1_b1,controls1)
	    test2 <- concs(c2,'Fosb',fiveminpulse_c2_b4,controls2)
	    test3 <- concs(c5,'Fosb', fiveminpulse_c5_b5_6,controls3)

normalize_conc <- function(x,c){
		mtx <- as.list(array(NA,dim=length(x)))
	for(i in 1:length(mtx)){ 	
		mtx[[i]] <- array(NA,dim=length(x[[i]]))
		}
      
        mtx[[1]] <- x[[1]] / c[1]        
		mtx[[2]] <- x[[2]] / c[2]	
		mtx[[3]] <- x[[3]] / c[3]	
		mtx[[4]] <- x[[4]] / c[4]	
		mtx[[5]] <- x[[5]] / c[5]	
		mtx[[6]] <- x[[6]] / c[6]
		mtx[[7]] <- x[[7]] / c[7]	
mtx
		}
			    	
		listc1 <- normalize_conc(test,controls1)
		listc2 <- normalize_conc(test2,controls2)
		listc3 <- normalize_conc(test3,contronls3)
		
        
		newlist <- mapply(c,listc1,listc2,listc3)
		means <- unlist(lapply(newlist,mean))
	    stderrs <- unlist(lapply(newlist,function(x){sd(x,na.rm=TRUE)/sqrt(length(x))}))
	    means <- means / means[1]
	    stderrs <- stderrs / means[1]
	    df <- data.frame(means=means,stderrs=stderrs,timepoints=c(0,20,40,60,120,240,360))
	    
	       ggplot(df,aes(x=timepoints,means)) + 
           geom_point() + 
           geom_line() +
           scale_x_continuous(name='Time (min)') +
           scale_y_continuous(name='Fold Induction') +
           geom_errorbar(aes(ymin=means-stderrs,ymax=means+stderrs),width=2) +
           opts(title = primerName)    
   }
  
  
   
        agg(chips,'Pcsk1',experiments,controls)   
        
        
  norms <-     function(x){	
	 arr <- array(NA,length(x))
	
	arr[1] <- x[[1]] / x[[1]]
	arr[2] <- x[[2]] / x[[1]]
	arr[3] <- x[[3]] / x[[1]]
	arr[4] <- x[[4]] / x[[1]]
	arr[5] <- x[[5]] / x[[1]]
	arr[6] <- x[[6]] / x[[1]]
	arr[7] <- x[[7]] / x[[1]]
	arr
}
  
     df <- data.frame(means=norms(means),stderrs=norms(stderrs))
  
     		
		
		
		lapply(means,function(x){df <- x[[1]] 
			x / df})
		
	normalize_conc	<- function(x,c){
		mtx <- as.list(array(NA,dim=length(x)))
	for(i in 1:length(mtx)){ 	
		mtx[[i]] <- array(NA,dim=length(x[[i]]))
		}
        mtx[[1]] <- x[[1]] / c[1]
		mtx[[2]] <- x[[2]] / c[2]	
		mtx[[3]] <- x[[3]] / c[3]	
		mtx[[4]] <- x[[4]] / c[4]	
		mtx[[5]] <- x[[5]] / c[5]	
		mtx[[6]] <- x[[6]] / c[6]	
		mtx[[7]] <- x[[7]] / c[7]	
mtx
		}
	
	
	}
mtx
}		
		
		

norm(test,controls1)
	  	
	norm	<- function(x,c){
			for(i in 1:length(x)){
				   # for(i in 1:length(c)){	
x[i] / c[1]			
			# }	
		   }		
	      }
	     
	     test[1][[1]] / controls1[1]
	       lapply()
     	     
 
	  	
	     norm(test,controls1)
	     		
		dim(test)	
			test
	test <- concs(c1_club,'Acan',fiveminpulse_c1,controls1)
	# concs_5min_c1_b2 <- concs(c1_club,'Acan',fiveminpulse_c1_b2,controls2_test)		

	head(c1_clu)
	
	min60_c1 <- c1_club[c1_club$sampleName == 'KClpulse-5min_55mM_60_b1_011311',]
	min60_c1 <- min60_c1[order(min60_c1$primer),]
######filter


	 filter_inductions <- function(x){
   x=x[,apply(x,2,function(x){any(is.finite(x))})]
     result <- matrix(data=NA,ncol=ncol(x),nrow=nrow(x))
        colnames(result) <- colnames(x)
   	       for(i in 1:ncol(x)){
   		      if (any(x[,i] >= 3)){
   	       result[,i] <- x[,i] 
   }
  } 
   	   result[,apply(result,2,function(x){any(!is.na(x))})]
 }
 \
 
 
      plotPulse('Egr1',c6,sixtyminpulse_55mM_b3,controls1)
     
     
      plot_gg <- function(primerName,x,exptName,c,replot=F,myColor="black",maxY=NA,lines=NULL){


   	   		# limit x to only the relevant rows; then set the x-values for plotting (designed for timecourses)
   			x <- x[which(x$primer==primerName),]

   	 	   exptName$exptname <- sapply(exptName$names,function(x){strsplit(x,split='_')[[1]][1]})
   	       exptName$concentrate <- sapply(exptName$names,function(x){strsplit(x,split='_')[[1]][2]})
  	       exptName$timepoints <- sapply(exptName$names,function(x){strsplit(x,split='_')[[1]][3]})
           exptName$day <- sapply(exptName$names,function(x){strsplit(x,split='_')[[1]][5]})
           

   	   	   x <- x[which(!is.na(match(x$timepoints,exptName$timepoints))),]
           x <- x[which(!is.na(match(x$exptname,exptName$exptname))),]
           x <- x[which(!is.na(match(x$concentrate,exptName$concentrate))),]
           x <- x[which(!duplicated(x$sampleName)),]

   			
   			if(any(duplicated(paste(x$exptname,x$primer,x$concentrate,x$timepoints)))){noBioRep <- F
   			} else {noBioRep <- TRUE}

			# across whatever replicates are available, calculate means and stderrs (normalized to the controls that were input)
   			stderr <- function(x){sd(x,na.rm=TRUE)/sqrt(length(x))}
   			use <- which(!is.na(x$timepoints))
   			means <- tapply(x$realConc[use],factor(as.numeric(x$timepoints[use])),mean,na.rm=T)
   			stderrs <- tapply(x$realConc[use],factor(as.numeric(x$timepoints[use])),stderr)
   			
   			means <- means/c
   			stderrs <- stderrs / c
   			stderrs <- stderrs / means[1]   	
   			means <- means/means[1]
           
			if(noBioRep == T){stderrs[1:length(stderrs)] <- 0} # there are no tech reps and therefore no std errs
 
            m = data.frame(means)
            colnames(m) = 'means'
             
            s = data.frame(stderrs)
            colnames(s) = 'stderrs'
              
            df <- cbind(m,s)
           
            if(all(!is.na(means))){

           ggplot(df,aes(x=as.numeric(rownames(df)),means),main=primerName) + 
           geom_point() + 
           geom_line() +
           scale_x_continuous(name='Time (min)') +
           scale_y_continuous(name='Fold Induction') +
           geom_errorbar(aes(ymin=means-stderrs,ymax=means+stderrs),width=2) +
           opts(title = primerName)     
           }
       }   
       
       
    plotAllPrimers_gg <- function(x,exptname,c,fN){

multiplot <- function(list,plotlist=NULL,cols) {
	require(grid)	
plots <- c(list,plotlist)
plotCols <- cols
plotRows <- 2
grid.newpage()
	pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x,y){
    	viewport(layout.pos.row = x,layout.pos.col=y)
    }
   for (i in 1:length(plots)){
      curRow = ceiling(i/plotCols)
  	  curCol = (i-1) %% plotCols + 1
  	  print(plots[[i]], vp = vplayout(curRow,curCol))  	
  }  
}   
  primes=unique(x$primer)
  plot_list <- lapply(primes,plotPulse,x,exptname,c)
  
   multiplot(plot_list[1:8],plotlist=NULL,cols=4)
   quartz.save(file=paste(fN,"_1-8.pdf",sep=""),type="pdf")
   multiplot(plot_list[9:16],plotlist=NULL,cols=4)
   quartz.save(file=paste(fN,"_9-16.pdf",sep=""),type="pdf")
   multiplot(plot_list[17:24],plotlist=NULL,cols=4)
   quartz.save(file=paste(fN,"_17-24.pdf",sep=""),type="pdf")
   multiplot(plot_list[25:32],plotlist=NULL,cols=4)
   quartz.save(file=paste(fN,"_25-32.pdf",sep=""),type="pdf")
   multiplot(plot_list[33:40],plotlist=NULL,cols=4)
   quartz.save(file=paste(fN,"_33-40.pdf",sep=""),type="pdf")
   multiplot(plot_list[41:48],plotlist=NULL,cols=4)
   quartz.save(file=paste(fN,"_41-48.pdf",sep=""),type="pdf")
   multiplot(plot_list[49:56],plotlist=NULL,cols=4)
   quartz.save(file=paste(fN,"_49-56.pdf",sep=""),type="pdf")
   multiplot(plot_list[57:64],plotlist=NULL,cols=4)
   quartz.save(file=paste(fN,"_57-64.pdf",sep=""),type="pdf")
   multiplot(plot_list[65:72],plotlist=NULL,cols=4)
   quartz.save(file=paste(fN,"_65-72.pdf",sep=""),type="pdf")
   multiplot(plot_list[73:80],plotlist=NULL,cols=4)
   quartz.save(file=paste(fN,"_73-80.pdf",sep=""),type="pdf")
   multiplot(plot_list[81:88],plotlist=NULL,cols=4)
   quartz.save(file=paste(fN,"_81-88.pdf",sep=""),type="pdf")
   multiplot(plot_list[89:length(plot_list)],plotlist=NULL,cols=4)
   quartz.save(file=paste(fN,"_89-end.pdf",sep=""),type="pdf")   
   }
       



##Experiment 1
##KCl Pulse length 60 min timepoints(min): 0,20,40,60,120,240,360
#######Just a note for this experiment and the rest, it doesn't matter that there is only one zero hour listed here, or one of the bio reps. The code will take in the other zero hr/biorep and incorporate that. This is essentially just giving it a backbone to build on.

   sixtyminpulse_55mM_b3 <- data.frame(names=c("0hr_55mM_0_t1_050312", "KClpulse-60min_55mM_20_b3_050312","KClpulse-60min_55mM_40_b3_050312","KClpulse-60min_55mM_60_b3_050312","KClpulse-60min_55mM_120_b3_050312","KClpulse-60min_55mM_240_b3_050312","KClpulse-60min_55mM_360_b3_050312"),timepoints=c(0,20,40,60,120,240,360)) 
   
   	controls1 <- plotCNTLs(c6,CNTLs_c5and6, sixtyminpulse_55mM_b3)
   	
   	
  fiveminpulse_c1 <- data.frame(names=c("0hr_55mM_0_b1_011311", "KClpulse-5min_55mM_20_b1_011311","KClpulse-5min_55mM_40_b1_011311","KClpulse-5min_55mM_60_b1_011311","KClpulse-5min_55mM_120_b1_011311","KClpulse-5min_55mM_240_b1_011311","KClpulse-5min_55mM_360_b1_011311"),timepoints=c(0,20,40,60,120,240,360))  
  
  	controls1c <- plotCNTLs(c1,CNTLs_c1,fiveminpulse_c1)
##############
#########
primes=unique(c5$sampleName)

  fiveminpulse_c5_10mM <- data.frame(names=c("0hr_55mM_0_t1_041812", "KClpulse-5min_10mM_20_b1_041812","KClpulse-5min_10mM_40_b1_041812","KClpulse-5min_10mM_60_b1_041812","KClpulse-5min_10mM_120_b1_041812","KClpulse-5min_10mM_240_b1_041812","KClpulse-5min_10mM_360_b1_041812"),timepoints=c(0,20,40,60,120,240,360))  

  	controls10mM_c5 <- plotCNTLs(c5,CNTLs_c5and6,fiveminpulse_c5_10mM)
  	
  	 mM_10 <-  sapply(primes,plotPulse,c5,fiveminpulse_c5_10mM,controls10mM_c5)
plotPulse

##############
#########
  fiveminpulse_c5_30mM <- data.frame(names=c("0hr_55mM_0_t1_041812", "KClpulse-5min_30mM_20_b1_041812","KClpulse-5min_30mM_40_b1_041812","KClpulse-5min_30mM_60_b1_041812","KClpulse-5min_30mM_120_b1_041812","KClpulse-5min_30mM_240_b1_041812","KClpulse-5min_30mM_360_b1_041812"),timepoints=c(0,20,40,60,120,240,360))  
  
    	controls_c5_30mM <- plotCNTLs(c5,CNTLs_c5and6,fiveminpulse_c5_30mM)
    	
    	mM_30 <-  sapply(primes,plotPulse,c5,fiveminpulse_c5_30mM,controls_c5_30mM)

##############
#########
    fiveminpulse_c5_55mM <- data.frame(names=c("0hr_55mM_0_t1_041812", "KClpulse-5min_55mM_20_b1_041812","KClpulse-5min_55mM_40_b1_041812","KClpulse-5min_55mM_60_b1_041812","KClpulse-5min_55mM_120_b1_041812","KClpulse-5min_55mM_240_b1_041812","KClpulse-5min_55mM_360_b1_041812"),timepoints=c(0,20,40,60,120,240,360))  

  	controls_c5_55mM <- plotCNTLs(c5,CNTLs_c5and6,fiveminpulse_c5_55mM)

    	mM_55 <-  sapply(primes,plotPulse,c5,fiveminpulse_c5_55mM,controls_c5_55mM)



   plotAllPrimers(x=c1,expt= fiveminpulse_c1,exptCntl= controlsz,primers=primerSet_c1,fN="maxt",newWindow=TRUE)     

##Experiment 1a
##KCl Pulse length 60 min timepoints(min): 0,20,40,60,120,240,360

   sixtyminpulse_55mM_b4 <- data.frame(names=c('0hr_55mM_0_t1_050412','KClpulse-60min_55mM_20_b4_050412','KClpulse-60min_55mM_40_b4_050412','KClpulse-60min_55mM_60_b4_050412','KClpulse-60min_55mM_120_b4_050412','KClpulse-60min_55mM_240_b4_050412','KClpulse-60min_55mM_360_b4_050412'),timepoints=c(0,20,40,60,120,240,360))
   
   controls1a <- plotCNTLs(c1,CNTLs_c1,sixtyminpulse_55mM_b4)
  
    plotAllPrimers(x=c1,expt= sixtyminpulse_55mM_b4,exptCntl= controls1a,primers=primerSet_c1,fN="biorep2",newWindow=TRUE)     
   
##Experiment 2
##KCl Pulse length 15 min timepoints(min): 0,20,40,60,120,240,360

   fifteenminpulse_55mM_b1 <- data.frame(names=c("0hr_55mM_0_t1_050312", "KClpulse-15min_55mM_20_b1_050312","KClpulse-15min_55mM_40_b1_050312","KClpulse-15min_55mM_60_b1_050312","KClpulse-15min_55mM_120_b1_050312","KClpulse-15min_55mM_240_b1_050312","KClpulse-15min_55mM_360_b1_050312"),timepoints=c(0,20,40,60,120,240,360)) 
   
      	(controls2 <- plotCNTLs(c6,CNTLs_c5and6, fifteenminpulse_55mM_b1))

   plotAllPrimers(x=c6,expt= fifteenminpulse_55mM_b1,exptCntl= controls2,primers=primerSet_c5and6,fN="f",newWindow=TRUE)     
      	 

 
##Experiment 3
##KCl Pulse length 120 min timepoints(min): 0,20,40,60,120,240,360   
  
   twohrpulse_55mM_b2 <- data.frame(names=c("KClpulse-1min_55mM_0min_t1_050412", "KClpulse-120min_55mM_20min_b2_050412","KClpulse-120min_55mM_40min_b2_050412","KClpulse-120min_55mM_60min_b2_050412","KClpulse-120min_55mM_120min_b2_050412","KClpulse-120min_55mM_240min_b2_050412","KClpulse-120min_55mM_360min_b2_050412"),timepoints=c(0,20,40,60,120,240,360)) 
   
        	(controls3 <- plotCNTLs(c1,CNTLs_c1, twohrpulse_55mM_b2))

    plotAllPrimers(x=c1,expt= twohrpulse_55mM_b2,exptCntl= controls3,primers=primerSet_c1,fN="twohrpulse_55mM_050312",newWindow=TRUE)     

   
##Experiment 4
##KCl Pulse length 30 min timepoints(min): 0,20,40,60,120,240,360 
   
  thirtyminpulse_55mM_b2 <- data.frame(names=c("KClpulse-1min_55mM_0min_t1_050412", "KClpulse-30min_55mM_20min_b2_050412","KClpulse-30min_55mM_40min_b2_050412","KClpulse-30min_55mM_60min_b2_050412","KClpulse-30min_55mM_120min_b2_050412","KClpulse-30min_55mM_240min_b2_050412","KClpulse-30min_55mM_360min_b2_050412"),timepoints=c(0,20,40,60,120,240,360)) 
  
        	(controls4 <- plotCNTLs(c1,CNTLs_c1, thirtyminpulse_55mM_b2))
      	
   plotAllPrimers(x=c1,expt= thirtyminpulse_55mM_b2,exptCntl= controls4,primers=primerSet_c1,fN="thirtyminpulse_55mM_050312",newWindow=TRUE)     
      	 

  
  
  
##Experiment 5
##KCl Pulse length 1 min timepoints(min): 0,20,40,60,120,240,360

  oneminpulse_55mM_b1 <- data.frame(names=c("KClpulse-1min_55mM_0min_t1_050412", "KClpulse-1min_55mM_20min_b1_050412","KClpulse-1min_55mM_40min_b1_050412","KClpulse-1min_55mM_60min_b1_050412","KClpulse-1min_55mM_120min_b1_050412","KClpulse-1min_55mM_240min_b1_050412","KClpulse-1min_55mM_360min_b1_050412"),timepoints=c(0,20,40,60,120,240,360)) 
  
        	(controls5 <- plotCNTLs(c1,CNTLs_c1, oneminpulse_55mM_b1))
      	
   plotAllPrimers(x=c1,expt= oneminpulse_55mM_b1,exptCntl= controls5,primers=primerSet_c1,fN="oneminpulse_55mM_050312",newWindow=TRUE)     
      	 

  
##Experiment 6
##KCl Sustained at 3mM timepoints(min): 0,60,180,360

  sustained_3mM_b2 <- data.frame(names=c("KClsustain_10mM_0min_t1_041112", "KClsustain_3mM_60min_b2_041112","KClsustain_3mM_180min_b2_041112","KClsustain_3mM_360min_b2_041112"),timepoints=c(0,60,180,360))
  
          	(controls6 <- plotCNTLs(c1,CNTLs_c1, sustained_3mM_b2 ))
      	
   plotAllPrimers(x=c1,expt= sustained_3mM_b2 ,exptCntl= controls6,primers=primerSet_c1,fN="sustained_3mM_b2_041112",newWindow=TRUE)     
  
##Experiment 7
##KCl Sustained at 10mM timepoints(min): 0,60,180,360

  sustained_10mM_b2 <- data.frame(names=c("KClsustain_10mM_0min_t1_041112", "KClsustain_10mM_60min_b2_041112","KClsustain_10mM_180min_b2_041112","KClsustain_10mM_360min_b2_041112"),timepoints=c(0,60,180,360))
  
            	(controls7 <- plotCNTLs(c1,CNTLs_c1, sustained_10mM_b2 ))
      	
   plotAllPrimers(x=c1,expt= sustained_10mM_b2 ,exptCntl= controls7,primers=primerSet_c1,fN="sustained_10mM_b2_041112",newWindow=TRUE)  
  
  
##Experiment 8
##KCl Sustained at 3mM timepoints(min): 0,20,40,60,120,240,360

  sustained_3mM_b2_moretime <- data.frame(names=c("KClsustain_3mM_0min_t1_051012", "KClsustain_3mM_20min_b2_051012","KClsustain_3mM_40min_b2_051012","KClsustain_3mM_60min_b2_051012","KClsustain_3mM_120min_b2_051012","KClsustain_3mM_240min_b2_051012","KClsustain_3mM_360min_b2_051012"),timepoints=c(0,20,40,60,120,240,360))
  
            	(controls8 <- plotCNTLs(c1,CNTLs_c1, sustained_3mM_b2_moretime ))
      	
   plotAllPrimers(x=c1,expt= sustained_3mM_b2_moretime ,exptCntl= controls8,primers=primerSet_c1,fN="sustained_3mM_b2_moretime_041112",newWindow=TRUE)  
   
##Experiment 9 
##KCl pulse length 5 min at 10mM timepoints(min): 0,20,40,60,120,240,360

  fiveminpulse_10mM_b3 <- data.frame(names=c("KClpulse-5min_30mM_0_t1_011811", "KClpulse-5min_10mM_20_b3_011811","KClpulse-5min_10mM_40_b3_011811","KClpulse-5min_10mM_60_b3_011811","KClpulse-5min_10mM_120_b3_011811","KClpulse-5min_10mM_240_b3_011811","KClpulse-5min_10mM_360_b3_011811"),timepoints=c(0,20,40,60,120,240,360))
  
            	controls9 <- plotCNTLs(c6,CNTLs_c5and6, fiveminpulse_10mM_b3)
      	
   plotAllPrimers(x=c1,expt= fiveminpulse_10mM_b2 ,exptCntl= controls9,primers=primerSet_c1,fN="fiveminpulse_10mM_011811",newWindow=TRUE)  
   
##Experiment 10
##KCl pulse length 5 min at 30mM timepoints(min): 0,20,40,50,120,240,360

  fiveminpulse_30mM_b3 <- data.frame(names=c("KClpulse-5min_30mM_0_t1_011811", "KClpulse-5min_30mM_20_b3_011811","KClpulse-5min_30mM_40_b3_011811","KClpulse-5min_30mM_60_b3_011811","KClpulse-5min_30mM_120_b3_011811","KClpulse-5min_30mM_240_b3_011811","KClpulse-5min_30mM_360_b3_011811"),timepoints=c(0,20,40,60,120,240,360))
  
            	controls10 <- plotCNTLs(c6,CNTLs_c5and6, fiveminpulse_30mM_b3)
      	
   plotAllPrimers(x=c1,expt= fiveminpulse_30mM_b2 ,exptCntl= controls10,primers=primerSet_c1,fN="fiveminpulse_30mM_b2_011812",newWindow=TRUE)  


 ###################

     make_plottable  <- function(x){

  y <- matrix(x,ncol=length(x),nrow=1)
  oddvals <- seq(1,ncol(y),by=2)
  evenvals <- seq(0,ncol(y),by=2)
  df <- data.frame(means=as.numeric(y[,oddvals]),stderrs=as.numeric(y[,evenvals]))
df
}


which(colnames(y)=='Lhx5')
y <- y[,-20]
listo <- apply(y,2,make_plottable)
plot_all_abstracts(plots,'abstract')


plots <- lapply(names(listo), function(x) plot_max_pulsewidths(listo[[x]],x))

 plot_max_pulsewidths <- function(x,t){
      ggplot(x,aes(x=pulse_lengths,means)) + 
           geom_point() + 
           geom_line() +
           scale_x_continuous(name='KCl duration (min)') +
           scale_y_continuous(name='Max Fold Induction') +
           geom_errorbar(aes(ymin=means-stderrs,ymax=means+stderrs),width=2) +
           opts(title = t)
        }  
 
 
 
 ################################
 ########### NEED THIS ##########
 ################################
 
 ##### x needs to be a made matrix of all maxes/stderrs of each desired timepoint for all primers
##### I'll leave my code for how I went about this       
       
        plot_abstract_pulsewidths <- function(x,fN){  

    pulse_lengths <- c(5,15,30,60,120,360)

     make_plottable  <- function(x){

  y <- matrix(x,ncol=length(x),nrow=1)
  oddvals <- seq(1,ncol(y),by=2)
  evenvals <- seq(0,ncol(y),by=2)
  df <- data.frame(means=as.numeric(y[,oddvals]),stderrs=as.numeric(y[,evenvals]))
df
}

listo <- apply(x,2,make_plottable)
 
 plot_max_pulsewidths <- function(x,t){
  ggplot(x,aes(x=pulse_lengths,means)) + 
           geom_point() + 
           geom_line() +
           scale_x_continuous(name='KCl duration (min)') +
           scale_y_continuous(name='Max Fold Induction') +
           geom_errorbar(aes(ymin=means-stderrs,ymax=means+stderrs),width=2) +
           opts(title = t)    
         } 

plots <- lapply(names(listo), function(x) plot_max_pulsewidths(listo[[x]],x))

 plot_all_abstracts(plots,fN)

    }

### I would like to just create a differention concentration function but there are different concentrations present for this particular experiment..probably could just change the aesthetics on the plot function

plot_abstract_5minconcentrations <- function(x,fN){  

### this will take the mess of matrix > list > list of maxes/stderrs into one df with two columns, making ggplot workable 

     make_plottable  <- function(x){

  y <- matrix(x,ncol=length(x),nrow=1)
  oddvals <- seq(1,ncol(y),by=2)
  evenvals <- seq(0,ncol(y),by=2)
  df <- data.frame(means=as.numeric(y[,oddvals]),stderrs=as.numeric(y[,evenvals]))
df
}

### convert mess matrix to workable data frames for ggplot for all primers (i.e. simple data frames with two columns, means/stderrs)
listo <- apply(x,2,make_plottable)


###this will make a simple scatter plot (stderrs included) for each data frame  
 plot_max_concentrations <- function(x,t){

 	    concs <- c(10,30,55)

  ggplot(x,aes(x=concs,means)) + 
           geom_point() + 
           geom_line() +
           scale_x_continuous(name='Conc (mM)') +
           scale_y_continuous(name='Max Fold Induction') +
           geom_errorbar(aes(ymin=means-stderrs,ymax=means+stderrs),width=2) +
           opts(title = t)    
         } 

   #### make list of plots w/ corresponding titles
plots <- lapply(names(listo), function(x) plot_max_concentrations(listo[[x]],x))

### plot_all_abstracts will simply take a list of plots and lay them out on a pdf file 
 plot_all_abstracts(plots,fN)

    }



##### like mentioned in above comments 

plot_abstract_sustain_concentrations <- function(x,fN){  

### this will take the mess of matrix > list > list of maxes/stderrs into one df with two columns, making ggplot workable 

     make_plottable  <- function(x){

  y <- matrix(x,ncol=length(x),nrow=1)
  oddvals <- seq(1,ncol(y),by=2)
  evenvals <- seq(0,ncol(y),by=2)
  df <- data.frame(means=as.numeric(y[,oddvals]),stderrs=as.numeric(y[,evenvals]))
df
}

### convert mess matrix to workable data frames for ggplot for all primers (i.e. simple data frames with two columns, means/stderrs)
listo <- apply(x,2,make_plottable)

 	    concs_sustained <- c(3,10,30,55)

###this will make a simple scatter plot (stderrs included) for each data frame  
 plot_max_sustained_concentrations <- function(x,t){
 	    concs_sustained <- c(3,10,30,55)

  ggplot(x,aes(x=concs_sustained,means)) + 
           geom_point() + 
           geom_line() +
           scale_x_continuous(name='Conc (mM)') +
           scale_y_continuous(name='Max Fold Induction') +
           geom_errorbar(aes(ymin=means-stderrs,ymax=means+stderrs),width=2) +
           opts(title = t)    
         } 

   #### make list of plots w/ corresponding titles
plots <- lapply(names(listo), function(x) plot_max_sustained_concentrations(listo[[x]],x))

### plot_all_abstracts will simply take a list of plots and lay them out on a pdf file 
 plot_all_abstracts(plots,fN)

    }


#########################
############


      maxes_sustained_3mM <- maxes(all_sustain_3mM)
      maxes_sustained_3mM$conc <- 3
      
      maxes_sustained_10mM <- maxes(all_sustain_10mM)
      maxes_sustained_10mM$conc <- 10
      
      maxes_sustained_30mM <- maxes(c5_30mMsustained_b1)
      maxes_sustained_30mM$conc <- 30
      
      maxes_sustained_55mM <- maxes(all_sustain_55mM)
      maxes_sustained_55mM$conc <- 55
      
      conc_sustain_maxes <- rbind(maxes_sustained_3mM,maxes_sustained_10mM,maxes_sustained_30mM,maxes_sustained_55mM)
      conc_sustain_maxes <- conc_sustain_maxes[order(conc_sustain_maxes$primer),]
      
      plotAbs_conc(primes_sustain_3mM,conc_sustain_maxes,'sustain_diff_concs')
           
      

Nr4a1$sde <- Nr4a1$sde / max(Nr4a1$mean)
Nr4a1$mean <- Nr4a1$mean / max(Nr4a1$mean)

Nr4a2$sde <- Nr4a2$sde / max(Nr4a2$mean)
Nr4a2$mean <- Nr4a2$mean / max(Nr4a2$mean)

Npas4$sde <- Npas4$sde / max(Npas4$mean)
Npas4$mean <- Npas4$mean / max(Npas4$mean)

Junb$sde <- Junb$sde / max(Junb$mean)
Junb$mean <- Junb$mean / max(Junb$mean)

Gadd45b$sde <- Gadd45b$sde / max(Gadd45b$mean)
Gadd45b$mean <- Gadd45b$mean / max(Gadd45b$mean)

Maff$sde <- Maff$sde / max(Maff$mean)
Maff$mean <- Maff$mean / max(Maff$mean)

iegs <- rbind(Nr4a1,Nr4a2,Npas4,Junb,Gadd45b,Maff)


one <- find_means(all_1min)

care1 <- one[one$primer %in% c('Nr4a1','Npas4'),]
caree1 <- one[one$primer %in% c('Nr4a1','Nr4a2'),]
care1a <- one[one$primer %in% c('Egr1','Nr4a2'),]
care1b <- one[one$primer %in% c('Egr2','Nr4a2'),]

five <- find_means(all_5min)

care5 <- five[five$primer %in% c('Nr4a1','Npas4'),]
caree5 <- five[five$primer %in% c('Nr4a1','Nr4a2'),]
care5a <- five[five$primer %in% c('Egr1','Nr4a2'),]
care5b <- five[five$primer %in% c('Egr2','Nr4a2'),]

sust <- find_means(all_sustain_55mM)

caresus <- sust[sust$primer %in% c('Nr4a1','Npas4'),]
careesus <- sust[sust$primer %in% c('Nr4a1','Nr4a2'),]
caresusa <- sust[sust$primer %in% c('Egr1','Nr4a2'),]
caresusb <- sust[sust$primer %in% c('Egr2','Nr4a2'),]

pdf('1min_egr1_compare.pdf',onefile=T)
 	    g <- ggplot(care1a,aes(x=timepoints,y= mean,ymax = mean + sde, ymin = mean  - sde,colour=primer)) + 
           geom_errorbar(width=2) + 
           geom_point() +
           geom_line() +
           xlab('KCl duration (min)') +
           ylab('Fold Induction') +          
           opts(title = '1 min KCl')    
        print(g)
dev.off()

 
 
 ################ MULTI PLOTS
 ##############################
 
 
 
pdf('PlotAA.pdf', onefile=T)
g <- ggplot(sus,aes(x=timepoints,y=mean, ymax = mean + sde, ymin = mean - sde,colour=primer)) +
   geom_point() + 
   geom_line() + 
   xlab('Time (min)') + 
   ylab('Fold Induction') + 
   geom_errorbar(width=2) +
   opts(title='6hr KCl')
print(g)  
dev.off()
   
lapply(names(your_list), function(x) plot2(your_list[[x]], x))


listo   <-  for (i in unique(fivemin$primer)){
  	g <- ggplot(fivemin[fivemin$primer == i,], aes(x = timepoints,y = mean, ymin = mean -sde, ymax = mean + sde)) +
  	geom_errorbar(width=2) + geom_point() + geom_line() +
  	facet_wrap(~ primer, ncol=5,nrow=1) +
  	xlab('Time (min)') + ylab('Fold Induction') 
  	print(g) 
  
  
         plot_multiple <- function(x){         
        x <- ddply(x, .(primer, timepoints), summarize, mean = mean(foldInduction, na.rm = TRUE), sde = sqrt(var(foldInduction, na.rm = TRUE)/length(foldInduction))) 
        x}
        
 	    ggplot(x,aes(x=timepoints,y= mean,ymax = mean + sde, ymin = mean  - sde)) + 
           geom_errorbar(width=2) + 
           geom_point() +
           geom_line() +
           xlab('Time (min)') +
           ylab('Fold Induction') +          
           opts(title = prime)    
           
           }  

#############
           
           
           ############### MORE RANDOM PLOTS FOR JESSE
           
           fos_new <- rbind(fos,fos_5)
colnames(fos_new) = c('KCl_Duration','timepoints','mean','sde')
fos_new <- fos_new[fos_new$timepoints != 180,]
egr4_new <- rbind(egr4,egr4_5)
colnames(egr4_new) = c('KCl_Duration','timepoints','mean','sde')
egr4_new <- egr4_new[egr4_new$timepoints != 180,]
npas4_new <- rbind(npas4,npas4_5)
colnames(npas4_new) = c('KCl_Duration','timepoints','mean','sde')
npas4_new <- npas4_new[npas4_new$timepoints != 180,]
junb_new <- rbind(junb,junb_5)
colnames(junb_new) = c('KCl_Duration','timepoints','mean','sde')
junb_new <- junb_new[junb_new$timepoints != 180,]
nr4a1_new <- rbind(nr4a1,nr4a1_5)
colnames(nr4a1_new) = c('KCl_Duration','timepoints','mean','sde')
nr4a1_new <- nr4a1_new[nr4a1_new$timepoints != 180,]
maff_new <- rbind(maff,maff_5)
colnames(maff_new) = c('KCl_Duration','timepoints','mean','sde')
maff_new <- maff_new[maff_new$timepoints != 180,]
egr1_new <- rbind(egr1,egr1_5)
colnames(egr1_new) = c('KCl_Duration','timepoints','mean','sde')
egr1_new <- egr1_new[egr1_new$timepoints != 180,]
arc_new <- rbind(arc,arc_5)
colnames(arc_new) = c('KCl_Duration','timepoints','mean','sde')
arc_new <- arc_new[arc_new$timepoints != 180,]
egr2_new <- rbind(egr2,egr2_5)
colnames(egr2_new) = c('KCl_Duration','timepoints','mean','sde')
egr2_new <- egr2_new[egr2_new$timepoints != 180,]
klf4_new <- rbind(klf4,klf4_5)
colnames(klf4_new) = c('KCl_Duration','timepoints','mean','sde')
klf4_new <- klf4_new[klf4_new$timepoints != 180,]


serpinb2_new <- rbind(serpinb2,serpinb2_5)
colnames(serpinb2_new) = c('KCl_Duration','timepoints','mean','sde')
gprc5a_new <- rbind(gprc5a,gprc5a_5)
colnames(gprc5a_new) = c('KCl_Duration','timepoints','mean','sde')
acan_new <- rbind(acan,acan_5)
colnames(acan_new) = c('KCl_Duration','timepoints','mean','sde')
acan_new <- acan_new[acan_new$timepoints != 180,]
pthlh_new <- rbind(pthlh,pthlh_5)
colnames(pthlh_new) = c('KCl_Duration','timepoints','mean','sde')
pthlh_new <- pthlh_new[pthlh_new$timepoints != 180,]
npffr2_new <- rbind(npffr2,npffr2_5)
colnames(npffr2_new) = c('KCl_Duration','timepoints','mean','sde')
npffr2_new[npffr2_new$timepoints != 180,]


rasd1_new <- rbind(rasd1,rasd1_5)
colnames(rasd1_new) = c('KCl_Duration','timepoints','mean','sde')
rasd1_new <- rasd1_new[rasd1_new$timepoints != 180,]
pcsk1_new <- rbind(pcsk1,pcsk1_5)
colnames(pcsk1_new) = c('KCl_Duration','timepoints','mean','sde')
pcsk1_new <- pcsk1_new[pcsk1_new$timepoints != 180,]
pax1_new <- rbind(pax1,pax1_5)
colnames(pax1_new) = c('KCl_Duration','timepoints','mean','sde')
pax1_new <- pax1_new[pax1_new$timepoints != 180,]
bdnf_new <- rbind(bdnf,bdnf_5)
colnames(bdnf_new) = c('KCl_Duration','timepoints','mean','sde')
bdnf_new <- bdnf_new[bdnf_new$timepoints != 180,]


pdf('fos.pdf',onefile=T)
 	    g <- ggplot(fos_new,aes(x=timepoints,y= mean,ymax = mean + sde, ymin = mean  - sde,colour=Exp)) + 
 	       scale_colour_manual(values=c('red','black')) +
           geom_errorbar(width=2) + 
           geom_point() +
           geom_line() +
           xlab('Time (min)') +
           ylab('Fold Induction') +          
           opts(title = '5 min v 6 hr KCl treatment (Fos)')    
        print(g)
dev.off()
 
 
 