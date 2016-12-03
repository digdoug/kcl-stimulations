# some generic stuff

options(stringsAsFactors = F)


# read in results from a specified file (with extra headers removed)
read_results <- function(fN,IEGs,LRGs,CNTLs){ 
	
   # troubleshooting
   # fN <- "/Volumes/Neurobio/GreenbergLab/Joe J Ling/IEG_project/Fluidigm/results/1361523246_2011_1_11_Taqman_Fluidigm_Assays_01142011_for_R.tsv"
   # fN <- "/Volumes/Neurobio/GreenbergLab/Joe J Ling/IEG_project/Fluidigm/results/1361523246_2011_1_11_Taqman_Fluidigm_Assays_01142011_for_R.tsv"

   head(x <- read.delim(file=fN,as.is=TRUE))

   x$sampleNum <- sapply(x$ID,function(x){strsplit(x,split="-")[[1]][1]})
   x$assayNum <- sapply(x$ID,function(x){strsplit(x,split="-")[[1]][2]})

   label_by_class <- function(myNames,IEGs,LRGs,CNTLs){
	 output <- myNames
	 output[which(!is.na(match(myNames,IEGs)))] <- "IEG"
	 output[which(!is.na(match(myNames,LRGs)))] <- "LRG"
	 output[which(!is.na(match(myNames,CNTLs)))] <- "CNTL"
	 output
   }
   head(x$gene_class <- label_by_class(x$Name.1,IEGs,LRGs,CNTLs))

   # reorder sensibly
   head(x <- x[order(x$assayNum,x$sampleNum),])

   # rename some columns sensibly
   colnames(x)[which(colnames(x) == "Name.1")] <- "primer"
   colnames(x)[which(colnames(x) == "Name")] <- "sampleName"
   
   # annotate the standard curve
   # date doesn't matter, can change this in the future, for now i'll just set an arbitrary number
   x$Conc <- 0.05
   x$Conc[which(x$sampleName == "Standard-5050mixof1and6hrKCl_1_NA_NA_123456")] <- 1
   x$Conc[which(x$sampleName == "Standard-5050mixof1and6hrKCl_0.2_NA_NA_123456")] <- 0.2   
   x$Conc[which(x$sampleName == "Standard-5050mixof1and6hrKCl_0.1_NA_NA_123456")] <- 0.1
   x$Conc[which(x$sampleName == "Standard-5050mixof1and6hrKCl_0.01_NA_NA_123456")] <- 0.01
   x$Conc[which(x$sampleName == "Standard-5050mixof1and6hrKCl_0.001_NA_NA_123456")] <- 0.001
   x$Conc[which(x$sampleName == "Standard-5050mixof1and6hrKCl_0.0001_NA_NA_123456")] <- 0.0001
   x$Conc[which(x$sampleName == "Standard-5050mixof1and6hrKCl_0.00001_NA_NA_123456")] <- 0.00001

   # label the standards
   x$SampleClass <- NA
   x$SampleClass[grep("Standard-5050mixof1and6hrKCl",x$sampleName)] <- "standard"
   
   
   
   # annotate biological and technical replicates
   x$sampleTechRep <- NA
   x$sampleTechRep[which(!is.na(match(x$sampleNum,c("S01","S13","S25","S37","S49"))))] <- 1 # this is sample technical replicates (look them up)
   x$sampleTechRep[which(!is.na(match(x$sampleNum,c("S07","S19","S31","S43","S55"))))] <- 2 # there could also be assay technical reps, ignored here
   
   primer_names <- unique(x$primer)
   primer_names <- primer_names[which(primer_names!="water")] 
   
   # if there isn't at least two points to draw a line, it results in "error in int_abline..."
   head(std_curves <- aperm(sapply(primer_names,plotStdCurve,x)))

   # add a concentration column based on standard curve
   head(x$realConc <- apply(x,1,getConc,std_curves))
   x$realConc <- as.numeric(x$realConc)

   #Add new columns containing all these things, which will make further analysis easier 
     x$exptname <- NA
     x$exptname <- sapply(x$sampleName,function(x){strsplit(x,split='_')[[1]][1]})
   	 x$concentrate <- NA
   	 x$concentrate <- sapply(x$sampleName,function(x){strsplit(x,split='_')[[1]][2]})
  	 x$timepoints <- NA
  	 x$timepoints <- sapply(x$sampleName,function(x){strsplit(x,split='_')[[1]][3]})
   	 x$replicate <- NA
   	 x$replicate <- sapply(x$sampleName,function(x){strsplit(x,split='_')[[1]][4]})
   	 x$day <- NA
   	 x$day <- sapply(x$sampleName,function(x){strsplit(x,split='_')[[1]][5]})


   # return the data
   x
   
   
}



# plot replicability across all primers and samples on a chip
# the Quad version assumes there are at least some 2x sample and 2x assay replicates 
plotReplicabilityQuad <- function(x,myWhich="both"){
#x <- p1   	
      a <- x$Value 
      b <- x$sampleNum
      c <- x$assayNum
      d <- x$primer
      e <- x$sampleName

      head(y <- data.frame(Value=a,sampleNum=b,assayNum=c,primer=d,sampleName=e))
      head(z <- y[order(y$sampleName,y$primer,y$assayNum),],n=20)
      z$Value[which(z$Value == 999)] <- NA
      
      
      # remove any samples /assays for which we don't have quadruplicates
      primer_sample <- paste(z$primer,z$sampleName)
      unique_primer_sample <- unique(primer_sample)
      test <- function(id,x){
         sum(!is.na(match(x,id)))	
      }
      passing <- unique_primer_sample[which(sapply(unique_primer_sample,test,primer_sample) == 4)]
      dim(z <- z[which(!is.na(match(primer_sample,passing))),])
      
      # this is really a BS way of doing this that may not work on all chips; assumes quadruplicate 2x sample 2x assay
      head(S1A1 <- z[seq(from=1,to=length(z$Value),by=4),]$Value)
      head(S2A1 <- z[seq(from=2,to=length(z$Value),by=4),]$Value)
      head(S1A2 <- z[seq(from=3,to=length(z$Value),by=4),]$Value)
      head(S2A2 <- z[seq(from=4,to=length(z$Value),by=4),]$Value)
      myMin <- min(z$Value,na.rm=T)
      myMax <- max(z$Value,na.rm=T)
      if(myWhich == "both"){
         plot(c(S1A1,S2A1),c(S1A2,S2A2),main="Fluidigm replicability",frame.plot=F,xlab="replicate 1",ylab="replicate 2",xlim=c(myMin,myMax),ylim=c(myMin,myMax),cex=0.1)     
         points(c(S1A1,S1A2),c(S2A1,S2A2),col="red",cex=0.1) 
         legend("topleft",legend=c("assay reps","sample reps"),text.col=c("black","red"),bty="n",cex=2)
      } else if(myWhich == "sample") {
      	 plot(c(S1A1,S1A2),c(S2A1,S2A2),main="Sample replicability",frame.plot=F,xlab="replicate 1",ylab="replicate 2",xlim=c(myMin,myMax),ylim=c(myMin,myMax),cex=0.1)
      } else if (myWhich == "assay"){
         plot(c(S1A1,S2A1),c(S1A2,S2A2),main="Assay replicability",frame.plot=F,xlab="replicate 1",ylab="replicate 2",xlim=c(myMin,myMax),ylim=c(myMin,myMax),cex=0.1)
      }
}







# plot replicability for a specific sample
plotReplicabilityBySample <- function(sampleName,x){

   	
      # look at replicability
      a <- x[which(x$sampleName==sampleName),]$Value 
      b <- x[which(x$sampleName==sampleName),]$sampleNum
      c <- x[which(x$sampleName==sampleName),]$assayNum
      d <- x[which(x$sampleName==sampleName),]$primer

      head(y <- data.frame(Value=a,sampleNum=b,assayNum=c,primer=d))
      head(z <- y[order(y$primer,y$assayNum),],n=20)
      z$Value[which(z$Value == 999)] <- NA
      
      head(S1A1 <- z[seq(from=1,to=length(z$Value),by=4),]$Value)
      head(S2A1 <- z[seq(from=2,to=length(z$Value),by=4),]$Value)
      head(S1A2 <- z[seq(from=3,to=length(z$Value),by=4),]$Value)
      head(S2A2 <- z[seq(from=4,to=length(z$Value),by=4),]$Value)
      myMin <- min(z$Value,na.rm=T)
      myMax <- max(z$Value,na.rm=T)
      plot(c(S1A1,S2A1),c(S1A2,S2A2),main="Fluidigm replicability",frame.plot=F,xlab="replicate 1",ylab="replicate 2",xlim=c(myMin,myMax),ylim=c(myMin,myMax))     
      points(c(S1A1,S1A2),c(S2A1,S2A2),col="red") 
      legend("topleft",legend=c("assay reps","sample reps"),text.col=c("black","red"),bty="n",cex=2)
        
	  assayR2 <- cor.test(c(S1A1,S2A1),c(S1A2,S2A2))$estimate
	  sampleR2 <- cor.test(c(S1A1,S1A2),c(S2A1,S2A2))$estimate
	  mean(assayR2,sampleR2)
}

# plot replicability for a specific primer
plotReplicabilityByPrimer <- function(primer,x,removeRT=T){


#primer <- "Fos"
#x <- p2

   	  if(removeRT == TRUE){x <- x[grep("RT-",x$sampleName,invert=T),]} # remove RT- samples
   	
      # look at replicability
      a <- x[which(x$primer==primer),]$Value 
      b <- x[which(x$primer==primer),]$sampleNum
      c <- x[which(x$primer==primer),]$assayNum
      d <- x[which(x$primer==primer),]$sampleName

      head(y <- data.frame(Value=a,sampleNum=b,assayNum=c,sample=d))
      head(z <- y[order(y$sample,y$assayNum),],n=20)
      z$Value[which(z$Value == 999)] <- NA
     
      head(S1A1 <- z[seq(from=1,to=length(z$Value),by=4),]$Value)
      head(S2A1 <- z[seq(from=2,to=length(z$Value),by=4),]$Value)
      head(S1A2 <- z[seq(from=3,to=length(z$Value),by=4),]$Value)
      head(S2A2 <- z[seq(from=4,to=length(z$Value),by=4),]$Value)

      myMin <- min(z$Value,na.rm=T)
      myMax <- max(z$Value,na.rm=T)
      plot(c(S1A1,S2A1),c(S1A2,S2A2),main=paste("Fluidigm replicability (",primer,")",sep=""),frame.plot=F,xlab="replicate 1",ylab="replicate 2",xlim=c(myMin,myMax),ylim=c(myMin,myMax))     
      points(c(S1A1,S1A2),c(S2A1,S2A2),col="red") 
      legend("topleft",legend=c("assay reps","sample reps"),text.col=c("black","red"),bty="n",cex=2)  

	  assayR2 <- cor.test(c(S1A1,S2A1),c(S1A2,S2A2))$estimate
	  sampleR2 <- cor.test(c(S1A1,S1A2),c(S2A1,S2A2))$estimate
	  mean(assayR2,sampleR2)
	  # to spit out the table: z
}


# for any given primer just plot the std curves by replicate
plotStdCurveByRep <- function(primer,x,printData=F){
   	  # sample replicate 1
      toPl <- x[which(x$primer == primer & x$SampleClass == "standard" & x$sampleTechRep == 1 & x$Value != 999 & x$Call == "Pass"),]
      if(dim(toPl)[1] > 1){
         plot(x=log10(toPl$Conc),y=toPl$Value,ylim=c(0,40),frame.plot=F,main=primer,xlab="standard dilution",ylab="Ct")
         myX <- log10(toPl$Conc); myY <- toPl$Value
         abline(myLM <- lm(myY ~ myX))
         r2.1 <- summary(myLM)$r.squared
         myRange.1 <- range(myY) 
         slope1 <- myLM$coefficients[2]
         int1 <- myLM$coefficients[1]
      } 
      
      # sample replicate 2
      toPl <- x[which(x$primer == primer & x$SampleClass == "standard" & x$sampleTechRep == 2 & x$Value != 999 & x$Call == "Pass"),]
      if(dim(toPl)[1] > 1){
         points(x=log10(toPl$Conc),y=toPl$Value,ylim=c(0,40),col="red")
         myX <- log10(toPl$Conc); myY <- toPl$Value
         abline(myLM <- lm(myY ~ myX),col="red")       
      	 r2.2 <- summary(myLM)$r.squared
         myRange.2 <- range(myY) 
         slope2 <- myLM$coefficients[2]
         int2 <- myLM$coefficients[1]
      } 
}
      

# for any given primer plot the standard curve values
# and spit out enough info to evaluate the standard curve later  
   plotStdCurve <- function(primer,x,printData=F){      

   	# aggregated analysis of both replicates
    toPl <- x[which(x$primer == primer & x$SampleClass == "standard" & x$Value != 999),]

    if(length(unique(toPl$sampleName)) > 1){
      if(dim(toPl)[1] > 1){
         plot(x=log10(toPl$Conc),y=toPl$Value,ylim=c(0,40),frame.plot=F,main=primer,xlab="standard dilution",ylab="Ct")
         myX <- log10(toPl$Conc); myY <- toPl$Value
         abline(myLM <- lm(myY ~ myX))
         r2.a <- summary(myLM)$r.squared
         myRange.a <- range(myY) 
         slopea <- myLM$coefficients[2]
         inta <- myLM$coefficients[1]
      } else {
      	 r2.a <- NA; myRange.a <- NA; slopea <- NA; inta <- NA
      }
      
      
      # processing both/either replicate
      # spit out a list with the results of the standard curve analysis
      if(printData == F){
      	list(r2=r2.a,meanL=myRange.a[1],meanH=myRange.a[2],slope=slopea,int=inta)
      } else {
      	x[which(x$primer == primer & x$SampleClass == "standard" & x$Value != 999 & x$Call == "Pass"),]
      }
   } else {
      if(printData == F){
      	list(r2=NA,meanL=NA,meanH=NA,slope=NA,int=NA)
      } else {
      	x[which(x$primer == primer & x$SampleClass == "standard" & x$Value != 999 & x$Call == "Pass"),]
      }   	        		
   }
  }
   
   

# for each sample, use the corresponding appropriate standard curve to get a concentration value
   getConc <- function(x,std_curves){
   	    # label the data
   		Ct <- as.numeric(x[7])
   		primerName <- x[5]
   		
   		# label/extract the std curve info
		minRange <- std_curves[which(rownames(std_curves) == primerName),]$meanL
		maxRange <- std_curves[which(rownames(std_curves) == primerName),]$meanH
		r2 <- std_curves[which(rownames(std_curves) == primerName),]$r2
		slope <- as.numeric(std_curves[which(rownames(std_curves) == primerName),]$slope)
		int <- as.numeric(std_curves[which(rownames(std_curves) == primerName),]$int)

		if(all(!is.na(c(minRange,maxRange,r2,slope,int)))){
			#if(Ct > minRange & Ct < maxRange & r2 > 0.2){
				conc <- 10^((Ct - int) / slope)
			#} else {
			#	conc <- NA	
			#}
		} else { conc <- NA}
		conc
   }
   
   
   
   
# across a particular experiment, plot the controls as a consistency check
   # for a timecourse
  
     getCNTLs_mod <- function(x,CNTLs,exptName){

   	   	   x <- x[which(!is.na(match(x$primer,CNTLs))),]           
   	   	   x <- x[which(!is.na(match(x$sampleName,exptName$names))),]
           x <- x[x$Value != 999,]
           
           agg <- c('sampleName','primer','exptname','concentrate','timepoints','replicate','day','realConc')
           x <- x[,agg]
           x <- aggregate(realConc ~ .,mean,data=x)

  		  means <- tapply(as.numeric(x$realConc),factor(as.numeric(x$timepoints)),mean,na.rm=T)
means
}

##normalize concentration value according to experiment and its control timepoin
  #### make individual experiments for this one

      normalize_con <- function(primers,x,exptName,c){


   	   	   # limit x to only the relevant rows; then set the x-values for plotting (designed for normal six hour timecourse plots)
   	   	   x <- x[which(!is.na(match(x$primer,primers))),]         	   	        
   	   	   x <- x[which(!is.na(match(x$sampleName,exptName$names))),]           
           x <- x[x$Value != 999,]
   
         
           ### If it is the case that technical replicates exist combine the technical replciates for all timepoints
           agg <- c('sampleName','primer','exptname','concentrate','timepoints','replicate','day','realConc')
           x <- x[,agg]
           x <- aggregate(realConc ~ .,mean,data=x)                                 	         	     
           class(x$timepoints) = 'numeric'

          ##divide by controls before filtering out genes
          ## now divide everythig by their controls/timepoints
          
          ### create a column which contains number which are the value
          x$control_match <- match(x$timepoints,names(c))
          
            for(i in 1:dim(x)[1]){   
          x$realConc[i] <- x$realConc[i] / c[[x$control_match[i]]]
          }  

        ## there are some primers which do not have any 0 hrs which passed the Ct. This filter out entire replicates for genes which do not
          
            
          filt <- ddply(x, .(primer,day), summarise, count=sum(timepoints==0)) # this will tell you all primers that have a 0 hr time point by giveing a 1 in the count column
             
             
        if (any(filt$count == 0)) { # this was the case once so I implemented this if else part
          
        filt <- filt[filt$count == 0,]
        include <-!(x$primer%in%filt$primer)&(x$day%in%filt$day) # all primers that have 0 hrs
        x <- x[include,] 
        ### for any given replicate, divide each timepoint by its zero hour 
        x <- ddply(x, .(primer),transform, foldInduction=realConc/realConc[timepoints==0])
         
        } else {
       	x <- ddply(x, .(primer), transform, foldInduction=realConc/realConc[timepoints==0]) 
        }
 
   x[,-9]

   } 	
  
  

   # for a timecourse with muultiple conditions
   plotCNTLsMulti <- function(x,CNTLs,exptName){
   
   			# troubleshooting

   	 		# create an average control profile
   	   		x <- x[which(!is.na(match(x$primer,CNTLs))),]

			myColors <- unique(exptName[,3])
			myLabels <- unique(exptName[,4])  	   		
   			x$myX <- NA
   			x$myColor <- NA
   			for(i in 1:dim(exptName)[1]){
   			   x$myX[which(exptName[i,1]==x$sampleName)] <- exptName[i,2]
   			   x$myColor[which(exptName[i,1]==x$sampleName)] <- exptName[i,3]
   			}

   			use <- which(!is.na(x$myX))

   			means <- by(as.numeric(x$realConc[use]),list(factor(x$myX[use]),factor(x$myColor[use])),mean,na.rm=T)[,]
		
   			if(length(dim(means)) > 1){
   			   plot(rownames(means),means[,1],xlab="Time (min)",ylab="Expression",frame.plot=F,main="control genes",ylim=c(0,1.5*max(means)))
   			   for(i in (match(myColors,colnames(means)))){
   			      points(rownames(means),means[,i],col=myColors[i])
   			   }
   			} else {
   			   names(means) <- myLabels[match(myColors,names(means))]
   			   barplot(means,xlab="Time (min)",ylab="Expression",main="control genes",ylim=c(0,1.5*max(means)))
				
   			}

   			means
   }   
 
 
   #####  get maximum values 
   # 
   
      maxes <- function(x){
   x <- ddply(x, .(primer,timepoints), summarize, mean = mean(foldInduction, na.rm = TRUE), sde = sqrt(var(foldInduction, na.rm = TRUE)/length(foldInduction))) 
   x <- ddply(x, .(primer), function(x) x[x$mean == max(x$mean),-2])
   x
}
   
##### This will actually give you means/stderrs 
   
   means <- function(x){
   	
   	        x <- ddply(x, .(primer, timepoints), summarize, mean = mean(foldInduction, na.rm = TRUE), sde = sqrt(var(foldInduction, na.rm = TRUE)/length(foldInduction))) 
   x
   }
 	   
 	   
## plot
 	   
 	 plot_sim <- function(x){
 	 	
 	    ggplot(Pcsk1,aes(x=timepoints,y= mean,ymin=mean-sde,ymax=mean+sde)) + 
           geom_errorbar(,width=2) + 
           geom_point(colour='black') +
           # stat_smooth(method='loess',colour='green') +
           geom_line(colour='black') +
           # geom_ribbon(aes(ymin=0,ymax=mean,colour='orange')) +
           xlab('Time (min)') +
           ylab('Fold Induction') +          
           opts(title = 'test')  
}

##
#  NORMAL PLOT


         plot_normal <- function(prime, x){   
        x <- ddply(x, .(primer, timepoints), summarize, mean = mean(foldInduction, na.rm = TRUE), sde = sqrt(var(foldInduction, na.rm = TRUE)/length(foldInduction))) 
        x <- x[x$primer == prime,]  			
 	    ggplot(x,aes(x=timepoints,y= mean,ymax = mean + sde, ymin = mean  - sde)) + 
           geom_errorbar(width=2) + 
           geom_point() +
           geom_line() +
           xlab('Time (min)') +
           ylab('Fold Induction') +          
           opts(title = prime)    
           
           }

     plot_2pulse <- function(prime, x){   
        x <- ddply(x, .(primer, timepoints), summarize, mean = mean(foldInduction, na.rm = TRUE), sde = sqrt(var(foldInduction, na.rm = TRUE)/length(foldInduction))) 
        x <- x[x$primer == prime,]  			
 	    ggplot(x,aes(x=timepoints,y= mean,ymax = mean + sde, ymin = mean  - sde)) + 
           geom_errorbar(width=2) + 
           geom_point() +
           geom_line() +
           geom_vline(xintercept = c(0,120),colour='orange') +
           xlab('Time (min)') +
           ylab('Fold Induction') +          
           opts(title = prime)    
           
           }

### If you want one plot w/ any number of genes for a fold induction v time plot, enter in any array of primers in the 'prime' argument


         plot_mult <- function(prime, x){   

        x <- x[x$primer %in% prime,]  
        		
 	    ggplot(x,aes(x=timepoints,y= mean,ymax = mean + sde, ymin = mean  - sde,colour=type)) + 
           geom_errorbar(width=2) + 
           geom_point() +
           geom_line() +
           xlab('Time (min)') +
           ylab('Fold Induction') +          
           opts(title = prime)    
           
           }
           
           plot_big_boy('Egr1',all_1min,all_5min,all_15min,all_30min,all_60min,all_120min,all_sus)
           
           plot_big_boy <- function(prime,ex1,ex2,ex3,ex4,ex5,ex6,ex7){   
      
      # I think one can optimize this funciton fairly simply by making it function(prime,...) but I'm not sure
      
        ex1 <- ex1[ex1$primer %in% prime,]
        ex1$primer <- '1 min'
        ex2 <- ex2[ex2$primer %in% prime,]
        ex2$primer <- '5 min'
        ex3 <- ex3[ex3$primer %in% prime,]
        ex3$primer <- '15 min'
        ex4 <- ex4[ex4$primer %in% prime,]
        ex4$primer <- '30 min'
        ex5 <- ex5[ex5$primer %in% prime,]
        ex5$primer <- '60 min'
        ex6 <- ex6[ex6$primer %in% prime,]
        ex6$primer <- '120 min'
        ex7 <- ex7[ex7$primer %in% prime,] 
        ex7$primer <- '360 min'
        
        x <- rbind(ex1,ex2,ex3,ex4,ex5,ex6,ex7)
        		
        		
 	    ggplot(x,aes(x=timepoints,y= mean,ymax = mean + sde, ymin = mean  - sde,colour=primer)) + 
           geom_errorbar(width=2) + 
           geom_point() +
           geom_line() +
           xlab('Time (min)') +
           ylab('Fold Induction') +          
           opts(title = prime)    
           
           }

## MAX PULSE


         plot_max_pulse <- function(prime, x){   
      
        x <- x[x$primer == prime,]  			
 	    ggplot(x,aes(x=kcl_duration,y= mean,ymax = mean + sde, ymin = mean  - sde)) + 
           geom_errorbar(width=2) + 
           geom_point() +
           geom_line() +
           xlab('KCl duration (min)') +
           ylab('Max Fold Induction') +          
           opts(title = prime)    
           
           }

### MAX CONC ####
####

         plot_max_conc <- function(prime, x){   
      
        x <- x[x$primer == prime,]  			
 	    ggplot(x,aes(x=conc,y= mean,ymax = mean + sde, ymin = mean  - sde)) + 
           geom_errorbar(width=2) + 
           geom_point() +
           geom_line() +
           xlab('KCl conc (mM)') +
           ylab('Max Fold Induction') +          
           opts(title = prime)    
           
           }


###          
## Using ggplot plotall primers for a given experiment, across several pdfs. Whoever runs this will need to manually stretch the quartz page in order for it LOOK good
      
    plotAllPrimers_gg <- function(primes,x,fN){

           
 plot_list <- lapply(primes,plot_normal,x)


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




   plotAll2pulse <- function(primes,x,fN){

           
 plot_list <- lapply(primes,plot_2pulse,x)


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






    plotchx_gg <- function(primes,x,fN){

           
 plot_list <- lapply(primes,plot_mult,x)


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
  
  
   multiplot(plot_list[1:6],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_1-6.pdf",sep=""),type="pdf")
   multiplot(plot_list[7:12],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_7-12.pdf",sep=""),type="pdf")
   multiplot(plot_list[13:18],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_13-18.pdf",sep=""),type="pdf")
   multiplot(plot_list[19:24],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_19-24.pdf",sep=""),type="pdf")
   multiplot(plot_list[25:30],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_25-30.pdf",sep=""),type="pdf")
   multiplot(plot_list[31:36],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_31-36.pdf",sep=""),type="pdf")
   multiplot(plot_list[37:42],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_37-42.pdf",sep=""),type="pdf")
   multiplot(plot_list[43:length(plot_list)],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_43-end.pdf",sep=""),type="pdf")

   }


#################
#############
#################


  plotAllrnas <- function(primes,x,fN){

           
 plot_list <- lapply(primes,rna_plot,x)


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
           
#### plot some abstractness

    plotAbs_pulse <- function(primes,x,fN){
           
 plot_list <- lapply(primes,plot_max_pulse,x)


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
  
   }


#### plot max FIs for the same experiment but each x value will be a different concentration

    plotAbs_conc <- function(primes,x,fN){
           
 plot_list <- lapply(primes,plot_max_conc,x)

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
   

    plotBIG <- function(primes,ex1,ex2,ex3,ex4,ex5,ex6,ex7,fN){
           
 plot_list <- lapply(primes,plot_big_boy,ex1,ex2,ex3,ex4,ex5,ex6,ex7)

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
  
  
   multiplot(plot_list[1:6],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_1-6.pdf",sep=""),type="pdf")
   multiplot(plot_list[7:12],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_7-12.pdf",sep=""),type="pdf")
   multiplot(plot_list[13:18],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_13-18.pdf",sep=""),type="pdf")
   multiplot(plot_list[19:24],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_19-24.pdf",sep=""),type="pdf")
   multiplot(plot_list[25:30],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_25-30.pdf",sep=""),type="pdf")
   multiplot(plot_list[31:36],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_31-36.pdf",sep=""),type="pdf")
   multiplot(plot_list[37:42],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_37-42.pdf",sep=""),type="pdf")
   multiplot(plot_list[43:48],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_43-48.pdf",sep=""),type="pdf")
   multiplot(plot_list[49:54],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_49-54.pdf",sep=""),type="pdf")
   multiplot(plot_list[55:60],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_55-60.pdf",sep=""),type="pdf")
   multiplot(plot_list[61:66],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_61-66.pdf",sep=""),type="pdf")
   multiplot(plot_list[67:72],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_67-72.pdf",sep=""),type="pdf")   
   multiplot(plot_list[73:78],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_73-78.pdf",sep=""),type="pdf") 
   multiplot(plot_list[79:84],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_79-84.pdf",sep=""),type="pdf") 
   multiplot(plot_list[85:length(plot_list)],plotlist=NULL,cols=3)
   quartz.save(file=paste(fN,"_85-end.pdf",sep=""),type="pdf") 
 
   }
 
 
 
 
 
 
 
   plotBIG_by2 <- function(primes,ex1,ex2,ex3,ex4,ex5,ex6,ex7,fN){
           
 plot_list <- lapply(primes,plot_big_boy,ex1,ex2,ex3,ex4,ex5,ex6,ex7)

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
  
  
   multiplot(plot_list[1:4],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_1-4.pdf",sep=""),type="pdf")
   multiplot(plot_list[5:8],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_5-8.pdf",sep=""),type="pdf")
   multiplot(plot_list[9:12],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_9-12.pdf",sep=""),type="pdf")
   multiplot(plot_list[13:16],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_13-16.pdf",sep=""),type="pdf")
   multiplot(plot_list[17:20],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_17-20.pdf",sep=""),type="pdf")
   multiplot(plot_list[21:24],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_21-24.pdf",sep=""),type="pdf")
   multiplot(plot_list[25:28],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_25-28.pdf",sep=""),type="pdf")
   multiplot(plot_list[29:32],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_29-32.pdf",sep=""),type="pdf")
   multiplot(plot_list[33:36],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_33-36.pdf",sep=""),type="pdf")
   multiplot(plot_list[37:40],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_37-40.pdf",sep=""),type="pdf")
   multiplot(plot_list[41:44],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_41-44.pdf",sep=""),type="pdf")
   multiplot(plot_list[45:48],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_45-48.pdf",sep=""),type="pdf")   
   multiplot(plot_list[49:52],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_49-52.pdf",sep=""),type="pdf") 
   multiplot(plot_list[53:56],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_53-56.pdf",sep=""),type="pdf") 
   multiplot(plot_list[57:60],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_57-60.pdf",sep=""),type="pdf") 
   multiplot(plot_list[61:64],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_61-64.pdf",sep=""),type="pdf")  
   multiplot(plot_list[65:68],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_65-68.pdf",sep=""),type="pdf")    
   multiplot(plot_list[69:72],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_69-72.pdf",sep=""),type="pdf")   
   multiplot(plot_list[73:76],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_73-76.pdf",sep=""),type="pdf")    
   multiplot(plot_list[77:80],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_77-80.pdf",sep=""),type="pdf") 
   multiplot(plot_list[81:84],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_81-84.pdf",sep=""),type="pdf") 
   multiplot(plot_list[85:length(plot_list)],plotlist=NULL,cols=2)
   quartz.save(file=paste(fN,"_85-end.pdf",sep=""),type="pdf") 

 }