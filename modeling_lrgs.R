#### KCl simulator 
### create functions for all stimulants, in this case KCl
#### the fact that there may be different volumes of KCl administered to the cells for any given experiment (even same ones but bio reps) could work out to our advantage, in determining constants 

### creat a function for a 55mM concentration of KCl
## for a 55mM treatment

    ## we don't want any negative bounds, hence < 0 in both u and l; upp is whatever timepoint (in appropriate units eg min,sec) the KCl is taken off
    #### for LRGs that end with a positive slope upp can be made up
    ### this function does not explain whatever mRNA inhibition the LRGs may be experiencing after X amount of time


# have to read in like 5 sheets 
# trained data
# all_sus
# these are the genes
Pcsk1 <- all_sus[all_sus$primer %in% 'Pcsk1',]
Serpinb2 <- all_sus[all_sus$primer %in% 'Serpinb2',]
Bdnf <- all_sus[all_sus$primer %in% 'Bdnf',]
Fosl2 <- all_sus[all_sus$primer %in% 'Fosl2',]
Pax1 <- all_sus[all_sus$primer %in% 'Pax1',]
Acan <- all_sus[all_sus$primer %in% 'Acan',]
Pthlh <- all_sus[all_sus$primer %in% 'Pthlh',]
Crh <- all_sus[all_sus$primer %in% 'Crh',]

# 60 min
# all_60min
Pcsk1_60min <- all_60min[all_60min$primer %in% 'Pcsk1',]
Serpinb2_60min <- all_60min[all_60min$primer %in% 'Serpinb2',]
Bdnf_60min <- all_60min[all_60min$primer %in% 'Bdnf',]
Fosl2_60min <- all_60min[all_60min$primer %in% 'Fosl2',]
Pax1_60min <- all_60min[all_60min$primer %in% 'Pax1',]
Acan_60min <- all_60min[all_60min$primer %in% 'Acan',]
Pthlh_60min <- all_60min[all_60min$primer %in% 'Pthlh',]
Crh_60min <- all_60min[all_60min$primer %in% 'Crh',]


# testing data 

# 5 min
# all_5min 
# I think this is garbage
Pcsk1_5min <- all_5min[all_5min$primer %in% 'Pcsk1',]
Serpinb2_5min <- all_5min[all_5min$primer %in% 'Serpinb2',]
Bdnf_5min <- all_5min[all_5min$primer %in% 'Bdnf',]
Fosl2_5min <- all_5min[all_5min$primer %in% 'Fosl2',]
Pax1_5min <- all_5min[all_5min$primer %in% 'Pax1',]
Acan_5min <- all_5min[all_5min$primer %in% 'Acan',]
Pthlh_5min <- all_5min[all_5min$primer %in% 'Pthlh',]
Crh_5min <- all_5min[all_5min$primer %in% 'Crh',]

# 15 min
Pcsk1_15min <- all_15min[all_15min$primer %in% 'Pcsk1',]
Serpinb2_15min <- all_15min[all_15min$primer %in% 'Serpinb2',]
Bdnf_15min <- all_15min[all_15min$primer %in% 'Bdnf',]
Fosl2_15min <- all_15min[all_15min$primer %in% 'Fosl2',]
Pax1_15min <- all_15min[all_15min$primer %in% 'Pax1',]
Acan_15min <- all_15min[all_15min$primer %in% 'Acan',]
Pthlh_15min <- all_15min[all_15min$primer %in% 'Pthlh',]
Crh_15min <- all_15min[all_15min$primer %in% 'Crh',]


# 30 min
# all_30min
Pcsk1_30min <- all_30min[all_30min$primer %in% 'Pcsk1',]
Serpinb2_30min <- all_30min[all_30min$primer %in% 'Serpinb2',]
Bdnf_30min <- all_30min[all_30min$primer %in% 'Bdnf',]
Fosl2_30min <- all_30min[all_30min$primer %in% 'Fosl2',]
Pax1_30min <- all_30min[all_30min$primer %in% 'Pax1',]
Acan_30min <- all_30min[all_30min$primer %in% 'Acan',]
Pthlh_30min <- all_30min[all_30min$primer %in% 'Pthlh',]
Crh_30min <- all_30min[all_30min$primer %in% 'Crh',]



# 120 min
# all_120min 
Pcsk1_120min <- all_120min[all_120min$primer %in% 'Pcsk1',]
Serpinb2_120min <- all_120min[all_120min$primer %in% 'Serpinb2',]
Bdnf_120min <- all_120min[all_120min$primer %in% 'Bdnf',]
Fosl2_120min <- all_120min[all_120min$primer %in% 'Fosl2',]
Pax1_120min <- all_120min[all_120min$primer %in% 'Pax1',]
Acan_120min <- all_120min[all_120min$primer %in% 'Acan',]
Pthlh_120min <- all_120min[all_120min$primer %in% 'Pthlh',]
Crh_120min <- all_120min[all_120min$primer %in% 'Crh',]


  # if it desired to use smaller KCl time scales, set up to time when KCl was taken off
  
  integg2 <- function(t,a,b,con,upp){  
    
    u <- ifelse(t - a < 0 ,0, t - a)
    u <- ifelse(u > upp,upp, t - a) 
    l <- ifelse(t - (a+b) < 0,0, t - (a+b))  
    integrate(Vectorize(function(foo,x){x}),lower=l,upper=u,x=con)      
   
     }


# tailor this to the real max time point of whatever gene you have
# and if we're going to try fitting using a scaling constant which is calculated by dividng the max fold I of whatever gene, all you need to do is pull out that individual 
# timepoint and then calculate the
# the only thing that really matters at this point are the integration times 

times1 <- c(0,20,40,60,120,180,240,360)
times2 <- c(0,20,40,60,120,240,360)
times3 <- c(0,60,180,360)
times4 <- c(0,60,360)

# beta is calculated using max fold induction

Pcsk1_train <-  train(times2,55,360,60,Pcsk1,Pcsk1_60min)
Serpinb2_train <- train(times2,55,360,60,Serpinb2,Serpinb2_60min)
Bdnf_train <- train(times2,55,360,60,Bdnf,Bdnf_60min)
Fosl2_train <- train(times2,55,360,60,Fosl2,Fosl2_60min)
Pax1_train <- train(times2,55,360,60,Pax1,Pax1_60min)
Acan_train <- train(times2,55,360,60,Acan,Acan_60min)
Pthlh_train <- train(times2,55,360,60,Pthlh,Pthlh_60min)
Crh_train <- train(times2,55,360,60,Crh,Crh_60min)

Pcsk1_train_single <- single_train(times2,55,360,Pcsk1)
Serpinb2_train_single <- single_train(times2,55,360,Serpinb2)
Bdnf_train_single <- single_train(times2,55,360,Bdnf)
Fosl2_train_single <- single_train(times2,55,360,Fosl2)
Pax1_train_single <- single_train(times2,55,360,Pax1)
Acan_train_single <- single_train(times2,55,360,Acan)


int_times <- read.delim('LRG_lag_stim_time.txt')
colnames(int_times) <- c('gene','lag','stim','beta')
int_times$num <- 1:nrow(int_times)
int_times$num


   ggplot(int_times,aes(x=-stim,y=num,colour=gene)) +
      scale_y_continuous(breaks=int_times$num,labels=int_times$gene) +
      # scale_x_reverse() +
      xlim(-400,0) +
      geom_segment(aes(xend=-lag,ystart=num,yend=num),size=5) +
      xlab('Integration Window (min)') +
      ylab('Genes') +
      opts(title='LRGs')

#train
Pcsk1_60_train <- test(46.53,173.53,times2,60,55,'test',0.002389883,Pcsk1_60min,'real','60 min train (Pcsk1)','60_min_pcsk1.pdf')
Pcsk1_360_train <- test(46.53,173.53,times2,360,55,'test',0.002389883,Pcsk1,'real','360 min train (Pcsk1)','360_min_psck1.pdf')
# everything that has been trained from 60 and 360 min
# Pcsk1 test
Psck1_5_test <- test(46.53,173.53,times2,5,55,'test',0.002389883,Pcsk1_5min,'real','5 min test (Pcsk1)','5_min_pcsk1.pdf')
Psck1_15_test <- test(46.53,173.53,times2,15,55,'test',0.002389883,Pcsk1_15min,'real','15 min test (Pcsk1)','15_min_pcsk1.pdf')
Psck1_30_test <- test(46.53,173.53,times2,30,55,'test',0.002389883,Pcsk1_30min,'real','30 min test (Pcsk1)','30_min_pcsk1.pdf')
Psck1_120_test <- test(46.53,173.53,times2,120,55,'test',0.002389883,Pcsk1_120min,'real','120 min test (Pcsk1)','120_min_pcsk1.pdf')


Pcsk1_360_train <- test(53,140.23,times2,360,55,'test',0.004482713,Pcsk1,'real','360 min train (Pcsk1)','360_min_psck1.pdf')
Pcsk1_60_train <- test(53,140.23,times2,60,55,'test',0.004482713,Pcsk1,'real','60 min train (Pcsk1)','60_min_psck1.pdf')

Pcsk1_5_train <- test(53,140.23,times2,5,55,'test',0.004482713,Pcsk1,'real','5 min train (Pcsk1)','5_min_psck1.pdf')
Pcsk1_15_train <- test(53,140.23,times2,15,55,'test',0.004482713,Pcsk1,'real','15 min train (Pcsk1)','15_min_psck1.pdf')
Pcsk1_30_train <- test(53,140.23,times2,30,55,'test',0.004482713,Pcsk1,'real','30 min train (Pcsk1)','30_min_psck1.pdf')
Pcsk1_120_train <- test(53,140.23,times2,120,55,'test',0.004482713,Pcsk1,'real','120 min train (Pcsk1)','120_min_psck1.pdf')

Pcsk1_60_train_single <- single_train(46.53)



      test <- function(a,b,t,upp,con,name,s,real,name2,titulo,bigt){

dt.a <- data.table(a,k=1,key='k')
dt.b <- data.table(b,k=1,key='k')
dt.t <- data.table(t,k=1,key='k')
cmb <- dt.a[dt.b[dt.t]]   
cmb[,u := ifelse(t - a < 0,0,t - a)] 
cmb[,u := ifelse(u > upp,upp,u)]
cmb[,l := ifelse(t - (a+b) < 0,0, t - (a+b))]
cmb[,l := ifelse(l > upp,upp,l)]
cmb[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))]
cmb <- as.data.frame(cmb)
cmb$value <- cmb$value * s
cmb <- cmb[,c('t','value')]
cmb$sde <- NA
cmb$KCl_duration <- name
cmb$value <- cmb$value + 1
real <- real[,-1]
colnames(real) = c('t','value','sde')
real$KCl_duration <- name2
yo <- rbind(cmb,real) 	    
    
    pdf(bigt,onefile=T)
 	  g <-  ggplot(yo,aes(x=t,y= value,ymax = value + sde, ymin = value  - sde,colour=KCl_duration)) + 
           geom_errorbar(width=2) + 
           geom_point() +
           geom_line() +
           xlab('Time (min)') +
           ylab('Fold Induction') +          
           opts(title = titulo)    
   print(g)
dev.off()   
}

      test_single <- function(a,b,t,upp,con,name,s,real,name2,titulo,bigt){

dt.a <- data.table(a,k=1,key='k')
dt.b <- data.table(b,k=1,key='k')
dt.t <- data.table(t,k=1,key='k')
cmb <- dt.a[dt.b[dt.t]]   
cmb[,u := ifelse(t - a < 0,0,t - a)] 
cmb[,u := ifelse(u > upp,upp,u)]
cmb[,l := ifelse(t - (a+b) < 0,0, t - (a+b))]
cmb[,l := ifelse(l > upp,upp,l)]
cmb[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))]
cmb <- as.data.frame(cmb)
cmb$value <- cmb$value * s
cmb <- cmb[,c('t','value')]
cmb$sde <- NA
cmb$KCl_duration <- name
cmb$value <- cmb$value + 1
real <- real[,-1]
colnames(real) = c('t','value','sde')
real$KCl_duration <- name2
yo <- rbind(cmb,real) 	    
    
    pdf(bigt,onefile=T)
 	  g <-  ggplot(yo,aes(x=t,y= value,ymax = value + sde, ymin = value  - sde,colour=KCl_duration)) + 
           geom_errorbar(width=2) + 
           geom_point() +
           geom_line() +
           xlab('Time (min)') +
           ylab('Fold Induction') +          
           opts(title = titulo)    
   print(g)
dev.off()   
}


              train <- function(t,con,upp,upp2,gene,gene2){ # it's about a 45 sec run time

##### 1st run
a <- seq(from=0,to=360,by=10) # 1st resolution will be 360 min (6 hrs) this is more than enough time
b <- seq(from=0,to=1440,by=10) # 1st resolution will be 1440 min (24 hrs) 
   

dt.a <- data.table(a,k=1,key='k')
dt.b <- data.table(b,k=1,key='k')
dt.t <- data.table(t,k=1,key='k')
cmb <- dt.a[dt.b[dt.t]]   
cmba <- dt.a[dt.b[dt.t]]
cmb[,u := ifelse(t - a < 0,0,t - a)]
cmba[,u := ifelse(t - a < 0,0,t - a)] 
cmb[,u := ifelse(u > upp,upp,u)]
cmba[,u := ifelse(u > upp2,upp2,u)]
cmb[,l := ifelse(t - (a+b) < 0,0, t - (a+b))]
cmba[,l := ifelse(t - (a+b) < 0,0, t - (a+b))]
cmb[,l := ifelse(l > upp,upp,l)]
cmba[,l := ifelse(l > upp2,upp2,l)]
cmb[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))] # 7.9 seconds
cmba[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))]
cmb <- as.data.frame(cmb)
cmba <- as.data.frame(cmba)
cmb <- cmb[order(cmb$a,cmb$b),] # I order here so I can make ordered matrices for the least squares fitting
cmba <- cmba[order(cmba$a,cmba$b),]
cmb <- cmb[,-1]
cmba <- cmba[,-1]

# everything else will have to be done twice on the training

# in the case of Pcsk1 360 time point is the maximum time point
# make a condensed "subset" of the population as far as newly fitted Ss go
# as mentioned above w/ the scaling constant, sub will only be set once, so use 'sub' again w/ increased resolution
genz <- rbind(gene,gene2)
sub <- cmb[cmb$t == genz[genz$mean == max(genz$mean),'timepoints'],] # you want the max fold induction time point for the real data, then solve for s using that time point
# can condense this part because only the integration times will affect the scaling constant values (s)
sub <- sub[!duplicated(sub$b),]
# find the s values

sub$s <- max(genz$mean) / sub$value
# match the s values to original 
   cmb$s <-  sub$s[match(cmb$b,sub$b)] 
   cmba$s <- sub$s[match(cmba$b,sub$b)]
# find the new values that are multiplied by the scaling constant
cmb$value <- cmb$value * cmb$s
cmba$value <- cmba$value * cmba$s
# here's a least squares function that will print lag times, int time periods, and scaling constants
# this is working to 
train <- rbind(cmb,cmba)
train <- train[order(train$a,train$b),]
first <- train[seq(from=1, to=nrow(train), by=length(t)*2),][which.min(colSums((genz$mean - matrix(train$value, nrow=length(t)*2))**2)),] # find minimized least squares


######### 2nd run 
# go 50 min in either direction for both lag and integration time window increasing the resolution now by 1 min

a2 <- seq(from=first$a - 50,to=first$a + 50,by=1)
a2 <- a2[a2 >= 0] # only positivity
b2 <- seq(from=first$b - 50,to=first$b + 50,by=1)
b2 <- b2[b2 >= 0] # only positivity
  

dt.a2 <- data.table(a2,k=1,key='k')
dt.b2 <- data.table(b2,k=1,key='k')
cmb2 <- dt.a2[dt.b2[dt.t]] 
cmb2a <- dt.a2[dt.b2[dt.t]]  
cmb2[,u := ifelse(t - a2 < 0,0,t - a2)]
cmb2a[,u := ifelse(t - a2 < 0,0,t - a2)]
cmb2[,u := ifelse(u > upp,upp,u)]
cmb2a[,u := ifelse(u > upp2,upp2,u)]
cmb2[,l := ifelse(t - (a2+b2) < 0,0, t - (a2+b2))]
cmb2a[,l := ifelse(t - (a2+b2) < 0,0, t - (a2+b2))]
cmb2[,l := ifelse(l > upp,upp,l)]
cmb2a[,l := ifelse(l > upp2,upp2,l)]
cmb2[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))] # 14.9 seconds
cmb2a[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))]
cmb2 <- as.data.frame(cmb2)
cmb2a <- as.data.frame(cmb2a)
cmb2 <- cmb2[order(cmb2$a2,cmb2$b2),]
cmb2a <- cmb2a[order(cmb2a$a2,cmb2a$b2),]
cmb2 <- cmb2[,-1]
cmb2a <- cmb2a[,-1]
sub2 <- cmb2[cmb2$t == genz[genz$mean == max(genz$mean),'timepoints'],]
sub2 <- sub2[!duplicated(sub2$b2),]
sub2$s <- max(genz$mean) / sub2$value
cmb2$s <-  sub2$s[match(cmb2$b2,sub2$b2)]
cmb2a$s <- sub2$s[match(cmb2a$b2,sub2$b2)] 
cmb2$value <- cmb2$value * cmb2$s
cmb2a$value <- cmb2a$value * cmb2a$s
train2 <- rbind(cmb2,cmb2a)
train2 <- train2[order(train2$a,train2$b),]
second <- train2[seq(from=1, to=nrow(train2), by=length(t)*2),1:2][which.min(colSums((genz$mean - matrix(train2$value, nrow=length(t)*2))**2)),] # real - data ^2
# for best fit go out 1 min in either direction now increasing resolution to seconds
# I'm not convinced that this is the best way to make 'second' arrays, feel free to edit 

    if (second$a2 < 1){ # otherwise you'll get some infinite bounds
a3 <- seq(from=strptime(paste((0)%/%60, (0)%%60, sep=':'), format='%H:%M'),to= strptime(paste((second$a2+1)%/%60, (second$a2+1)%%60, sep=':'), format='%H:%M'),by='sec')
a3_min <- as.numeric(format(a3,'%M.%S'))
a3_hr <- as.numeric(format(a3,'%H')) * 60
a3 <- a3_min + a3_hr
a3 <- a3[a3 >= 0]
 }
    else {
	a3 <- seq(from=strptime(paste((second$a2 - 1)%/%60, (second$a2 - 1)%%60, sep=':'), format='%H:%M'),to=strptime(paste((second$a2+1)%/%60, (second$a2+1)%%60, sep=':'), format='%H:%M'),by='sec')
a3_min <- as.numeric(format(a3,'%M.%S'))
a3_hr <- as.numeric(format(a3,'%H')) * 60
a3 <- a3_min + a3_hr
a3 <- a3[a3 >= 0]
 }
    if (second$b2 < 1){
b3 <- seq(from=strptime(paste((0)%/%60, (0)%%60, sep=':'), format='%H:%M'),to= strptime(paste((second$b2+1)%/%60, (second$b2+1)%%60, sep=':'), format='%H:%M'),by='sec')
b3_min <- as.numeric(format(b3,'%M.%S')) 
b3_hr <- as.numeric(format(b3,'%H')) * 60
b3 <- b3_min + b3_hr
b3 <- b3[b3 >= 0]
 }
    else {
b3 <- seq(from=strptime(paste((second$b2 - 1)%/%60, (second$b2 - 1)%%60, sep=':'), format='%H:%M'),to= strptime(paste((second$b2+1)%/%60, (second$b2+1)%%60, sep=':'), format='%H:%M'),by='sec')
b3_min <- as.numeric(format(b3,'%M.%S')) 
b3_hr <- as.numeric(format(b3,'%H')) * 60
b3 <- b3_min + b3_hr
b3 <- b3[b3 >= 0]	
 }
dt.a3 <- data.table(a3,k=1,key='k')
dt.b3 <- data.table(b3,k=1,key='k')
cmb3 <- dt.a3[dt.b3[dt.t]]   
cmb3a <- dt.a3[dt.b3[dt.t]]
cmb3[,u := ifelse(t - a3 < 0,0,t - a3)]
cmb3a[,u := ifelse(t - a3 < 0,0,t - a3)]
cmb3[,u := ifelse(u > upp,upp,u)]
cmb3a[,u := ifelse(u > upp2,upp2,u)]
cmb3[,l := ifelse(t - (a3+b3) < 0,0, t - (a3+b3))]
cmb3a[,l := ifelse(t - (a3+b3) < 0,0,t - (a3+b3))]
cmb3[,l := ifelse(l > upp,upp,l)]
cmb3a[,l := ifelse(l > upp2,upp2,l)]
cmb3[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))]
cmb3a[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))]
cmb3 <- as.data.frame(cmb3)
cmb3a <- as.data.frame(cmb3a)
cmb3 <- cmb3[order(cmb3$a3,cmb3$b3),]
cmb3a <- cmb3a[order(cmb3a$a3,cmb3$b3),]
cmb3 <- cmb3[,-1]
cmb3a <- cmb3a[,-1]
sub3 <- cmb3[cmb3$t == genz[genz$mean == max(genz$mean),'timepoints'],]
sub3 <- sub3[!duplicated(sub3$b3),]
sub3$s <- max(genz$mean) / sub3$value
cmb3$s <-  sub3$s[match(cmb3$b3,sub3$b3)] 
cmb3a$s <- sub3$s[match(cmb3a$b3,sub3$b3)]
cmb3$value <- cmb3$value * cmb3$s
cmb3a$value <- cmb3a$value * cmb3a$s
train3 <- rbind(cmb3,cmb3a)
train3 <- train3[order(train3$a,train3$b),]
colnames(train3) = c('lag_time','stim_intval','timepoints','u','l','value','s')
final <- train3[seq(from=1, to=nrow(train3), by=length(t)*2),1:2][which.min(colSums((genz$mean - matrix(train3$value, nrow=length(t)*2))**2)),]
sim <- train3[train3$lag_time == final$lag_time & train3$stim_intval == final$stim_intval,]
sim
}






  single_train <- function(t,con,upp,gene){ # it's about a 45 sec run time

##### 1st run
a <- seq(from=0,to=360,by=10) # 1st resolution will be 360 min (6 hrs) this is more than enough time
b <- seq(from=0,to=1440,by=10) # 1st resolution will be 1440 min (24 hrs) 
   

dt.a <- data.table(a,k=1,key='k')
dt.b <- data.table(b,k=1,key='k')
dt.t <- data.table(t,k=1,key='k')
cmb <- dt.a[dt.b[dt.t]]   
cmb[,u := ifelse(t - a < 0,0,t - a)]
cmb[,u := ifelse(u > upp,upp,u)]
cmb[,l := ifelse(t - (a+b) < 0,0, t - (a+b))]
cmb[,l := ifelse(l > upp,upp,l)]
cmb[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))] # 7.9 seconds
cmb <- as.data.frame(cmb)
cmb <- cmb[order(cmb$a,cmb$b),] # I order here so I can make ordered matrices for the least squares fitting
cmb <- cmb[,-1]

# everything else will have to be done twice on the training

# in the case of Pcsk1 360 time point is the maximum time point
# make a condensed "subset" of the population as far as newly fitted Ss go
# as mentioned above w/ the scaling constant, sub will only be set once, so use 'sub' again w/ increased resolution
sub <- cmb[cmb$t == gene[gene$mean == max(gene$mean),'timepoints'],] # you want the max fold induction time point for the real data, then solve for s using that time point
# can condense this part because only the integration times will affect the scaling constant values (s)
sub <- sub[!duplicated(sub$b),]
# find the s values

sub$s <- max(gene$mean) / sub$value
# match the s values to original 
   cmb$s <-  sub$s[match(cmb$b,sub$b)] 
# find the new values that are multiplied by the scaling constant
cmb$value <- cmb$value * cmb$s
# here's a least squares function that will print lag times, int time periods, and scaling constants
# this is working to 
first <- cmb[seq(from=1, to=nrow(cmb), by=length(t)),][which.min(colSums((gene$mean - matrix(cmb$value, nrow=length(t)))**2)),] # find minimized least squares


######### 2nd run 
# if 1st run best fit is within 10 min of parameter then go 80 min in that direction, otherwise
# go 50 min in either direction for both lag and integration time window increasing the resolution now by 1 min

a2 <- seq(from=first$a - 50,to=first$a + 50,by=1)
a2 <- a2[a2 >= 0] # only positivity
b2 <- seq(from=first$b - 50,to=first$b + 50,by=1)
b2 <- b2[b2 >= 0] # only positivity
  

dt.a2 <- data.table(a2,k=1,key='k')
dt.b2 <- data.table(b2,k=1,key='k')
cmb2 <- dt.a2[dt.b2[dt.t]] 
cmb2[,u := ifelse(t - a2 < 0,0,t - a2)]
cmb2[,u := ifelse(u > upp,upp,u)]
cmb2[,l := ifelse(t - (a2+b2) < 0,0, t - (a2+b2))]
cmb2[,l := ifelse(l > upp,upp,l)]
cmb2[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))] # 14.9 seconds
cmb2 <- as.data.frame(cmb2)
cmb2 <- cmb2[order(cmb2$a2,cmb2$b2),]
cmb2 <- cmb2[,-1]
sub2 <- cmb2[cmb2$t == gene[gene$mean == max(gene$mean),'timepoints'],]
sub2 <- sub2[!duplicated(sub2$b2),]
sub2$s <- max(gene$mean) / sub2$value
cmb2$s <-  sub2$s[match(cmb2$b2,sub2$b2)]
cmb2$value <- cmb2$value * cmb2$s
second <- cmb2[seq(from=1, to=nrow(cmb2), by=length(t)),1:2][which.min(colSums((gene$mean - matrix(cmb2$value, nrow=length(t)))**2)),] # real - data ^2
# for best fit go out 1 min in either direction now increasing resolution to seconds
# I'm not convinced that this is the best way to make 'second' arrays, feel free to edit 

    if (second$a2 < 1){ # otherwise you'll get some infinite bounds
a3 <- seq(from=strptime(paste((0)%/%60, (0)%%60, sep=':'), format='%H:%M'),to= strptime(paste((second$a2+1)%/%60, (second$a2+1)%%60, sep=':'), format='%H:%M'),by='sec')
a3_min <- as.numeric(format(a3,'%M.%S'))
a3_hr <- as.numeric(format(a3,'%H')) * 60
a3 <- a3_min + a3_hr
a3 <- a3[a3 >= 0]
 }
    else {
	a3 <- seq(from=strptime(paste((second$a2 - 1)%/%60, (second$a2 - 1)%%60, sep=':'), format='%H:%M'),to=strptime(paste((second$a2+1)%/%60, (second$a2+1)%%60, sep=':'), format='%H:%M'),by='sec')
a3_min <- as.numeric(format(a3,'%M.%S'))
a3_hr <- as.numeric(format(a3,'%H')) * 60
a3 <- a3_min + a3_hr
a3 <- a3[a3 >= 0]
 }
    if (second$b2 < 1){
b3 <- seq(from=strptime(paste((0)%/%60, (0)%%60, sep=':'), format='%H:%M'),to= strptime(paste((second$b2+1)%/%60, (second$b2+1)%%60, sep=':'), format='%H:%M'),by='sec')
b3_min <- as.numeric(format(b3,'%M.%S')) 
b3_hr <- as.numeric(format(b3,'%H')) * 60
b3 <- b3_min + b3_hr
b3 <- b3[b3 >= 0]
 }
    else {
b3 <- seq(from=strptime(paste((second$b2 - 1)%/%60, (second$b2 - 1)%%60, sep=':'), format='%H:%M'),to= strptime(paste((second$b2+1)%/%60, (second$b2+1)%%60, sep=':'), format='%H:%M'),by='sec')
b3_min <- as.numeric(format(b3,'%M.%S')) 
b3_hr <- as.numeric(format(b3,'%H')) * 60
b3 <- b3_min + b3_hr
b3 <- b3[b3 >= 0]	
 }
dt.a3 <- data.table(a3,k=1,key='k')
dt.b3 <- data.table(b3,k=1,key='k')
cmb3 <- dt.a3[dt.b3[dt.t]]   
cmb3[,u := ifelse(t - a3 < 0,0,t - a3)]
cmb3[,u := ifelse(u > upp,upp,u)]
cmb3[,l := ifelse(t - (a3+b3) < 0,0, t - (a3+b3))]
cmb3[,l := ifelse(l > upp,upp,l)]
cmb3[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))]
cmb3 <- as.data.frame(cmb3)
cmb3 <- cmb3[order(cmb3$a3,cmb3$b3),]
cmb3 <- cmb3[,-1]
sub3 <- cmb3[cmb3$t == gene[gene$mean == max(gene$mean),'timepoints'],]
sub3 <- sub3[!duplicated(sub3$b3),]
sub3$s <- max(gene$mean) / sub3$value
cmb3$s <-  sub3$s[match(cmb3$b3,sub3$b3)] 
cmb3$value <- cmb3$value * cmb3$s
colnames(cmb3) = c('lag_time','stim_intval','timepoints','u','l','value','s')
final <- cmb3[seq(from=1, to=nrow(cmb3), by=length(t)),1:2][which.min(colSums((gene$mean - matrix(cmb3$value, nrow=length(t)))**2)),]
sim <- cmb3[cmb3$lag_time == final$lag_time & cmb3$stim_intval == final$stim_intval,]
sim
}






 




 