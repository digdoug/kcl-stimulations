 
 
 times <- c(0,20,40,60,120,240,360)
 
train_360_60 <-  train(times,55,360,60,Pcsk1,Pcsk1_60min)
train_360_120 <- train(times,55,360,120,Pcsk1,Pcsk1_120min)

train_360_30 <-  train(times,55,360,30,Pcsk1,Pcsk1_30min) # my assumption is that the 'to' argument in the sequence is negative

train_360_15 <- train(times,55,360,15,Pcsk1,Pcsk1_15min)


           ggplot(int_times,aes(x=stim,y= num,colour=gene)) +
           scale_y_continuous(labels=c())
           geom_hline() + 
           xlab('x') +
           ylab('y') +          
           opts(title = 'Multiple Gs')    
           
           ggplot(int_times,aes(x=stim,y=num,colour=gene)) +
           scale_y_continuous(labels=int_times$gene) +
           geom_segment(aes(xend=length(stim),ystart=num,yend=num)) +
           xlab('x') +
           ylab('y') +
           opts(title='multi gs')
           
   
           
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
genes <- rbind(gene,gene2)
sub <- cmb[cmb$t == genes[genes$mean == max(genes$mean),'timepoints'],] # you want the max fold induction time point for the real data, then solve for s using that time point
suba <- cmba[cmba$t == genes[genes$mean == max(genes$mean),'timepoints'],]
# can condense this part because only the integration times will affect the scaling constant values (s)
sub <- sub[!duplicated(sub$b),]
suba <- suba[!duplicated(suba$b),] # need to make sure that ALL beta values only depend on b for this
# find the s values

sub$s <- max(genes$mean) / sub$value
suba$s <- max(genes$mean) / suba$value
# match the s values to original 
   cmb$s <-  sub$s[match(cmb$b,sub$b)] 
   cmba$s <- suba$s[match(cmba$b,suba$b)]

# find the new values that are multiplied by the scaling constant
cmb$value <- cmb$value * cmb$s
cmba$value <- cmba$value * cmba$s
# here's a least squares function that will print lag times, int time periods, and scaling constants
# this is working to 
train <- rbind(cmb,cmba)
train <- train[order(train$a,train$b),]
first <- train[seq(from=1, to=nrow(train), by=length(t)*2),1:2][which.min(colSums((genes$mean - matrix(train$value, nrow=length(t)*2))**2)),] # find minimized least squares

######### 2nd run 
# if 1st run best fit is within 10 min of parameter then go 80 min in that direction, otherwise
# go 50 min in either direction for both lag and integration time window increasing the resolution now by 1 min

a2 <- seq(from=first$a - 50,to=first$a + 50,by=1)
b2 <- seq(from=first$b - 50,to=first$b + 50,by=1)

  

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
cmb2[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))] # 14.9 seconds
cmb2a[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))]
cmb2 <- as.data.frame(cmb2)
cmb2a <- as.data.frame(cmb2a)
cmb2 <- cmb2[order(cmb2$a2,cmb2$b2),]
cmb2a <- cmb2a[order(cmb2a$a2,cmb2a$b2),]
cmb2 <- cmb2[,-1]
cmb2a <- cmb2a[,-1]

# going to have to do this again because we now have different time windows
sub2 <- cmb2[cmb2$t == genes[genes$mean == max(genes$mean),'timepoints'],]
sub2a <- cmb2a[cmb2a$t == genes[genes$mean == max(genes$mean),'timepoints'],]
# can condense this part because only the integration times will affect the scaling constant values (s)
sub2 <- sub2[!duplicated(sub2$b2),]
sub2a <- sub2a[!duplicated(sub2a$b2),]
# find the s values
sub2$s <- max(genes$mean) / sub2$value
sub2a$s <- max(genes$mean) / sub2a$value
# match the s values to 
   cmb2$s <-  sub2$s[match(cmb2$b2,sub2$b2)]
   cmb2a$s <- sub2a$s[match(cmb2a$b2,sub2a$b2)] 
# find the new values that are multiplied by the scaling constant
cmb2$value <- cmb2$value * cmb2$s
cmb2a$value <- cmb2a$value * cmb2a$s
# here's a least squares function that will print lag times and int time periods
train2 <- rbind(cmb2,cmb2a)
train2 <- train2[order(train2$a,train2$b),]
second <- train2[seq(from=1, to=nrow(train2), by=length(t)*2),1:2][which.min(colSums((genes$mean - matrix(train2$value, nrow=length(t)*2))**2)),] # real - data ^2

# i'm repeating all the code so that someone can follow me on this particular example w/ parameter selection

# for best fit go out 3 min in either direction now increasing resolution to seconds

# now we're going 3 minutes out on either side by seconds 
# I'm not convinced that this is the best way to make 'second' arrays, feel free to edit 

a3 <- seq(strptime(paste((second$a2 - 1)%/%60, (second$a2 - 1)%%60, sep=':'), format='%H:%M'), strptime(paste((second$a2+1)%/%60, (second$a2+1)%%60, sep=':'), format='%H:%M'),by='sec')
a3_min <- as.numeric(format(a3,'%M.%S'))
a3_hr <- as.numeric(format(a3,'%H')) * 60
a3 <- a3_min + a3_hr

b3 <- seq(strptime(paste((second$b2 - 1)%/%60, (second$b2 - 1)%%60, sep=':'), format='%H:%M'), strptime(paste((second$b2+1)%/%60, (second$b2+1)%%60, sep=':'), format='%H:%M'),by='sec')
b3_min <- as.numeric(format(b3,'%M.%S')) 
b3_hr <- as.numeric(format(b3,'%H')) * 60
b3 <- b3_min + b3_hr

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
cmb3[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))]
cmb3a[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))]
cmb3 <- as.data.frame(cmb3)
cmb3a <- as.data.frame(cmb3a)
cmb3 <- cmb3[order(cmb3$a3,cmb3$b3),]
cmb3a <- cmb3a[order(cmb3a$a3,cmb3$b3),]
cmb3 <- cmb3[,-1]
cmb3a <- cmb3a[,-1]

# going to have to do this again because we now have different time windows
sub3 <- cmb3[cmb3$t == genes[genes$mean == max(genes$mean),'timepoints'],]
sub3a <- cmb3a[cmb3a$t == genes[genes$mean == max(genes$mean),'timepoints'],]
# can condense this part because only the integration times will affect the scaling constant values (s)
sub3 <- sub3[!duplicated(sub3$b3),]
sub3a <- sub3a[!duplicated(sub3a$b3),]
# find the s values
sub3$s <- max(genes$mean) / sub3$value
sub3a$s <- max(genes$mean) / sub3a$value
# match the s values to 
   cmb3$s <-  sub3$s[match(cmb3$b3,sub3$b3)] 
   cmb3a$s <- sub3a$s[match(cmb3a$b3,sub3a$b3)]
# find the new values that are multiplied by the scaling constant
cmb3$value <- cmb3$value * cmb3$s
cmb3a$value <- cmb3a$value * cmb3a$s
# here's a least squares function that will print lag times and int time periods
# this is working to 
train3 <- rbind(cmb3,cmb3a)
train3 <- train3[order(train3$a,train3$b),]
colnames(train3) = c('lag_time','stim_intval','timepoints','u','l','value','s')

final <- train3[seq(from=1, to=nrow(train3), by=length(t)*2),1:2][which.min(colSums((genes$mean - matrix(train3$value, nrow=length(t)*2))**2)),]
sim <- train3[train3$lag_time == final$lag_time & train3$stim_intval == final$stim_intval,]
sim
}




