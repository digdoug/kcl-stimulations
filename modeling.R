fos <- all_chx[all_chx$primer == 'Fos',]
fos

structure(list(primer = c("Fos", "Fos", "Fos", "Fos", "Fos", 
"Fos", "Fos"), timepoints = c(0, 20, 40, 60, 120, 240, 360), 
    mean = c(1, 57.1037082843044, 138.4460661954, 129.729013143179, 
    175.237254456402, 226.101301528076, 230.734845935663), sde = c(0, 
    12.3355423113076, 28.2064805069725, 36.3619619198529, 60.2678462977631, 
    24.3781604232449, 82.9428307985539), type = c("KCl + CHX", 
    "KCl + CHX", "KCl + CHX", "KCl + CHX", "KCl + CHX", "KCl + CHX", 
    "KCl + CHX")), .Names = c("primer", "timepoints", "mean", 
"sde", "type"), row.names = 127:133, class = "data.frame")

# i realize this is not the equation I will be using in reality. This is just basic infra structure of least squares fitting

    train <- function(t,con,upp,upp2,gene,gene2){ # it's about a 45 sec run time


##### 1st run

# make possible arrays of whatever parameter, in this case we know it's minutes and these intervals should cover the reality
a <- seq(from=0,to=360,by=10) # 1st  resolution will be 360 min (6 hrs) this is more than enough time
b <- seq(from=0,to=1440,by=10) # 1st resolution will be 1440 min (24 hrs) 
   
# make data tables for each array of parameters
dt.a <- data.table(a,k=1,key='k')
dt.b <- data.table(b,k=1,key='k')
dt.t <- data.table(t,k=1,key='k')

# combine all parameter data tables. This will give you one final data frame with all possibilites of k parameters
cmb <- dt.a[dt.b[dt.t]]   
cmba <- dt.a[dt.b[dt.t]] # second data frame is from the second experiment
cmb[,u := ifelse(t > upp,upp,t)] # u will be the upper bound in the integral, and obviously we don't want negative bounds
cmba[,u := ifelse(t > upp2,upp2,t)] 
cmb[,l := 0]
cmba[,l := 0]
cmb[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,amax,ka,a,n,beta,data){amax * (ka)}),a=con))] # in R you have to Vectorize constants if you want to integrate them
cmba[,value := mapply(function(...){integrate(...)$value},lower=l,upper=u,MoreArgs=list(f=Vectorize(function(x,constant){constant}),constant=con))]
cmb <- as.data.frame(cmb)
cmba <- as.data.frame(cmba)
cmb <- cmb[order(cmb$a,cmb$b),] # I order here so I can make ordered matrices for the least squares fitting
cmba <- cmba[order(cmba$a,cmba$b),]
cmb <- cmb[,-1]
cmba <- cmba[,-1]
################# Previously I would solve for s (the scaling constant) here, but now I wont be doing that now. 

train <- rbind(cmb,cmba)
train <- train[order(train$a,train$b),]
# here's a least squares line that will print k parameters
first <- train[seq(from=1, to=nrow(train), by=length(t)*2),][which.min(colSums((genz$mean - matrix(train$value, nrow=length(t)*2))**2)),] # find minimized least squares


######### 2nd run 
# go 50 min in either direction for both lag and integration time window increasing the resolution now by 1 min

# intervals will vary depending on the parameter
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
train3 <- rbind(cmb3,cmb3a)
train3 <- train3[order(train3$a,train3$b),]
colnames(train3) = c('lag_time','stim_intval','timepoints','u','l','value','s')
final <- train3[seq(from=1, to=nrow(train3), by=length(t)*2),1:2][which.min(colSums((genz$mean - matrix(train3$value, nrow=length(t)*2))**2)),]
sim <- train3[train3$lag_time == final$lag_time & train3$stim_intval == final$stim_intval,]
sim
}
