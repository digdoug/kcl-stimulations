  # for now the lower bound will be 0 but the lower bound can easily be more generalized
  
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
Pcsk1_60min <- all_60min[all_60min$primer %in% 'Pcsk1',]
Serpinb2_60min <- all_60min[all_60min$primer %in% 'Serpinb2',]
Bdnf_60min <- all_60min[all_60min$primer %in% 'Bdnf',]
Fosl2_60min <- all_60min[all_60min$primer %in% 'Fosl2',]
Pax1_60min <- all_60min[all_60min$primer %in% 'Pax1',]
Acan_60min <- all_60min[all_60min$primer %in% 'Acan',]
Pthlh_60min <- all_60min[all_60min$primer %in% 'Pthlh',]
Crh_60min <- all_60min[all_60min$primer %in% 'Crh',]

Pcsk1[Pcsk1$timepoints == t, 'mean']  
  
  integrand <- Vectorize(function(x){scale * ((al * (t - tau)) - b * dat)})
  al <- 55
  t <- 60
  tau <- 50 
  b <- 0.01
  dat <- Pcsk1[Pcsk1$timepoints == t, 'mean']
  scale <-  0.0001
  
  integrate(integrand,0,t)
  times <- c(0,20,40,60,120,240,360)
 
 # Pcsk1 360 min
 structure(list(primer = c("Pcsk1", "Pcsk1", "Pcsk1", "Pcsk1", 
"Pcsk1", "Pcsk1", "Pcsk1"), timepoints = c(0, 20, 40, 60, 120, 
240, 360), mean = c(1, 0.720264228850627, 0.800030483103659, 
1.7430986348151, 16.5172242381516, 33.9140055666552, 34.57359836929
), sde = c(0, 0.0301589738399932, 0.167577682741619, 0.534133256883202, 
3.63347764961647, 3.14731189930687, 3.48532358368367)), .Names = c("primer", 
"timepoints", "mean", "sde"), row.names = c(413L, 414L, 415L, 
416L, 417L, 419L, 420L), class = "data.frame")


# Pcsk1 60 min
structure(list(primer = c("Pcsk1", "Pcsk1", "Pcsk1", "Pcsk1", 
"Pcsk1", "Pcsk1", "Pcsk1"), timepoints = c(0, 20, 40, 60, 120, 
240, 360), mean = c(1, 0.754077948594051, 0.764135807165231, 
0.869415006051656, 16.5431465447981, 16.749349274373, 6.21452582124456
), sde = c(0, 0.0254695904347977, 0.076465889809583, 0.198010119220461, 
14.0651396033096, 0.718440892568028, 0.973327784130113)), .Names = c("primer", 
"timepoints", "mean", "sde"), row.names = 450:456, class = "data.frame")
 
 timez <- c(0,20,40,60,120,240,360)
 
 
 test <- save[order(save$scalz,save$tauz,save$betaz,save$V1),]
 head(test,n=21L)
  
save <- train(timez,55,Pcsk1,Pcsk1_60min,360,60)

  Pcsk1_60min

  train <- function(tz,con,data,data2,up,up2){ # it's about a 45 sec run time


##### 1st run
scalz <- seq(from=0.00001,to=0.01,by=0.001)
tauz <- seq(from=0,to=360,by=10) # 1st resolution will be 360 min (6 hrs) this is more than enough time
betaz <- seq(from=0.00001,to=0.01,by=0.001) # 1st resolution will be 1440 min (24 hrs) 
   
dt.scalz <- data.table(scalz,k=1,key='k')
dt.tauz <- data.table(tauz,k=1,key='k')
dt.betaz <- data.table(betaz,k=1,key='k')
dt.tz <- data.table(tz,k=1,key='k')
dt.dataz <- data.table(data$mean,k=1,key='k')
dt.dataz2 <- data.table(data2$mean,k=1,key='k')
cmb <- dt.scalz[dt.tauz[dt.betaz[dt.dataz[dt.tz]]]]
cmb2 <- dt.scalz[dt.tauz[dt.betaz[dt.dataz2[dt.tz]]]]
  
cmb[,u := ifelse(tz > up,up,tz)] # up is when stimulus is taken off
cmb2[,u := ifelse(tz > up2,up2,tz)]
cmb[,value := mapply(function(...){integrate(...)$value},lower=0,upper=u,scale=scalz,concentrate=con,time=tz,lag_time=tauz,degrade=betaz,mrna=V1,MoreArgs=list(f=Vectorize(function(x,scale,concentrate,time,lag_time,degrade,mrna){scale * (concentrate * (time - lag_time) - degrade * mrna)})))] # 7.9 
# if applicable (lag time * stim duration) will have to set equal to whatever function lag time is over time (ie, f(lag time v time)) 

cmb[,value := mapply(function(...){integrate(...)$value},lower=0,upper=u,scale=scalz2,concentrate=con,time=tz,lag_time=tauz2,degrade=betaz2,mrna=V1,MoreArgs=list(f=Vectorize(function(x,scale,concentrate,time,lag_time,degrade,mrna){scale * (concentrate * (time - lag_time) - degrade * mrna)})))] # 7.9 

# everything else will have to be done twice on the training

cmb <- as.data.frame(cmb)
cmb$train <- 1
cmb2 <- as.data.frame(cmb2)
cmb2$train <- 2
genez <- rbind(data,data2)

sub <- cmb[cmb$tz == genez[genez$mean == max(genez$mean),'timepoints'],]
sub2 <- cmb2[cmb2$tz == genez[genez$mean == max(genez$mean),'timepoints'],]
sub$s <- max(genez$mean) / sub$value
sub2$s <- max(genez$mean) / sub2$value

train <- rbind(cmb,cmb2)
train <- train[order(train$tauz),]
cmb$s <-  sub$s[match(cmb$tauz,sub$tauz)] 

cmb

train <- rbind(cmb,cmb2)
train}
train <- train[order(train$scalz,train$tauz,train$betaz),]
first <- train[seq(from=1, to=nrow(train), by=length(tz)*2),2:4][which.min(colSums((genez$mean - matrix(train$value, nrow=length(tz)*2))**2)),] # find minimized least squares
first}





######### 2nd run 
# if 1st run best fit is within 10 min of parameter then go 80 min in that direction, otherwise
# go 50 min in either direction for both lag and integration time window increasing the resolution now by 1 min

scalz2 <- seq(from=first$scalz - 0.00001, to=first$scalz + 0.00001, by=0.000001)
scalz2 <- scalz2[scalz2 >= 0] # only positivity
tauz2 <- seq(from=first$tauz - 10,to=first$tauz + 10,by=1)
tauz2 <- tauz2[tauz2 >= 0] # only positivity
betaz2 <- seq(from=first$betaz - 0.001,to=first$betaz + 0.001,by=0.0001)
betaz2 <- betaz2[betaz2 >= 0]
  

dt.scalz2 <- data.table(scalz2,k=1,key='k')
dt.tauz2 <- data.table(tauz2,k=1,key='k')
dt.betaz2 <- data.table(betaz2,k=1,key='k')
cmb3 <- dt.scalz2[dt.tauz2[dt.betaz2[dt.dataz[dt.tz]]]]
cmb4 <- dt.scalz2[dt.tauz2[dt.betaz2[dt.dataz[dt.tz]]]]
cmb3[,u := ifelse(tz > up,up,tz)]
cmb4[,u := ifelse(tz > up2,up2,tz)]
cmb3[,value := mapply(function(...){integrate(...)$value},lower=0,upper=u,scale=scalz2,concentrate=con,time=tz,lag_time=tauz2,degrade=betaz2,mrna=V1,MoreArgs=list(f=Vectorize(function(x,scale,concentrate,time,lag_time,degrade,mrna){scale * (concentrate * (time - lag_time) - degrade * mrna)})))] # 7.9 
cmb4[,value := mapply(function(...){integrate(...)$value},lower=0,upper=u,scale=scalz2,concentrate=con,time=tz,lag_time=tauz2,degrade=betaz2,mrna=V1,MoreArgs=list(f=Vectorize(function(x,scale,concentrate,time,lag_time,degrade,mrna){scale * (concentrate * (time - lag_time) - degrade * mrna)})))] # 7.9 



cmb3 <- as.data.frame(cmb3)
cmb4 <- as.data.frame(cmb4)
train2 <- rbind(cmb3,cmb4)
train2 <- train2[order(train2$scalz2,train2$tauz2,train2$betaz2,train2$V1),]
second <- train2[seq(from=1, to=nrow(train2), by=length(tz)*2),2:4][which.min(colSums((data$mean - matrix(train2$value, nrow=length(tz)*2))**2)),] # real - data ^2

second}

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