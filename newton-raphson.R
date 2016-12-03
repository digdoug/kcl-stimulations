library('rootSolve')
library('maxLik')





loglik <- function(param){
	
	scale <- param[1]
	lag_time <- param[2]
	degrade <- param[3]
	con <- 55
	time <- 360
	
    integ <- Vectorize(function(x){scale * (con * (time - lag_time)) - degrade * 34.5735984})

           inin <-  integrate(integ,lower=0,upper=time)
inin$value
}



gaglik <- function(param){
	
	scale <- param[1]
	lag_time <- param[2]
	degrade <- param[3]
	con <- 55 
	time <- 60 
	
	integ <- Vectorize(function(x){scale * (con * (time - lag_time)) - degrade * })
	
	inin <- integrate(integ,lower=0,upper=time)
	
	inin$value
}

try <- c(0.0001,50,0.01)
x <- rnorm(1000,1,2)
N <- length(xgoo)
res <- maxNR(loglik, start=c(0.0001,40,0.0001))
loglik
?maxNR
summary(res)


 ?maxNR
function(...){integrate(...)$value},lower=0,upper=u,scale=scalz,concen=con,time=tz,lag_time=tauz,degrade=betaz,mrna=V1,MoreArgs=list(f=Vectorize(function(x,scale,concentrate,time,lag_time,degrade,mrna){scale * (concen * (time - lag_time)) - degrade * mrna}))

maxNR(,)


?maxNR








