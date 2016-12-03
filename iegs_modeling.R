# 5 min 

fos_5min <- all_5min[all_5min$primer %in% 'Fos',]
egr4_5min <- all_5min[all_5min$primer %in% 'Egr4',]
npas4_5min <- all_5min[all_5min$primer %in% 'Npas4',]
junb_5min <- all_5min[all_5min$primer %in% 'Junb',] 

# sus

fos <- all_sus[all_sus$primer %in% 'Fos',]
egr4 <- all_sus[all_sus$primer %in% 'Egr4',]
npas4 <- all_sus[all_sus$primer %in% 'Npas4',]
junb <- all_sus[all_sus$primer %in% 'Junb',]

fos_5min
fos


# ieg equation
rt <- function(w,t){
	w * exp(-d * t) * integrate(exp(d * tau) * function(...),lower=0,upper=t)  
}

yo <- deriv(~ a^2,"a")
class(yo)

yo[.grad]

install.packages('stats')
library('stats')

fit <- loess(mean ~ timepoints,data=fos)
fit.points <- predict(fit,newdata=data.frame(speed=seq(min(timepoints),max(timepoints),length=100)),se=F)
fitdf <- data.frame(x=seq(min(timepoints),max(timepoints), length=100),y=fit.points) 	
 	predict.loess()
 	    ggplot(fos,aes(x=timepoints,y= mean,ymax = mean + sde, ymin = mean  - sde)) + 
           geom_errorbar(width=2) +
           geom_point() +
           geom_line() +
           stat_smooth(method='loess',se=F) + 
           xlab('Time (min)') +
           ylab('Fold Induction') +          
           opts(title = 'yo')   
           
