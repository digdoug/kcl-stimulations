rna <- read.delim('RNA-seq/RNA-seq.txt')
rna <- rna[rna$Gene %in% lrgs,]

rna$intden_0 <- (rna$INT_Density.Rep1.0hr + rna$INT_Density.Rep2.0hr) / 2
rna$intden_1 <- (rna$INT_Density.Rep1.1hr + rna$INT_Density.Rep2.1hr) / 2
rna$intden_6 <- (rna$INT_Density.Rep1.6hr + rna$INT_Density.Rep2.6hr) / 2

rna$exden_0 <- (rna$EXN_Density.Rep1.0hr + rna$EXN_Density.Rep2.0hr) / 2
rna$exden_1 <- (rna$EXN_Density.Rep1.1hr + rna$EXN_Density.Rep2.1hr) / 2
rna$exden_6 <- (rna$EXN_Density.Rep1.6hr + rna$EXN_Density.Rep2.6hr) / 2


# INT density
for(i in 1:nrow(rna)){
rna$intden_0_sde[i] <- std.error(c(rna$INT_Density.Rep1.0hr[i],rna$INT_Density.Rep2.0hr[i]))
	}
	
for(i in 1:nrow(rna)){
rna$intden_1_sde[i] <- std.error(c(rna$INT_Density.Rep1.1hr[i],rna$INT_Density.Rep2.1hr[i]))
	}

for(i in 1:nrow(rna)){
rna$intden_6_sde[i] <- std.error(c(rna$INT_Density.Rep1.6hr[i],rna$INT_Density.Rep2.6hr[i]))
	}


# EXN density
for(i in 1:nrow(rna)){
rna$exden_0_sde[i] <- std.error(c(rna$EXN_Density.Rep1.0hr[i],rna$EXN_Density.Rep2.0hr[i]))
	}
	
for(i in 1:nrow(rna)){
rna$exden_1_sde[i] <- std.error(c(rna$EXN_Density.Rep1.1hr[i],rna$EXN_Density.Rep2.1hr[i]))
	}

for(i in 1:nrow(rna)){
rna$exden_6_sde[i] <- std.error(c(rna$EXN_Density.Rep1.6hr[i],rna$EXN_Density.Rep2.6hr[i]))
	}

## INTRON EXON density


in0 <- rna[,c('Gene','intden_0','intden_0_sde')]
in0$num <- 0
colnames(in0) <- c('gene','mean','sde','num')
in1 <- rna[,c('Gene','intden_1','intden_1_sde')]
in1$num <- 60
colnames(in1) <- c('gene','mean','sde','num')
in6 <- rna[,c('Gene','intden_6','intden_6_sde')]
in6$num <- 360
colnames(in6) <- c('gene','mean','sde','num')


ex0 <- rna[,c('Gene','exden_0','exden_0_sde')]
ex0$num <- 0
colnames(ex0) <- c('gene','mean','sde','num')
ex1 <- rna[,c('Gene','exden_1','exden_1_sde')]
ex1$num <- 60
colnames(ex1) <- c('gene','mean','sde','num')
ex6 <- rna[,c('Gene','exden_6','exden_6_sde')]
ex6$num <- 360
colnames(ex6) <- c('gene','mean','sde','num')


intron <- rbind(in0,in1,in6)
intron$type <- 'intron'
exon <- rbind(ex0,ex1,ex6)
exon$type <- 'exon'
intro <- ddply(intron, .(gene), summarise, mean = max(mean))
intro$sde <- intron$sde[match(intro$mean,intron$mean)]
intro$time <- intron$num[match(intro$mean,intron$mean)]
intro_max <- intro[intro$gene != c('Ier2','Junb'),]
intro_last <- intron[intron$gene %in% unique(intro_max$gene) & intron$num == 360,]
intro_last <- intro_last[,-5]
colnames(intro_last)[4] <- 'time'
intron_density <- rbind(intro_max,intro_last)


intron_density <- rbind(intro_max,intro_last)

new <- ddply(intron,.(gene), transform, sde=sde/max(mean),mean=mean/max(mean))
new1 <- ddply(exon, .(gene), transform, sde=sde/max(mean),mean=mean/max(mean))

comb <- rbind(new,new1)
y <- ddply(new, .(gene), summarise, con=!is.na(any(mean)))
y1 <- ddply(new1, .(gene), summarise, con=!is.na(any(mean)))
yii <- rbind(y,y1)
check <- ddply(yii, .(gene), summarise, which(!any(con==F)))
yo <- setdiff(all_primers,check$gene)
new_genes <- setdiff(all_primers,yo)
comb[order(comb$gene),]
comb <- comb[comb$gene %in% new_genes,]

normed <- ddply(intron_density, .(gene), transform, norm = mean / max(mean))


rna_plot <- function(prime,x){	
        x <- x[x$gene %in% prime,]
 	    ggplot(x,aes(x=num,y= mean, ymax= mean-sde,ymin=mean+sde,colour=type)) +
 	       geom_errorbar(width=2) +  
           geom_point() +
           geom_line() +
           xlab('Time (min)') +
           ylab('Intron/Exon Density') +          
           opts(title = prime)             
}


    ggplot(intron,aes(x=num,y=mean,colour=gene)) + 
           geom_point() +
           geom_line() +
           # scale_x_continuous(breaks=c(60,360),label=c('','Final')) +
           xlab('Time Point') +
           ylab('Intron Density') +          
           opts(title = 'Transcription Increase')   

yo <- comb[comb$gene %in% 'Egr1',]

rna_plot('Egr1',yo)

plotAllrnas(new_genes,comb,'In-ex')





