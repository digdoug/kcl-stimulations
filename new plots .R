### 1min @ 55min

           c1_1min_b1 <- normalize_con(primerSet_c1,c1,oneminpulse_c1_b1,controls1min_b1)
           c6_1min_b2 <- normalize_con(primerSet_c5and6,c6,oneminpulse_c6_b2,controls1min_b2)

### 5min @ 55mM            
           c1_5min_b1 <-  normalize_con(primerSet_c1,c1,fiveminpulse_c1_b1,controls1)
           c1_5min_b2 <-  normalize_con(primerSet_c1,c1,fiveminpulse_c1_b2,controls1a)
           c1_5min_b3 <-  normalize_con(primerSet_c1,c1,fiveminpulse_c1_b3,controls1b)           
           c2_5min_b4 <-  normalize_con(primerSet_c2to4,c2,fiveminpulse_c2_b4,controls2)
           c5_5min_b5 <-  normalize_con(primerSet_c5and6,c5,fiveminpulse_c5_b5,controls3)
           c5_5min_b6 <-  normalize_con(primerSet_c5and6,c5,fiveminpulse_c5_b6,controls4)
                    
### 5min @ 30mM
		   c5_5min_30mM_b1 <- normalize_con(primerSet_c5and6,c5,fiveminpulse_30mM_c5_b1,controls_5min_c5_30mM_b1)
		   c5_5min_30mM_b2 <- normalize_con(primerSet_c5and6,c5,fiveminpulse_30mM_c5_b2,controls_5min_c5_30mM_b2)		   
		   c6_5min_30mM_b3 <- normalize_con(primerSet_c5and6,c6,fiveminpulse_30mM_c6_b3,controls_5min_c6_30mM_b3)

### 5min @ 10mM 
            c5_5min_10mM_b1 <- normalize_con(primerSet_c5and6,c5,fiveminpulse_10mM_b1,controls_5min_10mM_b1)
            c5_5min_10mM_b2 <- normalize_con(primerSet_c5and6,c5,fiveminpulse_10mM_b2,controls_5min_10mM_b2)            
            c6_5min_10mM_b3 <- normalize_con(primerSet_c5and6,c6,fiveminpulse_10mM_b3,controls_5min_10mM_b3)
           
### sustained @ 55mM
           c4_sustained_b1 <-  normalize_con(primerSet_c2to4,c4,sustained_55mM_c4_b1,controls_sustained_b1)
           c4_sustained_b2 <-  normalize_con(primerSet_c2to4,c4,sustained_55mM_c4_b2,controls_sustained_b2)
           c5_sustained_b1_lt <- normalize_con(primerSet_c5and6,c5,sustained_55mM_c5_lt,controls_sustained_c5_lt)

### sustained @ 3mM
           c5_3mMsustained_b1 <- normalize_con(primerSet_c5and6,c5,sustained_3mM_c5_b1,controls_sustained_c5_3mM_b1)
           c6_3mMsustained_b2 <- normalize_con(primerSet_c5and6,c6,sustained_3mM_c6_b2,controls_sustained_c6_3mM_b2)   
           c6_3mMsustained_b1_mt <- normalize_con(primerSet_c5and6,c6,sustained_3mM_c6_moretime,controls_sustained_c6_3mM_mt)
                     
### sustained @ 10mM           
           c5_10mMsustained_b1 <- normalize_con(primerSet_c5and6,c5,sustained_10mM_c5_b1,controls_sustained_c5_10mM_b1) 
           c6_10mMsustained_b2 <- normalize_con(primerSet_c5and6,c6,sustained_10mM_c6_b2,controls_sustained_c6_10mM_b2)

### sustained @ 30mM          
           c5_30mMsustained_b1 <- normalize_con(primerSet_c5and6,c5,sustain_30mM_b1,controls_sustain_30mM_b1)
           all_30mM <- means(c5_30mMsustained_b1)
           primes_sustain_30mM <- all_primers[!is.na(match(all_primers,unique(c5_30mMsustained_b1$primer)))] 

### 15min @ 55mM                  
           c6_15min_b1 <- normalize_con(primerSet_c5and6,c6,fifteenminpulse_55mM_b1,controls_15min)
           all_15min <- means(c6_15min_b1)
           primes_15min <- all_primers[!is.na(match(all_primers,unique(c6_15min_b1$primer)))] 
            
### 30min @ 55mM          
           c5_30min_b1 <- normalize_con(primerSet_c5and6,c5,thirtyminpulse_55mM_b1,controls_30min_b1)        
           c6_30min_b2 <- normalize_con(primerSet_c5and6,c6,thirtyminpulse_55mM_b2,controls_30min_b2)
                     
### 120min @ 55mM           
           c5_120min_b1 <- normalize_con(primerSet_c5and6,c5,twohrpulse_55mM_b1,controls_120min_b1)       
           c6_120min_b2 <- normalize_con(primerSet_c5and6,c6,twohrpulse_55mM_b2,controls_120min_b2)
           
### 60min @ 55mM           
           c6_60min_b1 <- normalize_con(primerSet_c5and6,c6,sixtyminpulse_55mM_b1,controls_60min_b1)
           c6_60min_b2 <- normalize_con(primerSet_c5and6,c6,sixtyminpulse_55mM_b2,controls_60min_b2)
           
### chx chip 4
           c4_chx_b1 <- normalize_con(primerSet_c2to4,c4_chx,chx_55mM_b1,controls_chx_b1)           
           c4_chx_b2 <- normalize_con(primerSet_c2to4,c4_chx,chx_55mM_b2,controls_chx_b2)           

### 2 pulse @ 55mM 
           twopulse_55mM_b1 <- normalize_con(primerSet_c2to4,doub_pulse,kcl_2pulse_55mM_b1,controls_2pulse_b1)           
           twopulse_55mM_b2 <- normalize_con(primerSet_c2to4,doub_pulse,kcl_2pulse_55mM_b2,controls_2pulse_b2)   
    
    
           all_2pulse <- rbind(twopulse_55mM_b1,twopulse_55mM_b2)
           primes_2pulse <- all_primers[!is.na(match(all_primers,unique(all_2pulse$primer)))]
    
           all_chx <- rbind(c4_chx_b1,c4_chx_b2)
           all_chx$type <- 'KCl + CHX'
           primes_chx <- all_primers[!is.na(match(all_primers,unique(all_chx$primer)))]

           all_1min <- rbind(c1_1min_b1,c6_1min_b2)
           primes_1min <- all_primers[!is.na(match(all_primers,unique(all_1min$primer)))]
        
           all_5min <- rbind(c1_5min_b1,c1_5min_b2,c1_5min_b3,c2_5min_b4,c5_5min_b5,c5_5min_b6)
           primes_5min <- all_primers[!is.na(match(all_primers,unique(all_5min$primer)))] 
        
           all_5min_30mM <- rbind(c5_5min_30mM_b1,c5_5min_30mM_b2,c6_5min_30mM_b3)
           primes_5min_30mM <- all_primers[!is.na(match(all_primers,unique(all_5min_30mM$primer)))] 

           all_5min_10mM <- rbind(c5_5min_10mM_b1,c5_5min_10mM_b2,c6_5min_10mM_b3)
           primes_5min_10mM <- all_primers[!is.na(match(all_primers,unique(all_5min_10mM$primer)))] 

           all_sustain_55mM <- means(rbind(c4_sustained_b1,c4_sustained_b2,c5_sustained_b1_lt))
           all_sustain_55mM$type <- 'KCl'
           all_sus <- all_sustain_55mM[all_sustain_55mM$timepoints != 180,]         
           all_sus <- all_sus[all_sus$primer %in% primes_chx,]
           primes_sustain <- all_primers[!is.na(match(all_primers,unique(all_sustain_55mM$primer)))] 

           all_sustain_3mM <- rbind(c5_3mMsustained_b1,c6_3mMsustained_b2,c6_3mMsustained_b1_mt)
           primes_sustain_3mM <- all_primers[!is.na(match(all_primers,unique(all_sustain_3mM$primer)))] 

           all_sustain_10mM <- rbind(c5_10mMsustained_b1,c6_10mMsustained_b2)
           primes_sustain_10mM <- all_primers[!is.na(match(all_primers,unique(all_sustain_10mM$primer)))] 

           all_30min <- rbind(c5_30min_b1,c6_30min_b2)
           primes_30min <- all_primers[!is.na(match(all_primers,unique(all_30min$primer)))] 

           all_120min <- rbind(c5_120min_b1,c6_120min_b2)
           primes_120min <- all_primers[!is.na(match(all_primers,unique(all_120min$primer)))] 
           
           all_60min <- rbind(c6_60min_b1,c6_60min_b2)         
           primes_60min <- all_primers[!is.na(match(all_primers,unique(all_60min$primer)))] 
 

           kcl_chx_df <- rbind(all_sus,all_chx)
           
           d_prgs <- all_chx[all_chx$primer %in% c('Rasd1','Pax1','Pcsk1','Crh','Bdnf'),]
           d_prgs_max <- d_prgs[d_prgs$timepoints == 240,] 
           d_prgs_max$timepoints <- 0
           d_prgs_fin <- d_prgs[d_prgs$timepoints == 360,]
           d_prgs_fin$timepoints <- 50
           
           dprgs <- rbind(d_prgs_max,d_prgs_fin)
           dprgs
           test <- ddply(dprgs,.(primer),transform, fi = mean/min(mean),fi_sde = sde/min(mean))
           test 
           
           ggplot(test,aes(x=timepoints,y= fi,colour=primer)) + 
           # geom_errorbar(width=.2) + 
           geom_point() +
           geom_line() +
           xlab('Time (min)') +
           ylab('Fold Induction') +          
           opts(title = 'dPRG Super Induction') 

###################
########
#### normal plots for all our data (all bio reps included)
########
###################
#### normal plots were made using this



plotAll2pulse(primes_2pulse,all_2pulse,'2pulse')

plotchx_gg(primes_chx,kcl_chx_df,'CHX')

plotAllPrimers_gg(primes_1min,all_1min,'1min')

plotAllPrimers_gg(primes_5min,all_5min,'5min')

plotAllPrimers_gg(primes_5min_30mM,all_5min_30mM,'5min_30mM')

plotAllPrimers_gg(primes_5min_10mM,all_5min_10mM,'5min_10mM')

plotAllPrimers_gg(primes_sustain,all_sus,'test')

plotAllPrimers_gg(primes_sustain_3mM,all_sustain_3mM,'sustain_3mM')

plotAllPrimers_gg(primes_sustain_10mM,all_sustain_10mM,'sustain_10mM')

plotAllPrimers_gg(primes_sustain_30mM,c5_30mMsustained_b1,'sustain_30mM')

plotAllPrimers_gg(primes_15min,c6_15min_b1,'15min')

plotAllPrimers_gg(primes_30min,all_30min,'30min')

plotAllPrimers_gg(primes_120min,all_120min,'120min')

plotAllPrimers_gg(primes_60min,all_60min,'60min')

plotBIG_by2(primez,all_1min,all_5min,all_15min,all_30min,all_60min,all_120min,all_sus,'BIG')

cel <- setdiff(all_primers,unique(all_1min$primer))
primez <- setdiff(all_primers,cel)
yo1 <- setdiff(primez,unique(all_15min$primer))
primez <- setdiff(primez,yo1)
yo1 <- setdiff()
yo2 <- setdiff(primes_big,unique(all_sus$primer))


################################################################################
################################################################################
###maxes for differnt pulse widths #############################################

max1min <- maxes(all_1min)
max1min$kcl_duration <- 1

max5min <- maxes(all_5min)
max5min$kcl_duration <- 5

max15min <- maxes(c6_15min_b1)
max15min$kcl_duration <- 15

max30min <- maxes(all_30min)
max30min$kcl_duration <- 30

max60min <- maxes(all_60min)
max60min$kcl_duration <- 60

max120min <- maxes(all_120min)
max120min$kcl_duration <- 120

max360min <- maxes(all_sustain_55mM)
max360min$kcl_duration <- 360

pulses <-  rbind(max1min,max5min,max15min,max30min,max60min,max120min,max360min)
 
##### this is just to make it easier to look at genes
 
 pulses <- pulses[order(pulses$primer),]
 
 ### find out which filters don't pass the cycle threshold for any zero hours
 setdiff(primerSet_c5and6,unique(pulses$primer))
 primes_pulses <- primerSet_c5and6[primerSet_c5and6 != 'Fosb']
 
plotAbs_pulse(primes_pulses,pulses,'Max_FI_test') 

###############################################################################################
###############################################################################################
##### maxes for different concentrations on 5 min pulse
## all of these were run on the last two chips, so the primer vector will contain ~90 primers depending on what is filtered
###############################################################################################


#### in this next part you could recieve a 'no non-missing arguments to max; returning -Inf' error message
#### if this is the case the filter didn't work and non-existent values go included in the array (eg 'Numeric')

      maxes_5min_10mM <- maxes(all_5min_10mM)
      maxes_5min_10mM$conc <- 10
      
      maxes_5min_30mM <- maxes(all_5min_30mM)
      maxes_5min_30mM$conc <- 30
      
      maxes_5min_55mM <- maxes(all_5min)
      maxes_5min_55mM$conc <- 55
      
      conc_5min_maxes <- rbind(maxes_5min_10mM,maxes_5min_30mM,maxes_5min_55mM)
      
      conc_5min_maxes <- conc_5min_maxes[order(conc_5min_maxes$primer),]

      plotAbs_conc(primes_5min_30mM,conc_5min_maxes,'5min_diff_concs')
      
      
################################################################################################
################################################################################################
##### maxes for different concentrations on a sustained 6 hour experiment 

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

