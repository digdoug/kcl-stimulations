######## 1min pulse 55mM bio reps 1 and 1/2  ##################################
###############################################################################
###############################################################################

 oneminpulse_c1_b1 <- data.frame(names=c("0hr_55mM_0_b1_011211", "KClpulse-1min_55mM_20_b1_011211","KClpulse-1min_55mM_40_b1_011211","KClpulse-1min_55mM_60_b1_011211","KClpulse-1min_55mM_120_b1_011211"),timepoints=c(0,20,40,60,120))  
 
   	controls1min_b1 <- getCNTLs_mod(c1,CNTLs_c1,oneminpulse_c1_b1) 

 oneminpulse_c6_b2 <- data.frame(names=c("0hr_55mM_0_b1_050412", "KClpulse-1min_55mM_20_b1_050412","KClpulse-1min_55mM_40_b1_050412","KClpulse-1min_55mM_60_b1_050412","KClpulse-1min_55mM_120_b1_050412","KClpulse-1min_55mM_240_b1_050412","KClpulse-1min_55mM_360_b1_050412"),timepoints=c(0,20,40,60,120,240,360))  

   	controls1min_b2 <- getCNTLs_mod(c6,CNTLs_c5and6,oneminpulse_c6_b2) 


######## 5min pulse 55mM bio reps 1 - 6 #######################################
###############################################################################
###############################################################################

## I was matching on three criteria, (exptname,concentration,and replicate number) there will be additional zero hours based on this code solely, I could probably write something just for the zero hours. But I want to abandon this paradigm all together


  fiveminpulse_c1_b1 <- data.frame(names=c("0hr_55mM_0_b1_011311", "KClpulse-5min_55mM_20_b1_011311","KClpulse-5min_55mM_40_b1_011311","KClpulse-5min_55mM_60_b1_011311","KClpulse-5min_55mM_120_b1_011311","KClpulse-5min_55mM_240_b1_011311","KClpulse-5min_55mM_360_b1_011311"),timepoints=c(0,20,40,60,120,240,360))  
 
   	controls1 <- getCNTLs_mod(c1,CNTLs_c1,fiveminpulse_c1_b1) 
  
  fiveminpulse_c1_b2 <- data.frame(names=c("0hr_55mM_0_b2_011411", "KClpulse-5min_55mM_20_b2_011411","KClpulse-5min_55mM_40_b2_011411","KClpulse-5min_55mM_60_b2_011411","KClpulse-5min_55mM_120_b2_011411","KClpulse-5min_55mM_240_b2_011411","KClpulse-5min_55mM_360_b2_011411"),timepoints=c(0,20,40,60,120,240,360))
  
   controls1a <- getCNTLs_mod(c1,CNTLs_c1,fiveminpulse_c1_b2)  
  
  fiveminpulse_c1_b3 <- data.frame(names=c("0hr_55mM_0_b3_011511", "KClpulse-5min_55mM_20_b3_011511","KClpulse-5min_55mM_40_b3_011511","KClpulse-5min_55mM_60_b3_011511","KClpulse-5min_55mM_120_b3_011511","KClpulse-5min_55mM_240_b3_011511","KClpulse-5min_55mM_360_b3_011511"),timepoints=c(0,20,40,60,120,240,360))      
  
  	controls1b <- getCNTLs_mod(c1,CNTLs_c1,fiveminpulse_c1_b3)
  	  	
  fiveminpulse_c2_b4 <- data.frame(names=c("0hr_55mM_0_b4_031111", "KClpulse-5min_55mM_20_b4_031111","KClpulse-5min_55mM_40_b4_031111","KClpulse-5min_55mM_60_b4_031111","KClpulse-5min_55mM_120_b4_031111","KClpulse-5min_55mM_240_b4_031111","KClpulse-5min_55mM_360_b4_031111"),timepoints=c(0,20,40,60,120,240,360))  
  
  	controls2 <- getCNTLs_mod(c2,CNTLs_chips2to4,fiveminpulse_c2_b4)
  
  
### consequently the zero hour will match w/ multiple other zero hrs because the only distinguishing factor amongst the zero hours, while still parsing out bio reps, will be the concentration
### for chips 1,2 and 4 this isn't a problem since there is essentially the same experiments on all chips, so I can probably keep this for know all those chips, but will need to create new data frames for each experiment on chips 5 and 6   
  
    fiveminpulse_c5_b5 <- data.frame(names=c("0hr_55mM_0_b1_041812", "KClpulse-5min_55mM_20_b5_041812","KClpulse-5min_55mM_40_b5_041812","KClpulse-5min_55mM_60_b5_041812","KClpulse-5min_55mM_120_b5_041812","KClpulse-5min_55mM_240_b5_041812","KClpulse-5min_55mM_360_b5_041812"),timepoints=c(0,20,40,60,120,240,360))  

  	controls3 <- getCNTLs_mod(c5,CNTLs_c5and6,fiveminpulse_c5_b5)

    fiveminpulse_c5_b6 <- data.frame(names=c("0hr_55mM_0_b1_050312", "KClpulse-5min_55mM_20_b6_050312","KClpulse-5min_55mM_40_b6_050312","KClpulse-5min_55mM_60_b6_050312","KClpulse-5min_55mM_120_b6_050312","KClpulse-5min_55mM_240_b6_050312","KClpulse-5min_55mM_360_b6_050312"),timepoints=c(0,20,40,60,120,240,360))  
    
  	controls4 <- getCNTLs_mod(c5,CNTLs_c5and6,fiveminpulse_c5_b6)    
  	
  	
######## 5min 30mM bio reps 1 - 3 #############################################
###############################################################################
###############################################################################      
      
  fiveminpulse_30mM_c5_b1 <- data.frame(names=c("0hr_55mM_0_b1_041812", "KClpulse-5min_30mM_20_b1_041812","KClpulse-5min_30mM_40_b1_041812","KClpulse-5min_30mM_60_b1_041812","KClpulse-5min_30mM_120_b1_041812","KClpulse-5min_30mM_240_b1_041812","KClpulse-5min_30mM_360_b1_041812"),timepoints=c(0,20,40,60,120,240,360))  
  
     controls_5min_c5_30mM_b1 <- getCNTLs_mod(c5,CNTLs_c5and6,fiveminpulse_30mM_c5_b1)

  fiveminpulse_30mM_c5_b2 <- data.frame(names=c("0hr_55mM_0_b1_050312", "KClpulse-5min_30mM_20_b2_050312","KClpulse-5min_30mM_40_b2_050312","KClpulse-5min_30mM_60_b2_050312","KClpulse-5min_30mM_120_b2_050312","KClpulse-5min_30mM_240_b2_050312","KClpulse-5min_30mM_360_b2_050312"),timepoints=c(0,20,40,60,120,240,360))  

     controls_5min_c5_30mM_b2 <- getCNTLs_mod(c5,CNTLs_c5and6,fiveminpulse_30mM_c5_b2)
    
  fiveminpulse_30mM_c6_b3 <- data.frame(names=c("0hr_30mM_0_b1_011811", "KClpulse-5min_30mM_20_b3_011811","KClpulse-5min_30mM_40_b3_011811","KClpulse-5min_30mM_60_b3_011811","KClpulse-5min_30mM_120_b3_011811","KClpulse-5min_30mM_240_b3_011811","KClpulse-5min_30mM_360_b3_011811"),timepoints=c(0,20,40,60,120,240,360))  

     controls_5min_c6_30mM_b3 <- getCNTLs_mod(c6,CNTLs_c5and6,fiveminpulse_30mM_c6_b3)


############ 5min 10mM bio rep 1-3 ############################################
######################################################################################       
###################################################################################### 
  	
  	  fiveminpulse_10mM_b1 <- data.frame(names=c("0hr_55mM_0_b1_041812", "KClpulse-5min_10mM_20_b1_041812","KClpulse-5min_10mM_40_b1_041812","KClpulse-5min_10mM_60_b1_041812","KClpulse-5min_10mM_120_b1_041812",'KClpulse-5min_10mM_240_b1_041812',"KClpulse-5min_10mM_360_b1_041812"),timepoints=c(0,20,40,60,120,240,360))
  
            	controls_5min_10mM_b1 <- getCNTLs_mod(c5,CNTLs_c5and6,fiveminpulse_10mM_b1)

  	  fiveminpulse_10mM_b2 <- data.frame(names=c("0hr_55mM_0_b1_050312", "KClpulse-5min_10mM_20_b2_050312","KClpulse-5min_10mM_40_b2_050312","KClpulse-5min_10mM_60_b2_050312","KClpulse-5min_10mM_120_b2_050312",'KClpulse-5min_10mM_240_b2_050312',"KClpulse-5min_10mM_360_b2_050312"),timepoints=c(0,20,40,60,120,240,360))
  	 
            	controls_5min_10mM_b2 <- getCNTLs_mod(c5,CNTLs_c5and6,fiveminpulse_10mM_b2)
  	  	 
  	  fiveminpulse_10mM_b3 <- data.frame(names=c("0hr_30mM_0_b1_011811", "KClpulse-5min_10mM_20_b3_011811","KClpulse-5min_10mM_40_b3_011811","KClpulse-5min_10mM_60_b3_011811","KClpulse-5min_10mM_120_b3_011811","KClpulse-5min_10mM_240_b3_011811","KClpulse-5min_10mM_360_b3_011811"),timepoints=c(0,20,40,60,120,240,360))
  
            	controls_5min_10mM_b3 <- getCNTLs_mod(c6,CNTLs_c5and6,fiveminpulse_10mM_b3)



########Sustained 55mM bio reps 1 - 2 (3rd included) ##########################
###############################################################################
###############################################################################

  	
 sustained_55mM_c4_b1 <- data.frame(names=c("0hr_55mM_0_b1_032411", "KClsustain_55mM_20_b1_032411","KClsustain_55mM_40_b1_032411","KClsustain_55mM_60_b1_032411","KClsustain_55mM_120_b1_032411","KClsustain_55mM_240_b1_032411","KClsustain_55mM_360_b1_032411"),timepoints=c(0,20,40,60,120,240,360))

    controls_sustained_b1 <- getCNTLs_mod(c4,CNTLs_chips2to4, sustained_55mM_c4_b1)

 sustained_55mM_c4_b2 <- data.frame(names=c("0hr_55mM_0_b2_033111", "KClsustain_55mM_20_b2_033111","KClsustain_55mM_40_b2_033111","KClsustain_55mM_60_b2_033111","KClsustain_55mM_120_b2_033111","KClsustain_55mM_240_b2_033111","KClsustain_55mM_360_b2_033111"),timepoints=c(0,20,40,60,120,240,360))    

    controls_sustained_b2 <- getCNTLs_mod(c4,CNTLs_chips2to4, sustained_55mM_c4_b2)
    
  sustained_55mM_c5_lt <- data.frame(names=c("0hr_55mM_0_b1_050312", "KClsustain_55mM_60_b1_050312","KClsustain_55mM_180_b1_050312","KClsustain_55mM_360_b1_050312"),timepoints=c(0,60,180,360))

    controls_sustained_c5_lt <- getCNTLs_mod(c5,CNTLs_c5and6, sustained_55mM_c5_lt)


########Sustained 3mM bio reps 1 - 2 (3rd included) ###########################
###############################################################################
###############################################################################

   
  	  sustained_3mM_c5_b1 <- data.frame(names=c("0hr_55mM_0_b1_050312", "KClsustain_3mM_60_b1_050312","KClsustain_3mM_180_b1_050312","KClsustain_3mM_360_b1_050312"),timepoints=c(0,60,180,360))

      controls_sustained_c5_3mM_b1 <- getCNTLs_mod(c5,CNTLs_c5and6, sustained_3mM_c5_b1)

  	  sustained_3mM_c6_b2 <- data.frame(names=c("0hr_10mM_0_b1_041112", "KClsustain_3mM_60_b1_041112","KClsustain_3mM_180_b1_041112","KClsustain_3mM_360_b1_041112"),timepoints=c(0,60,180,360))
  	
      controls_sustained_c6_3mM_b2 <- getCNTLs_mod(c6,CNTLs_c5and6, sustained_3mM_c6_b2)

  	  sustained_3mM_c6_moretime <- data.frame(names=c("0hr_3mM_0_b1_051012", "KClsustain_3mM_20_b1_051012","KClsustain_3mM_40_b1_051012","KClsustain_3mM_60_b1_051012",'KClsustain_3mM_120_b1_051012','KClsustain_3mM_240_b1_051012','KClsustain_3mM_360_b1_051012'),timepoints=c(0,20,40,60,120,240,360))
  	
      controls_sustained_c6_3mM_mt <- getCNTLs_mod(c6,CNTLs_c5and6, sustained_3mM_c6_moretime)
      
######## Sustained 10mM bio reps 1 - 2 ########################################
###############################################################################
###############################################################################      
      
      sustained_10mM_c5_b1 <- data.frame(names=c("0hr_55mM_0_b1_050312", "KClsustain_10mM_60_b1_050312","KClsustain_10mM_180_b1_050312","KClsustain_10mM_360_b1_050312"),timepoints=c(0,60,180,360))

      controls_sustained_c5_10mM_b1 <- getCNTLs_mod(c5,CNTLs_c5and6, sustained_10mM_c5_b1)

  	  sustained_10mM_c6_b2 <- data.frame(names=c("0hr_10mM_0_b1_041112", "KClsustain_10mM_60_b1_041112","KClsustain_10mM_180_b1_041112","KClsustain_10mM_360_b1_041112"),timepoints=c(0,60,180,360))
  	
      controls_sustained_c6_10mM_b2 <- getCNTLs_mod(c6,CNTLs_c5and6, sustained_10mM_c6_b2)
      

############ sustain 30 min @ 55mM  ##################################################
######################################################################################       
###################################################################################### 
  	
  	  sustain_30mM_b1 <- data.frame(names=c("0hr_55mM_0_b1_050312","KClsustain_30mM_60_b1_050312","KClsustain_30mM_180_b1_050312","KClsustain_30mM_360_b1_050312"),timepoints=c(0,60,180,360))

   controls_sustain_30mM_b1 <- getCNTLs_mod(c5,CNTLs_c5and6,sustain_30mM_b1)
      


############# 15 min pulse bio rep 1 #################################################
######################################################################################
######################################################################################
     
   fifteenminpulse_55mM_b1 <- data.frame(names=c("0hr_55mM_0_b1_050312", "KClpulse-15min_55mM_20_b1_050312","KClpulse-15min_55mM_40_b1_050312","KClpulse-15min_55mM_60_b1_050312","KClpulse-15min_55mM_120_b1_050312","KClpulse-15min_55mM_240_b1_050312","KClpulse-15min_55mM_360_b1_050312"),timepoints=c(0,20,40,60,120,240,360)) 

       controls_15min <- getCNTLs_mod(c6,CNTLs_c5and6,fifteenminpulse_55mM_b1)


############ 30 min pulse bio rep 1 - 2 ##############################################
######################################################################################       
######################################################################################       
c6[c6$primer == 'Rasd1' & c6$exptname == 'KClpulse-15min',]

  thirtyminpulse_55mM_b1 <- data.frame(names=c("0hr_55mM_0_b1_040711", "KClpulse-30min_55mM_20_b1_040711","KClpulse-30min_55mM_40_b1_040711","KClpulse-30min_55mM_60_b1_040711","KClpulse-30min_55mM_120_b1_040711","KClpulse-30min_55mM_240_b1_040711","KClpulse-30min_55mM_360_b1_040711"),timepoints=c(0,20,40,60,120,240,360)) 

       controls_30min_b1 <- getCNTLs_mod(c5,CNTLs_c5and6,thirtyminpulse_55mM_b1)
       
  thirtyminpulse_55mM_b2 <- data.frame(names=c("0hr_55mM_0_b1_050412", "KClpulse-30min_55mM_20_b2_050412","KClpulse-30min_55mM_40_b2_050412","KClpulse-30min_55mM_60_b2_050412","KClpulse-30min_55mM_120_b2_050412","KClpulse-30min_55mM_240_b2_050412","KClpulse-30min_55mM_360_b2_050412"),timepoints=c(0,20,40,60,120,240,360)) 

       controls_30min_b2 <- getCNTLs_mod(c6,CNTLs_c5and6,thirtyminpulse_55mM_b2)

############ 120 min pulse bio rep 1 - 2 #############################################
######################################################################################       
######################################################################################  

   twohrpulse_55mM_b1 <- data.frame(names=c("0hr_55mM_0_b1_040711", "KClpulse-120min_55mM_20_b1_040711","KClpulse-120min_55mM_40_b1_040711","KClpulse-120min_55mM_60_b1_040711","KClpulse-120min_55mM_120_b1_040711","KClpulse-120min_55mM_240_b1_040711","KClpulse-120min_55mM_360_b1_040711"),timepoints=c(0,20,40,60,120,240,360)) 

       controls_120min_b1 <- getCNTLs_mod(c5,CNTLs_c5and6,twohrpulse_55mM_b1)
       
   twohrpulse_55mM_b2 <- data.frame(names=c("0hr_55mM_0_b1_050412", "KClpulse-120min_55mM_20_b2_050412","KClpulse-120min_55mM_40_b2_050412","KClpulse-120min_55mM_60_b2_050412","KClpulse-120min_55mM_120_b2_050412","KClpulse-120min_55mM_240_b2_050412","KClpulse-120min_55mM_360_b2_050412"),timepoints=c(0,20,40,60,120,240,360)) 

       controls_120min_b2 <- getCNTLs_mod(c6,CNTLs_c5and6,twohrpulse_55mM_b2)

############ 60 min pulse bio rep 1 -2 ###############################################
######################################################################################       
######################################################################################  

  	   sixtyminpulse_55mM_b1 <- data.frame(names=c("0hr_55mM_0_b1_050312", "KClpulse-60min_55mM_20_b1_050312","KClpulse-60min_55mM_40_b1_050312","KClpulse-60min_55mM_60_b1_050312","KClpulse-60min_55mM_120_b1_050312","KClpulse-60min_55mM_240_b1_050312","KClpulse-60min_55mM_360_b1_050312"),timepoints=c(0,20,40,60,120,240,360)) 

       controls_60min_b1 <- getCNTLs_mod(c6,CNTLs_c5and6,sixtyminpulse_55mM_b1)

  	   sixtyminpulse_55mM_b2 <- data.frame(names=c("0hr_55mM_0_b1_050412", "KClpulse-60min_55mM_20_b2_050412","KClpulse-60min_55mM_40_b2_050412","KClpulse-60min_55mM_60_b2_050412","KClpulse-60min_55mM_120_b2_050412","KClpulse-60min_55mM_240_b2_050412","KClpulse-60min_55mM_360_b2_050412"),timepoints=c(0,20,40,60,120,240,360)) 

       controls_60min_b2 <- getCNTLs_mod(c6,CNTLs_c5and6,sixtyminpulse_55mM_b2)


############ CHX/ANS bio rep 1 -2 ####################################################
######################################################################################       
######################################################################################  

  	   chx_55mM_b1 <- data.frame(names=c("0hr-CHX_55mM_0_b1_032411", "KClsustain-CHX_55mM_20_b1_032411","KClsustain-CHX_55mM_40_b1_032411","KClsustain-CHX_55mM_60_b1_032411","KClsustain-CHX_55mM_120_b1_032411","KClsustain-CHX_55mM_240_b1_032411","KClsustain-CHX_55mM_360_b1_032411"),timepoints=c(0,20,40,60,120,240,360)) 

       controls_chx_b1 <- getCNTLs_mod(c4_chx,CNTLs_chips2to4,chx_55mM_b1)

  	   chx_55mM_b2 <- data.frame(names=c("0hr-CHX_55mM_0_b2_033111", "KClsustain-CHX_55mM_20_b2_033111","KClsustain-CHX_55mM_40_b2_033111","KClsustain-CHX_55mM_60_b2_033111","KClsustain-CHX_55mM_120_b2_033111","KClsustain-CHX_55mM_240_b2_033111","KClsustain-CHX_55mM_360_b2_033111"),timepoints=c(0,20,40,60,120,240,360)) 

       controls_chx_b2 <- getCNTLs_mod(c4_chx,CNTLs_chips2to4,chx_55mM_b2)
       
       
############ 2 pulses rep 1 -2 ####################################################
######################################################################################       
######################################################################################         


  	   kcl_2pulse_55mM_b1 <- data.frame(names=c("0hr-2pulse_55mM_0_b1_030711", "KCl2pulse-120min_55mM_20_b1_030711","KCl2pulse-120min_55mM_40_b1_030711","KCl2pulse-120min_55mM_60_b1_030711","KCl2pulse-120min_55mM_80_b1_030711","KCl2pulse-120min_55mM_100_b1_030711","KCl2pulse-120min_55mM_120_b1_030711","KCl2pulse-120min_55mM_140_b1_030711","KCl2pulse-120min_55mM_160_b1_030711",'KCl2pulse-120min_55mM_180_b1_030711','KCl2pulse-120min_55mM_240_b1_030711','KCl2pulse-120min_55mM_360_b1_030711','KCl2pulse-120min_55mM_480_b1_030711'),timepoints=c(0,20,40,60,80,100,120,140,160,180,240,360,480)) 

       controls_2pulse_b1 <- getCNTLs_mod(doub_pulse,CNTLs_chips2to4,kcl_2pulse_55mM_b1)
       
  	   kcl_2pulse_55mM_b2 <- data.frame(names=c("0hr-2pulse_55mM_0_b2_030811", "KCl2pulse-120min_55mM_20_b2_030811","KCl2pulse-120min_55mM_40_b2_030811","KCl2pulse-120min_55mM_60_b2_030811","KCl2pulse-120min_55mM_80_b2_030811","KCl2pulse-120min_55mM_100_b2_030811","KCl2pulse-120min_55mM_120_b2_030811","KCl2pulse-120min_55mM_140_b2_030811","KCl2pulse-120min_55mM_160_b2_030811",'KCl2pulse-120min_55mM_180_b2_030811','KCl2pulse-120min_55mM_240_b2_030811','KCl2pulse-120min_55mM_360_b2_030811','KCl2pulse-120min_55mM_480_b2_030811'),timepoints=c(0,20,40,60,80,100,120,140,160,180,240,360,480)) 
       
       controls_2pulse_b2 <- getCNTLs_mod(doub_pulse,CNTLs_chips2to4,kcl_2pulse_55mM_b2)
       