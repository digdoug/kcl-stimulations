#####The only thing that may be different for each chip is their genes, depending on the size, or Jesse's selection
#####REq-funcitons_supermod.R should be in memory before reading any of this. All of these raw data files have been modified, and should all contain the modified_lucas.R tag on all of them

#### upload these packages, they make life worth living

library('ggplot2') 
library('reshape2')
library('plyr')
library('grid')
library('data.table')
library('robustbase')
library('chron')
library('plotrix')
library('micEcon')

#####Set working Directory, this is mine  

   setwd('/Users/ta/Documents/R/Fluidigm_data_and_sample_sub_sheets')
  
#####Define genes to respective classes
### part of the project is to be able to identify these

   IEGs_c1 <- c("Fosb","Npas4","Nr4a1","Junb","Fos","Egr2","Egr4","Pcsk1","Egr1","Gadd45b","Maff","Ier2","Klf4","Vgf","Arc","Atf7ip2","Nr4a2","Egr3")
   LRGs_c1 <- c("Col10a1","Gprc5a","Bves","Pax1","Bdnf","Pdlim3","Acan","Npffr2","Crh","Rasd1","Serpinb2","Ankrd56","Corin","Pthlh","Fosl2")
   CNTLs_c1 <- c("Casp6","Prkcn","Ssna1","Hs2st1","Mrpl50","Nnmt","Mllt11","Eif4h","Atp2c1","Tubb3","Actb","Gapdh")
   primerSet_c1 <- c(IEGs_c1,LRGs_c1,CNTLs_c1)
        
   
   IEGs_chips2to4 <- c("Fosb","Npas4","Nr4a1","Junb","Fos","Egr2","Egr4","Pcsk1","Egr1","Gadd45b","Maff","Ier2","Klf4","Vgf","Arc","Atf7ip2","Nr4a2","Egr3")
   LRGs_chips2to4 <- c("Col10a1","Gprc5a","Bves","Pax1","Bdnf","Pdlim3","Acan","Npffr2","Crh","Rasd1","Serpinb2","Ankrd56","Corin","Pthlh","Fosl2","Fosl1")
   CNTLs_chips2to4 <- c("Casp6","Prkcn","Ssna1","Hs2st1","Mrpl50","Nnmt","Mllt11","AI597468","Eif4h","Atp2c1","Tubb3","Actb","Gapdh")
   primerSet_c2to4 <- c(IEGs_chips2to4,LRGs_chips2to4,CNTLs_chips2to4) 

###I'm going to use the same gene classes for chip 6 as were used for chip 5   
###Just because they are the same size, not a good excuse, should ask Jesse about this


  IEGs_c5and6 <- c("Fosb","Egr1","Egr2","Egr3","Egr4","Nr4a1","Nr4a2","Nr4a3","Arc","Atf7ip2","Npas4","Gadd45b","Junb","Ier2","Klf4","Fos","Maff","Cyr61","Ccno","Lhx5","Amigo3","Dusp1","Fbxo33","Gltscr2")
   LRGs_c5and6 <- c("Pdlim3","Nnmt","Acan","Fosl2","Pthlh","Csrnp1","Trib1","Popdc3","Sertad1","Dusp5","Wif1","Crhbp","Pnoc","Nppc","Npffr2","Rel","Vgf","Bdnf","Corin","Pam","Rbpms","Vamp4","Tnn","Crh","Mlf1","Sik1","Penk","Gprc5a","Rasd1","Ptger4","Cartpt","Hspb3","Iqsec3","Bves","Serpinb2","Ptgs2","Nfkbiz","Ell2","Pax1","Ankrd56","Fosl1","Dscl","Fap","Hcn1","Dkk2","Gch1","Pcsk1","Gpr3","Rgs2","Adcyap1","Coq10b","Crem")
   CNTLs_c5and6 <- c("Eif4h","Mllt11","Casp6","Atp2c1","Col1a1","Prkcb","Tubb3","Auh","Ncoa6","Ssna1","Actb","Rab7","Hs2st1","Gapdh","Map2k4","Mrpl50","Glra2","AI597468","Nnmt","Rala")
   primerSet_c5and6 <- c(IEGs_c5and6,LRGs_c5and6,CNTLs_c5and6)  
   

#### Newly defined gene classes

  iegs <- c("Fosb","Fos","Egr4","Npas4","Junb","Nr4a1","Maff","Egr1","Arc","Egr2","Klf4","Gadd45b","Egr3","Nr4a2","Ier2","Nr4a3",'Csrnp1',"Cyr61","Dusp1","Ccno","Fbxo33",'Adcyap1','Dusp5','Sik1','Fap',"Trib1","Amigo3","Gltscr2","Nfkbiz")   
  dprgs <- c('Rasd1','Fosl2','Pcsk1','Pax1','Vgf','Crh','Bdnf','Pdlim3')
  lrgs <- c('Cartpt','Serpinb2','Ptgs2','Acan','Gprc5a','Sertad1','Pthlh','Npffr2','Hspb3','Rel','Crhbp','Tnn','Bves','Gch1','Popdc3','Ptger4','Gpr3','Rgs2','Nppc','Pnoc','Ell2','Dkk2','Wif1','Mlf1','Pam','Penk','Iqsec3') 
   cntls <- c("Lhx5",'Rbpms','Crem','Hcn1','Atf7ip2','Dscl','Ankrd56',"Vamp4","Coq10b","Eif4h","Mllt11","Casp6","Atp2c1","Col1a1","Prkcb","Tubb3","Auh","Ncoa6","Ssna1","Actb","Rab7","Hs2st1","Gapdh","Map2k4","Mrpl50","Glra2","AI597468","Nnmt","Rala")
 
 all_primers <- c(iegs,dprgs,lrgs,cntls)
 
 
   iejs <- c("Fosb","Fos","Egr4","Npas4","Junb","Nr4a1","Maff","Egr1","Arc","Egr2","Klf4","Gadd45b","Egr3","Nr4a2","Ier2","Nr4a3")   

 
###Read in all chips that have been run until this point   
###if this file does not contain a modified_lucas.txt tag then it wont work using the given functions file that should already be in memory


   c1 <- read_results(fN="2011_01_11_Athar_and_5min_pulses_b1-b3/chip1_modified_lucas.txt",IEGs=IEGs_c1,LRGs=LRGs_c1,CNTLs=CNTLs_c1)

   c2 <- read_results(fN="2011_03_07/chip2_modified_lucas.txt",IEGs=IEGs_chips2to4,LRGs=LRGs_chips2to4,CNTLs=CNTLs_chips2to4)
   
   c4 <- read_results(fN="2011_04_11_sustained_KCl_CHX_ANS/chip4_modified_lucas.txt",IEGs=IEGs_chips2to4,LRGs=LRGs_chips2to4,CNTLs=CNTLs_chips2to4)
       
   c4_chx <- read_results(fN='2011_04_11_sustained_KCl_CHX_ANS/chx_edit.txt',IEGs=IEGs_chips2to4,LRGs=LRGs_chips2to4,CNTLs=CNTLs_chips2to4)

   doub_pulse <- read_results(fN='double pulse/2_pulse_edit_lucas.txt',IEGs=IEGs_chips2to4,LRGs=LRGs_chips2to4,CNTLs=CNTLs_chips2to4)
        
   c5 <- read_results(fN="2012_05_17_30minKCl_2hrKCl_WY1_and_10_30_55mM_5minPulses/missing_01_041812.txt",IEGs=IEGs_c5and6,LRGs=LRGs_c5and6,CNTLs=CNTLs_c5and6)
  
   c6 <- read_results(fN="2012_07_11/chip6_modified_lucas.txt",IEGs=IEGs_c5and6,LRGs=LRGs_c5and6,CNTLs=CNTLs_c5and6)

   


   

