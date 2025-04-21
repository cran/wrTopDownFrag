## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment = "#>")

## ----setup1, echo=TRUE, warnings=FALSE----------------------------------------
library(wrMisc)
library(wrProteo)
library(wrTopDownFrag)

## ----AAfragSettings1, echo=TRUE-----------------------------------------------
## common settings
str(AAfragSettings())

## ----massDeFormula1, echo=TRUE------------------------------------------------
# Standard way to obtain the (monoisotopic) mass of water (H20) or a phosphorylation (PO3)
# Note that the number a molecule appears must be written in front of the molecule (no number means one occurance)
massDeFormula(c("2HO", "P3O"))

# Undereith this runs (for H20):
2*.atomicMasses()["H","mono"] +.atomicMasses()["O","mono"]   # H2O


## ----convAASeq2mass1, echo=TRUE-----------------------------------------------
# Let's define two small amino-acid sequences 
protP <- c(pepK="KPEPTIDRPEP", pepC="CEPEPTRT", pepC2="PECEPTRT")

# The sequence converted to mass  
convAASeq2mass(protP)   


## ----nFragm1, out.width="150%", out.heigth="80%", echo=TRUE-------------------
marks <- data.frame(name=c("Ubiquitin\n76aa","Glutamate dehydrogenase 1\n501aa"), length=c(76,501))
layout(matrix(1:2,ncol=2))
plotNTheor(x=20:750, log="", mark=marks)
plotNTheor(x=20:750, log="xy", mark=marks)
mtext("log/log scale", cex=0.8, line=0.1)
  

## ----fragmentSeq1, echo=TRUE--------------------------------------------------
protP <- c(pepK="KPEPTIDRPEP", pepC="CEPEPTRT")
## Basic output 
fragmentSeq(protP[1], minSize=3, internFragments=TRUE, pref="pepK")

## ----makeFragments1, echo=TRUE------------------------------------------------
## Elaborate output
protP2 <- cbind(na=names(protP), se=protP, ma=wrProteo::convAASeq2mass(protP, seqName=TRUE))
pepT1 <- makeFragments(protTab=protP2, minFragSize=3, maxFragSize=9, internFra=TRUE)

## ----makeFragments2, echo=TRUE------------------------------------------------
head(pepT1)

dim(pepT1)

## The repartition between types of fragments :
table(pepT1[,"ty"])

## Types of ambiguities encountered
table(pepT1[,"ambig"])


## ----example1a, echo=TRUE-----------------------------------------------------

Alb.seq <- "HLVDEPQNLIK" 
(Alb.mass <- wrProteo::convAASeq2mass("HLVDEPQNLIK", massTy="mono")) 
## Now the mass for the MH+ ion
(Alb.MHp <- wrProteo::convAASeq2mass("HLVDEPQNLIK", massTy="mono") + wrProteo::.atomicMasses()["H","mono"] - wrProteo::.atomicMasses()["e","mono"]) 

## ----example1b, echo=TRUE-----------------------------------------------------
## delta prospector
Alb.MHp - 1305.7161

## ----example1c, echo=TRUE-----------------------------------------------------
AlbFrag <- makeFragments(Alb.seq, minFragSize=2, maxFragSize=17, internFra=FALSE)
head(AlbFrag)
AlbFrag[which(AlbFrag[,"ty"]=="Nter"), c("seq","mass","beg","ty")]

## ----example1d, echo=TRUE-----------------------------------------------------
Alb.b <- as.numeric(AlbFrag[which(AlbFrag[,"ty"]=="Nter"),c("mass")]) - 1*wrProteo::.atomicMasses()["O","mono"] - wrProteo::.atomicMasses()["H","mono"] - wrProteo::.atomicMasses()["e","mono"] 
Alb.y <- as.numeric(AlbFrag[which(AlbFrag[,"ty"]=="Cter"),c("mass")]) + wrProteo::.atomicMasses()["H","mono"] - 1*wrProteo::.atomicMasses()["e","mono"] 

Alb.a <- as.numeric(AlbFrag[which(AlbFrag[,"ty"]=="Nter"),c("mass")]) - wrProteo::.atomicMasses()["C","mono"] - 2*wrProteo::.atomicMasses()["O","mono"] - wrProteo::.atomicMasses()["H","mono"] - 1*wrProteo::.atomicMasses()["e","mono"] 
Alb.x <- as.numeric(AlbFrag[which(AlbFrag[,"ty"]=="Cter"),c("mass")]) + 1*wrProteo::.atomicMasses()["C","mono"] + 1*wrProteo::.atomicMasses()["O","mono"] - 1*wrProteo::.atomicMasses()["H","mono"] - 1*wrProteo::.atomicMasses()["e","mono"] 

AlbFrag.int <- makeFragments(Alb.seq, minFragSize=2, maxFragSize=17, internFra=TRUE)
dim(AlbFrag.int)

## ----example1e, echo=TRUE-----------------------------------------------------
## Masses from Prospector for a,x,b- and y-ions
Alb.prs.a <- c(223.1553, 322.2238,	437.2507,	566.2933,	663.3461,	791.4046,	905.4476,	1018.5316,	1131.6157)
Alb.prs.x <- rev(c(1194.6365,	1081.5524, 982.4840, 867.4571, 738.4145, 641.3617, 513.3031, 399.2602, 286.1761)) #, 173.0921))

Alb.prs.b <- c(251.1503, 350.2187, 465.2456, 594.2882, 691.3410, 819.3995, 933.4425, 1046.5265, 1159.6106)
Alb.prs.y <- rev(c(1168.6572, 1055.5732, 956.5047, 841.4778, 712.4352, 615.3824, 487.3239, 373.2809, 260.1969))

cbind(wr.b=Alb.b, d.b=Alb.b - Alb.prs.b, wr.y=Alb.y, d.y= Alb.y - Alb.prs.y)

## ----makeFragmentsB1, echo=TRUE-----------------------------------------------
prot1 <- "RTVAAPSVFIFPPSDEQLKSG"

## ----makeFragmentsB2, echo=TRUE-----------------------------------------------
(prot1.MHp <- wrProteo::convAASeq2mass(prot1, massTy="mono") + wrProteo::.atomicMasses()["H","mono"] -wrProteo::.atomicMasses()["e","mono"])
 
## Besides, lateron we'll also need the mass of water
(H2O.mass <- 2*wrProteo::.atomicMasses()["H","mono"]  + wrProteo::.atomicMasses()["O","mono"])

## ----makeFragmentsB3, echo=TRUE-----------------------------------------------
##      b			     	      y
##     ---	  1 	R	21	   ---
##  258.1561	2 	T	20	2090.0804
##  357.2245	3 	V	19	1989.0328
##  428.2616	4 	A	18	1889.9644  <== use b4 as example
##  499.2987	5 	A	17	1818.9272  <=    and corresp y17
##       ...  .    ...        ...   

## ----makeFragmentsB4, echo=TRUE-----------------------------------------------
(prot1.b4.MHp <- wrProteo::convAASeq2mass("RTVA", massTy="mono") - wrProteo::.atomicMasses()["H","mono"] - wrProteo::.atomicMasses()["O","mono"] - wrProteo::.atomicMasses()["e","mono"]) 

## ----makeFragmentsB5, echo=TRUE-----------------------------------------------
(prot1.y18.MHp <- wrProteo::convAASeq2mass("APSVFIFPPSDEQLKSG", massTy="mono") + wrProteo::.atomicMasses()["H","mono"] - wrProteo::.atomicMasses()["e","mono"])  

## ----makeFragmentsB6, echo=TRUE-----------------------------------------------
(prot1.MHp) - (prot1.b4.MHp + prot1.y18.MHp - wrProteo::.atomicMasses()["H","mono"] + wrProteo::.atomicMasses()["e","mono"])     

## ----makeFragmentsB7, echo=TRUE-----------------------------------------------
prot1Frag <- makeFragments(prot1, minFragSize=4, maxFragSize=17, internFra=FALSE)
prot1Frag[c(2,27),]            # our prevous b- & y-ion (as neutral peptide mass)

## ----makeFragmentsB8, echo=TRUE-----------------------------------------------
as.numeric(prot1Frag[2,"mass"]) - H2O.mass + wrProteo::.atomicMasses()["H","mono"] - wrProteo::.atomicMasses()["e","mono"]  

## ----makeFragmentsB9, echo=TRUE-----------------------------------------------
as.numeric(prot1Frag[27,"mass"]) + wrProteo::.atomicMasses()["H","mono"] - wrProteo::.atomicMasses()["e","mono"]  

## ----makeFragmentsC1, echo=TRUE-----------------------------------------------
## Prospector internal fragments for "PSDEQL":
## Internal Sequence      	b	        a	    b-NH3	    b-H2O
##             ...        ...       ...       ...       ...
##          IFPPSD   657.3243  629.3293       ---  639.3137
##          PSDEQL	 670.3042	 642.3093  653.2777  652.2937

## ----makeFragmentsC2, echo=TRUE-----------------------------------------------
prot1FragInt <- makeFragments(prot1, minFragSize=4, maxFragSize=17, internFra=TRUE)
dim(prot1FragInt)

## ----makeFragmentsC3, echo=TRUE-----------------------------------------------
PSDEQL.idx <- which(prot1FragInt[,2] =="PSDEQL")                                   
prot1FragInt[PSDEQL.idx -(1:0),]

## ----makeFragmentsC4, echo=TRUE-----------------------------------------------
(PSDEQL.MH <- as.numeric(prot1FragInt[PSDEQL.idx,"mass"]) - H2O.mass + wrProteo::.atomicMasses()["H","mono"])
(PSDEQL.MHp <- as.numeric(prot1FragInt[PSDEQL.idx,"mass"]) - H2O.mass + wrProteo::.atomicMasses()["H","mono"] - wrProteo::.atomicMasses()["e","mono"])

c(Prospector=670.3042) - c(PSDEQL.MH, PSDEQL.MHp)

## ----toyData1, echo=TRUE------------------------------------------------------
# The toy peptide/protein sequnce
protP <- c(pepK="KPEPTIDRPEP", pepC="CEPEPTRT")

obsMass1 <- cbind(a=c(424.2554,525.3031,638.3872,753.4141,909.5152,1006.5680,1135.6106),
  b=c(452.2504,553.2980,666.3821,781.4090,937.5102,1034.5629,1163.6055),
  x=c(524.2463,639.2733,752.3573,853.4050,950.4578,1079.5004,1176.5531),
  y=c(498.2671,613.2940,726.3781,827.4258,924.4785,1053.5211,1150.5739),
  bdH=c(434.2398,535.2875,648.3715,763.3985,919.4996,1016.5524,1145.5949),
  ydH=c(480.2565,1132.5633,595.2835,708.3675,809.4152,906.4680,1035.5106),
  bi=c(498.2307,583.3198,583.3198,611.3148,680.3726,712.3624,712.3624),
  bidH=c(662.3620,694.3519,694.3519,791.4046,791.4046,791.4046,888.4574),
  bidN=c(663.3461,695.3359,695.3359,792.3886,792.3886,792.3886,889.4414),
  ai=c(652.3777,684.3675,684.3675,781.4203,781.4203,781.4203,878.4730) )
rownames(obsMass1) <- c("P","T","I","D","R","P","E")      # only for N-term (a & b)

## This example contains several iso-mass cases
length(obsMass1)
length(unique(as.numeric(obsMass1)))


## ----toyData2, echo=TRUE------------------------------------------------------
## have same mass   
##       480.2201 	DRPE-H2O & y4-H2O
##       583.3198 	PTIDR & TIDRP ;         565.3093 	PTIDR-H2O & TIDRP-H2O 
##       809.4152 	y7-H2O  & EPTIDRP  
##       684.3675 	TIDRPE-CO & EPTIDR-CO ; 694.3519 	EPTIDR-H2O & TIDRPE-H2O
##       712.3624 	TIDRPE & EPTIDR   


## ----toyData3, echo=TRUE------------------------------------------------------
# Now we'll add some random noise (to mimick experimental data)
set.seed(2020)
obsMass2 <- as.numeric(obsMass1)
obsMass2 <- obsMass2 + runif(length(obsMass2), min=-2e-4, max=2e-4)

## ----identifFixedModif2, echo=TRUE--------------------------------------------
identP1 <- identifFixedModif(prot=protP[1], expMass=obsMass2, minFragSize=3, maxFragSize=7, modTy=list(basMod=c("b","y")))     # should find 10term +10inter
 
## This function returns a list   
str(identP1)  
## $masMatch identifies each 
  

## ----identifFixedModif2b, echo=TRUE-------------------------------------------
identP1$preMa[ which(identP1$preMa[,"no"] %in% (names(identP1$massMatch))),c("seq","orig","ty","beg","end","precAA","finMass","mod")] 

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

