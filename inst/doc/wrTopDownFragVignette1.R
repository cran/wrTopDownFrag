## ---- include = FALSE---------------------------------------------------------
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
marks <- data.frame(name=c("Ubiquitin\n76aa","Glutamate dehydrogenase 1\n501aa"),length=c(76,501))
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
protP2 <- cbind(na=names(protP), se=protP, ma=wrProteo::convAASeq2mass(protP,seqName=TRUE))
pepT1 <- makeFragments(protTab=protP2, minFragSize=3, maxFragSize=9, internFra=TRUE)

## ----makeFragments2, echo=TRUE------------------------------------------------
head(pepT1)

dim(pepT1)

## The repartition between types of fragments :
table(pepT1[,"ty"])

## Types of ambiguities encountered
table(pepT1[,"ambig"])


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
# Now we'll add some random noise
set.seed(2020)
obsMass2 <- as.numeric(obsMass1)
obsMass2 <- obsMass2 + runif(length(obsMass2), min=-2e-4, max=2e-4)

## ----identifFixedModif1, echo=TRUE--------------------------------------------

identP1 <- identifFixedModif(prot=protP[1], expMass=obsMass2, minFragSize=3, 
  maxFragSize=7, modTy=list(basMod=c("b","y")))     # should find 10term +10inter
 
## This function returns a list   
str(identP1)  
## $masMatch identifies each 
  

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

