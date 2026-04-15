library(ggplot2)
library(tidyverse)
library(broom)
library(lme4)
library(emmeans)
library(dplyr)

library(stringr)
library(readxl)
library(ggpubr)

setwd("/home/taslima/data/WittkoppLab/Research_Projects/PlasticNoise/")

# REF Strain descriptions: 
# Y1617 : EXP WT
# Y1002 : EXP REF
# Y2675 : EXP URA3 REF.2X
# Y2676 : EXP URA3 REF.2X-DUP
# Y2679 : EXP URA3 REF.2X-DUP
# Y978  : EXP NEGATIVE
# Y1189 : Fitness REF single
# Y2682: Fitness REF for double
# Y1177: Fitness deletion

#### Reference values
datExpGluRef<-read.delim("Data/GLUCOSE.WHEEL.DATA.txt") %>% 
  dplyr::filter(STRAIN %in% c("Y1617", "Y978"))  
datExpGalRef<-read.delim("Data/GALACTOSE.ROTOR.DATA.txt") %>% 
  dplyr::filter(STRAIN %in% c("Y1617", "Y978"))  
datExpGlyRef<-read.delim("Data/GLYCEROL.WHEEL.DATA.txt") %>% 
  dplyr::filter(STRAIN %in% c("Y1617", "Y978"))  
datExpEthRef<-read.delim("Data/ETHANOL.WHEEL.DATA.txt") %>% 
  dplyr::filter(STRAIN %in% c("Y1617", "Y978"))  

datREF<-rbind(datExpGluRef, datExpGalRef, 
              datExpGlyRef, datExpEthRef) %>% 
  mutate(ENVIRONMENT = factor(ENVIRONMENT, 
                              levels=c("GLUCOSE", "GALACTOSE","GLYCEROL","ETHANOL"))) %>% 
  select(STRAIN,MUTATION,ENVIRONMENT,YFP.MEDIAN.FINAL,YFP.SD.FINAL,YFP.SD.SCALED) %>% 
  group_by(STRAIN,MUTATION,ENVIRONMENT) %>% 
  summarise_all(funs(mean,sd, se =sd(.)/sqrt(n())))


######### Result Section 1 #############
###### A. Test the effect of TDH3 mutant alleles on expression noise 
# This one requires RAW data
# will use first YFP.SD.FINAL which is corrected for flow rotation, random positional variation
# After testing for GxE since there is only a small% variance explained by GXE we 
# will use YFP.SD.SCALED which is corrected for flow rotation, random positional variation
# and scaled to noise of WT

#Glucose
datExpGlu<-read.delim("Data/GLUCOSE.WHEEL.DATA.txt") %>% 
  dplyr::filter(!STRAIN %in% c("Y1002", "Y2675", "Y2676", "Y978"))  


mod<-anova( lm(YFP.SD.FINAL ~ STRAIN, data = datExpGlu))

cat("In Glucose", "df=", mod$Df[1], 
    "F-value=", round(mod$`F value`[1],1), 
    "p-value=", mod$`Pr(>F)`[1])

mod<-anova( lm(YFP.SD.FINAL/YFP.MEDIAN.FINAL ~ STRAIN, data = datExpGlu))
mod<-anova( lm(YFP.SD.FINAL^2/YFP.MEDIAN.FINAL ~ STRAIN, data = datExpGlu))

#Galactose
datExpGal<-read.delim("Data/GALACTOSE.ROTOR.DATA.txt") %>% 
  dplyr::filter(!STRAIN %in% c("Y1002", "Y2675", "Y2676", "Y978"))

mod<-anova( lm(YFP.SD.FINAL ~ STRAIN, data = datExpGal))

cat("In Galactose", "df=", mod$Df[1], 
    "F-value=", round(mod$`F value`[1],1), 
    "p-value=", mod$`Pr(>F)`[1])

mod<-anova( lm(YFP.SD.FINAL/YFP.MEDIAN.FINAL ~ STRAIN, data = datExpGal))
mod<-anova( lm(YFP.SD.FINAL^2/YFP.MEDIAN.FINAL ~ STRAIN, data = datExpGal))

#Glycerol
datExpGly<-read.delim("Data/GLYCEROL.WHEEL.DATA.txt") %>% 
  dplyr::filter(!STRAIN %in% c("Y1002", "Y2675", "Y2676", "Y978"))

mod<-anova( lm(YFP.SD.FINAL ~ STRAIN, data = datExpGly))

cat("In Glycerol", "df=", mod$Df[1], 
    "F-value=", round(mod$`F value`[1],1), 
    "p-value=", mod$`Pr(>F)`[1])

mod<-anova( lm(YFP.SD.FINAL/YFP.MEDIAN.FINAL ~ STRAIN, data = datExpGly))
mod<-anova( lm(YFP.SD.FINAL^2/YFP.MEDIAN.FINAL ~ STRAIN, data = datExpGly))

#Ethanol
datExpEth<-read.delim("Data/ETHANOL.WHEEL.DATA.txt") %>% 
  dplyr::filter(!STRAIN %in% c("Y1002", "Y2675", "Y2676", "Y978"))

mod<-anova( lm(YFP.SD.FINAL ~ STRAIN, data = datExpEth))

cat("In Ethanol", "df=", mod$Df[1], 
    "F-value=", round(mod$`F value`[1],1), 
    "p-value=", mod$`Pr(>F)`[1])

mod<-anova( lm(YFP.SD.FINAL/YFP.MEDIAN.FINAL ~ STRAIN, data = datExpEth))
mod<-anova( lm(YFP.SD.FINAL^2/YFP.MEDIAN.FINAL ~ STRAIN, data = datExpEth))

# pairwise correlation: SD
cor.test(datExpGlu$YFP.SD.FINAL, datExpGal$YFP.SD.FINAL)
cor.test(datExpGlu$YFP.SD.FINAL, datExpGly$YFP.SD.FINAL)
cor.test(datExpGlu$YFP.SD.FINAL, datExpEth$YFP.SD.FINAL)
cor.test(datExpGal$YFP.SD.FINAL, datExpGly$YFP.SD.FINAL)
cor.test(datExpGal$YFP.SD.FINAL, datExpEth$YFP.SD.FINAL)
cor.test(datExpGly$YFP.SD.FINAL, datExpEth$YFP.SD.FINAL)

# pairwise correlation: CV
cor.test(datExpGlu$YFP.SD.FINAL/datExpGlu$YFP.MEDIAN.FINAL, 
         datExpGal$YFP.SD.FINAL/datExpGal$YFP.MEDIAN.FINAL)
cor.test(datExpGlu$YFP.SD.FINAL/datExpGlu$YFP.MEDIAN.FINAL, 
         datExpGly$YFP.SD.FINAL/datExpGly$YFP.MEDIAN.FINAL)
cor.test(datExpGlu$YFP.SD.FINAL/datExpGlu$YFP.MEDIAN.FINAL, 
         datExpEth$YFP.SD.FINAL/datExpEth$YFP.MEDIAN.FINAL)
cor.test(datExpGal$YFP.SD.FINAL/datExpGal$YFP.MEDIAN.FINAL, 
         datExpGly$YFP.SD.FINAL/datExpGly$YFP.MEDIAN.FINAL)
cor.test(datExpGal$YFP.SD.FINAL/datExpGal$YFP.MEDIAN.FINAL, 
         datExpEth$YFP.SD.FINAL/datExpEth$YFP.MEDIAN.FINAL)
cor.test(datExpGly$YFP.SD.FINAL/datExpGly$YFP.MEDIAN.FINAL, 
         datExpEth$YFP.SD.FINAL/datExpEth$YFP.MEDIAN.FINAL)

# pairwise correlation: FANO
cor.test(datExpGlu$YFP.SD.FINAL^2/datExpGlu$YFP.MEDIAN.FINAL, 
         datExpGal$YFP.SD.FINAL^2/datExpGal$YFP.MEDIAN.FINAL)
cor.test(datExpGlu$YFP.SD.FINAL^2/datExpGlu$YFP.MEDIAN.FINAL, 
         datExpGly$YFP.SD.FINAL^2/datExpGly$YFP.MEDIAN.FINAL)
cor.test(datExpGlu$YFP.SD.FINAL^2/datExpGlu$YFP.MEDIAN.FINAL, 
         datExpEth$YFP.SD.FINAL^2/datExpEth$YFP.MEDIAN.FINAL)
cor.test(datExpGal$YFP.SD.FINAL^2/datExpGal$YFP.MEDIAN.FINAL, 
         datExpGly$YFP.SD.FINAL^2/datExpGly$YFP.MEDIAN.FINAL)
cor.test(datExpGal$YFP.SD.FINAL^2/datExpGal$YFP.MEDIAN.FINAL, 
         datExpEth$YFP.SD.FINAL^2/datExpEth$YFP.MEDIAN.FINAL)
cor.test(datExpGly$YFP.SD.FINAL^2/datExpGly$YFP.MEDIAN.FINAL, 
         datExpEth$YFP.SD.FINAL^2/datExpEth$YFP.MEDIAN.FINAL)

##### Get the regression coeff
mod1<-lm(datExpGlu$YFP.SD.FINAL~ datExpGal$YFP.SD.FINAL)
mod.emt1 <- emtrends(mod1, ~1, var="datExpGal$YFP.SD.FINAL")
test(mod.emt1, null=1)

mod2<-lm(datExpGlu$YFP.SD.FINAL ~ datExpGly$YFP.SD.FINAL)
mod.emt2 <- emtrends(mod2, ~1, var="datExpGly$YFP.SD.FINAL")
test(mod.emt2, null=1)

mod3<-lm(datExpGlu$YFP.SD.FINAL ~ datExpEth$YFP.SD.FINAL)
mod.emt3 <- emtrends(mod3, ~1, var="datExpEth$YFP.SD.FINAL")
test(mod.emt3, null=1)

mod4<-lm(datExpGal$YFP.SD.FINAL ~ datExpGly$YFP.SD.FINAL)
mod.emt4 <- emtrends(mod4, ~1, var="datExpGly$YFP.SD.FINAL")
test(mod.emt4, null=1)

mod5<-lm(datExpGal$YFP.SD.FINAL ~ datExpEth$YFP.SD.FINAL)
mod.emt5 <- emtrends(mod5, ~1, var="datExpEth$YFP.SD.FINAL")
test(mod.emt5, null=1)

mod6<-lm(datExpGly$YFP.SD.FINAL ~ datExpEth$YFP.SD.FINAL)
mod.emt6 <- emtrends(mod6, ~1, var="datExpEth$YFP.SD.FINAL")
test(mod.emt6, null=1)

#### Supplementary Figure 1: plot pairwise correlation of noise 
defaultmar<-c(5.1, 4.1 ,4.1, 2.1)
newmar<-c(5.1, 5.1 ,3.1, 2.1)

#tiff("Plots/SupFig_S1.tiff",width=8,height=6,units="in",res=300)
tiff("Plots/Fig1_part2.tiff",width=10,height=6,units="in",res=300)
par(mfrow=c(2,3), mar=newmar)
plot(datExpGlu$YFP.SD.FINAL, datExpGal$YFP.SD.FINAL,col="firebrick",
     xlab="Glucose",ylab="Galactose",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.6), ylim = c(0,0.6))+
  abline(lm(datExpGal$YFP.SD.FINAL ~ datExpGlu$YFP.SD.FINAL), col = "black", lty = 2)+
  #text(paste("r = ", round(cor(datExpGlu$YFP.SD.FINAL, datExpGal$YFP.SD.FINAL), 2)), 
  #     x = 0.05, y = 0.4, col = "navy", cex=1.5)
  text(paste("slope = ", round(as.numeric(lm(datExpGlu$YFP.SD.FINAL ~datExpGal$YFP.SD.FINAL)$coefficients[2]), 2)), 
       x = 0.35, y = 0.15, col = "navy", cex=1.5)
plot(datExpGlu$YFP.SD.FINAL, datExpGly$YFP.SD.FINAL,col="firebrick",
     xlab="Glucose",ylab="Glycerol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.6), ylim = c(0,0.6))+
  abline(lm(datExpGly$YFP.SD.FINAL ~ datExpGlu$YFP.SD.FINAL), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGlu$YFP.SD.FINAL~ datExpGly$YFP.SD.FINAL)$coefficients[2]), 2)), 
       x = 0.35, y = 0.15, col = "navy", cex=1.5)
plot(datExpGlu$YFP.SD.FINAL, datExpEth$YFP.SD.FINAL,col="firebrick",
     xlab="Glucose",ylab="Ethanol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.6), ylim = c(0,0.6))+
  abline(lm(datExpEth$YFP.SD.FINAL ~ datExpGlu$YFP.SD.FINAL), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGlu$YFP.SD.FINAL ~ datExpEth$YFP.SD.FINAL)$coefficients[2]), 2)), 
       x = 0.35, y = 0.15, col = "navy", cex=1.5)
plot(datExpGal$YFP.SD.FINAL, datExpGly$YFP.SD.FINAL,col="firebrick",
     xlab="Galactose",ylab="Glycerol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.6), ylim = c(0,0.6))+
  abline(lm(datExpGly$YFP.SD.FINAL ~ datExpGal$YFP.SD.FINAL), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGal$YFP.SD.FINAL ~ datExpGly$YFP.SD.FINAL)$coefficients[2]), 2)), 
       x = 0.35, y = 0.15, col = "navy", cex=1.5)
plot(datExpGal$YFP.SD.FINAL, datExpEth$YFP.SD.FINAL,col="firebrick",
     xlab="Galactose",ylab="Ethanol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.6), ylim = c(0,0.6))+
  abline(lm(datExpEth$YFP.SD.FINAL ~ datExpGal$YFP.SD.FINAL), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGal$YFP.SD.FINAL ~ datExpEth$YFP.SD.FINAL)$coefficients[2]), 2)), 
       x = 0.35, y = 0.15, col = "navy", cex=1.5)
plot(datExpGly$YFP.SD.FINAL, datExpEth$YFP.SD.FINAL,col="firebrick",
     xlab="Glycerol",ylab="Ethanol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.6), ylim = c(0,0.6))+
  abline(lm(datExpEth$YFP.SD.FINAL ~ datExpGly$YFP.SD.FINAL), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGly$YFP.SD.FINAL ~ datExpEth$YFP.SD.FINAL)$coefficients[2]), 2)), 
       x = 0.35, y = 0.15, col = "navy", cex=1.5)
dev.off()

## FANO
datExpGlu$NoiseF<-datExpGlu$YFP.SD.FINAL^2/datExpGlu$YFP.MEDIAN.FINAL
datExpGal$NoiseF<-datExpGal$YFP.SD.FINAL^2/datExpGal$YFP.MEDIAN.FINAL
datExpGly$NoiseF<-datExpGly$YFP.SD.FINAL^2/datExpGly$YFP.MEDIAN.FINAL
datExpEth$NoiseF<-datExpEth$YFP.SD.FINAL^2/datExpEth$YFP.MEDIAN.FINAL

## CV
datExpGlu$NoiseCV<-datExpGlu$YFP.SD.FINAL/datExpGlu$YFP.MEDIAN.FINAL
datExpGal$NoiseCV<-datExpGal$YFP.SD.FINAL/datExpGal$YFP.MEDIAN.FINAL
datExpGly$NoiseCV<-datExpGly$YFP.SD.FINAL/datExpGly$YFP.MEDIAN.FINAL
datExpEth$NoiseCV<-datExpEth$YFP.SD.FINAL/datExpEth$YFP.MEDIAN.FINAL

#### FANO
#tiff("Plots/SupFig_S1.tiff",width=10,height=6,units="in",res=300)
png("Plots/SupFig_S1.png",width=10,height=6,units="in",res=300)
par(mfrow=c(2,3), mar=newmar)
plot(datExpGlu$NoiseF, datExpGal$NoiseF,col="firebrick",
     xlab="Glucose",ylab="Galactose",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.2), ylim = c(0,0.2))+
  abline(lm(datExpGal$NoiseF ~ datExpGlu$NoiseF), col = "black", lty = 2)+
  #text(paste("r = ", round(cor(datExpGlu$YFP.SD.FINAL, datExpGal$YFP.SD.FINAL), 2)), 
  #     x = 0.05, y = 0.4, col = "navy", cex=1.5)
  text(paste("slope = ", round(as.numeric(lm(datExpGlu$NoiseF ~datExpGal$NoiseF)$coefficients[2]), 2)), 
       x = 0.15, y = 0.05, col = "navy", cex=1.5)
plot(datExpGlu$NoiseF, datExpGly$NoiseF,col="firebrick",
     xlab="Glucose",ylab="Glycerol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.2), ylim = c(0,0.2))+
  abline(lm(datExpGly$NoiseF ~ datExpGlu$NoiseF), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGlu$NoiseF~ datExpGly$NoiseF)$coefficients[2]), 2)), 
       x = 0.15, y = 0.05, col = "navy", cex=1.5)
plot(datExpGlu$NoiseF, datExpEth$NoiseF,col="firebrick",
     xlab="Glucose",ylab="Ethanol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.2), ylim = c(0,0.2))+
  abline(lm(datExpEth$NoiseF ~ datExpGlu$NoiseF), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGlu$NoiseF ~ datExpEth$NoiseF)$coefficients[2]), 2)), 
       x = 0.15, y = 0.05, col = "navy", cex=1.5)
plot(datExpGal$NoiseF, datExpGly$NoiseF,col="firebrick",
     xlab="Galactose",ylab="Glycerol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.2), ylim = c(0,0.2))+
  abline(lm(datExpGly$NoiseF ~ datExpGal$NoiseF), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGal$NoiseF ~ datExpGly$NoiseF)$coefficients[2]), 2)), 
       x = 0.15, y = 0.05, col = "navy", cex=1.5)
plot(datExpGal$NoiseF, datExpEth$NoiseF,col="firebrick",
     xlab="Galactose",ylab="Ethanol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.2), ylim = c(0,0.2))+
  abline(lm(datExpEth$NoiseF ~ datExpGal$NoiseF), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGal$NoiseF ~ datExpEth$NoiseF)$coefficients[2]), 2)), 
       x = 0.15, y = 0.05, col = "navy", cex=1.5)
plot(datExpGly$NoiseF, datExpEth$NoiseF,col="firebrick",
     xlab="Glycerol",ylab="Ethanol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.2), ylim = c(0,0.2))+
  abline(lm(datExpEth$NoiseF ~ datExpGly$NoiseF), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGly$NoiseF ~ datExpEth$NoiseF)$coefficients[2]), 2)), 
       x = 0.15, y = 0.05, col = "navy", cex=1.5)
dev.off()


### CV
#tiff("Plots/SupFig_S2.tiff",width=10,height=6,units="in",res=300)
png("Plots/SupFig_S2.png",width=10,height=6,units="in",res=300)
par(mfrow=c(2,3), mar=newmar)
plot(datExpGlu$NoiseCV, datExpGal$NoiseCV,col="firebrick",
     xlab="Glucose",ylab="Galactose",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.3), ylim = c(0,0.3))+
  abline(lm(datExpGal$NoiseCV ~ datExpGlu$NoiseCV), col = "black", lty = 2)+
  #text(paste("r = ", round(cor(datExpGlu$YFP.SD.FINAL, datExpGal$YFP.SD.FINAL), 2)), 
  #     x = 0.05, y = 0.4, col = "navy", cex=1.5)
  text(paste("slope = ", round(as.numeric(lm(datExpGlu$NoiseCV ~datExpGal$NoiseCV)$coefficients[2]), 2)), 
       x = 0.18, y = 0.05, col = "navy", cex=1.5)
plot(datExpGlu$NoiseCV, datExpGly$NoiseCV,col="firebrick",
     xlab="Glucose",ylab="Glycerol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.3), ylim = c(0,0.3))+
  abline(lm(datExpGly$NoiseCV ~ datExpGlu$NoiseCV), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGlu$NoiseCV~ datExpGly$NoiseCV)$coefficients[2]), 2)), 
       x = 0.18, y = 0.05, col = "navy", cex=1.5)
plot(datExpGlu$NoiseCV, datExpEth$NoiseCV,col="firebrick",
     xlab="Glucose",ylab="Ethanol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.3), ylim = c(0,0.3))+
  abline(lm(datExpEth$NoiseCV ~ datExpGlu$NoiseCV), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGlu$NoiseCV ~ datExpEth$NoiseCV)$coefficients[2]), 2)), 
       x = 0.18, y = 0.05, col = "navy", cex=1.5)
plot(datExpGal$NoiseCV, datExpGly$NoiseCV,col="firebrick",
     xlab="Galactose",ylab="Glycerol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.3), ylim = c(0,0.3))+
  abline(lm(datExpGly$NoiseCV ~ datExpGal$NoiseCV), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGal$NoiseCV ~ datExpGly$NoiseCV)$coefficients[2]), 2)), 
       x = 0.18, y = 0.05, col = "navy", cex=1.5)
plot(datExpGal$NoiseCV, datExpEth$NoiseCV,col="firebrick",
     xlab="Galactose",ylab="Ethanol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.3), ylim = c(0,0.3))+
  abline(lm(datExpEth$NoiseCV ~ datExpGal$NoiseCV), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGal$NoiseCV ~ datExpEth$NoiseCV)$coefficients[2]), 2)), 
       x = 0.18, y = 0.05, col = "navy", cex=1.5)
plot(datExpGly$NoiseCV, datExpEth$NoiseCV,col="firebrick",
     xlab="Glycerol",ylab="Ethanol",lty=1,lwd=0.75, cex.axis=2,cex.lab=2,
     xlim = c(0,0.3), ylim = c(0,0.3))+
  abline(lm(datExpEth$NoiseCV ~ datExpGly$NoiseCV), col = "black", lty = 2)+
  text(paste("slope = ", round(as.numeric(lm(datExpGly$NoiseCV ~ datExpEth$NoiseCV)$coefficients[2]), 2)), 
       x = 0.18, y = 0.05, col = "navy", cex=1.5)
dev.off()

# back to default
par(mfrow=c(1,1), mar=defaultmar)
# Joint analysis

dat<-rbind(datExpGlu, datExpGal, 
           datExpGly, datExpEth) %>% 
  mutate(ENVIRONMENT = factor(ENVIRONMENT, 
                              levels=c("GLUCOSE", "GALACTOSE","GLYCEROL","ETHANOL")))

### test for interaction of strains and environments on expression level and noise
library(effectsize)
# GxE
#Expression
modEXP <- aov(YFP.MEDIAN.FINAL ~ STRAIN*ENVIRONMENT, data = dat)
outEXP<-anova(modEXP)
eta_squared(modEXP, partial = F)

#Noise
modNoise <- aov(YFP.SD.FINAL ~ STRAIN*ENVIRONMENT, data = dat)
outNoise<-anova(modNoise)
eta_squared(modNoise, partial = F)


t.test(datExpGlu$YFP.SD.FINAL-datExpGal$YFP.SD.FINAL, mu=0,alternative = "less")$p.value
t.test(datExpGlu$YFP.SD.FINAL-datExpGly$YFP.SD.FINAL, mu=0,alternative = "less")$p.value
t.test(datExpGlu$YFP.SD.FINAL-datExpEth$YFP.SD.FINAL, mu=0,alternative = "less")$p.value
t.test(datExpGal$YFP.SD.FINAL-datExpGly$YFP.SD.FINAL, mu=0,alternative = "less")$p.value
t.test(datExpGal$YFP.SD.FINAL-datExpEth$YFP.SD.FINAL, mu=0,alternative = "less")$p.value
t.test(datExpGlu$YFP.SD.FINAL-datExpGly$YFP.SD.FINAL, mu=0,alternative = "less")$p.value
#to test for two-way env pairs: emmeans
# Fit the linear model (including interaction between Genotype and Environment)
#YFP.MEDIAN.FINAL
model <- lm(YFP.SD.FINAL ~ STRAIN * ENVIRONMENT, data=dat)
#model <- lm(YFP.SD.FINAL/YFP.MEDIAN.FINAL ~ STRAIN * ENVIRONMENT, data=dat)
#model <- lm(YFP.SD.FINAL^2/YFP.MEDIAN.FINAL ~ STRAIN * ENVIRONMENT, data=dat)

# Use emmeans to get estimated marginal means
emm <- emmeans(model, ~ ENVIRONMENT | STRAIN)

# Perform pairwise comparisons between environments for each genotype
pairwise_results <- pairs(emm, adjust = "none")

summary(pairwise_results)

# Optionally, to make the results more readable, you can convert them into a data frame
pairwise_df <- as.data.frame(pairwise_results) 

head(pairwise_df)

#Select strains which at least vary in one comparisons

pairwise_df_sig <- pairwise_df %>% 
  as_tibble() %>% 
  #dplyr::filter(`p.value` < 0.05) #this is already bonferroni corrected for 6 ways
  dplyr::filter(`p.value` < 0.05/(48*6)) #this is already bonferroni corrected for 6 ways; adjust fr 48 strains

length(unique(pairwise_df_sig$STRAIN))

write.csv(as_tibble(pairwise_df), "Results/Supplementary_Table7.csv", row.names = F)

cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", 
                "#0072B2", "#D55E00", "#CC79A7",
                "#999999", "#F0E442")

##### Supplementary Figure S1
(Fig1_3<-ggplot(pairwise_df)+
    geom_histogram(aes(x=estimate, color=contrast, fill=contrast), 
                   binwidth = 0.01, alpha=0.7,
                   show.legend = F)+
    geom_vline(xintercept = 0, color="grey50", lty=2)+
    facet_wrap(~contrast)+theme_classic()+
    scale_x_continuous(limits = c(-0.35,0.35))+
    scale_color_manual(values = cbPalette)+
    scale_fill_manual(values = cbPalette)+
    labs(x="Noise plasticity", y="Frequency")+
    theme_classic(base_size = 16)+
    theme(text = element_text(family = "sans",size=18),
          strip.background = element_rect(colour = NA,
                                          fill = "white", size = 1)))

tiff("Plots/Fig1_3.tiff",width=9,height=6,units="in",res=300)
Fig1_3
dev.off()

png("Plots/Fig1_3.png",width=8,height=6,units="in",res=300)
Fig1_3
dev.off()