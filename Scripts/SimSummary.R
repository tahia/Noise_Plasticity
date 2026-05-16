library(ggplot2)
library(tidyverse)

setwd("/home/taslima/data/WittkoppLab/Research_Projects/PlasticNoise/")
source("Scripts/SimparseFunctions.R")

######## Figure 4
seed(123456)
#low noise
ggplot()+ 
  geom_density(aes(x=rnorm(200000,mean = 1,sd=0.05)), size=8, color="#A7DE70")+
  theme_classic()+
  labs(x="", y="")+ 
  scale_x_continuous(limits = c(0,2))+
  scale_y_continuous(limits = c(0,8))+
  theme(axis.text = element_blank(), 
        #axis.text.x = element_text(size=40, color="black", face="bold"), 
        axis.ticks.y = element_blank(),
        axis.line= element_line(size = 4, color="grey50"), 
        panel.background = element_rect(fill = "transparent"))

#high noise
ggplot()+ 
  geom_density(aes(x=rnorm(200000,mean = 1,sd=0.25)), size=8, color="#C76E00")+
  theme_classic()+
  labs(x="", y="")+ 
  scale_x_continuous(limits = c(0,2))+
  scale_y_continuous(limits = c(0,8))+
  theme(axis.text = element_blank(), 
        #axis.text.x = element_text(size=40, color="black", face="bold"), 
        axis.ticks.y = element_blank(), 
        axis.line= element_line(size = 4, color="grey50"), 
        panel.background = element_rect(fill = "transparent"))

SumFitness_NM<- process_sims_comp("Data/Simultions/Python/DEAPSimOut/Final/Normal_Comp_twoD_G_c112_Fitvar1_1_Fitvar2_0.4.csv",
                                  Relmean = 1, RelSD = 0.05) %>%
  #filter(Heritability==0) %>%
  mutate(Function="Gaussian")

SumFitness_LN<-process_sims_comp("Data/Simultions/Python/DEAPSimOut/Final/LogNormal_Comp_twoD_G_c90_Fitvar1_0.5_Fitvar2_0.7.csv",
                                 Relmean = 1, RelSD = 0.05) %>%
  #filter(Heritability==0) %>% 
  mutate(Function="Lognormal")

SumFitness_MN<- process_sims_comp("Data/Simultions/Python/DEAPSimOut/Final/MixedNormal_Comp_twoD_G_c90_Fitvar1_0.5_0.005_Fitvar2_0.25_0.5_w_0.6.csv",
                                  Relmean = 0.5, RelSD = 0.05) %>%
  #filter(Heritability==0) %>% 
  mutate(Function="Mixed-Gaussian")


SumFitness_MLN<- process_sims_comp("Data/Simultions/Python/DEAPSimOut/Final/MixedLognormal_Comp_twoD_G_c90_Fitvar1_0.5_0.07_Fitvar2_0.25_0.75_w_0.5.csv",
                                   Relmean = 1.5, RelSD = 0.05) %>%
  #filter(Heritability==0) %>% 
  mutate(Function="Mixed-Lognormal")


SumFitness <-rbind(SumFitness_NM, SumFitness_LN, 
                   SumFitness_MN, SumFitness_MLN) %>% 
  mutate(Function=factor(Function, 
                         levels=c("Gaussian", "Lognormal", "Mixed-Gaussian", "Mixed-Lognormal"))) %>% 
  mutate(Mean=100*Mean) %>% 
  mutate(ExNoise = as.numeric(as.character(Noise))) %>% 
  filter(ExNoise <= 1600)


##### Figure 4
ggplot(SumFitness_NM, aes(y=Fitness, x=Mean))+
  geom_line(aes(y=Fitness, x=Mean, group=Noise, color=Noise),size=1.0,
            show.legend = F)+
  geom_point(aes(y=Fitness, x=Mean, color=Noise), shape=21,size=2,
             show.legend = F)+
  geom_linerange(aes(color=Noise,
                     ymin=Fitness-1.96*Fitness_se, ymax=Fitness+1.96*Fitness_se),
                 show.legend = F)+
  geom_hline(yintercept = c(0.5,0.75), color="grey50",lty=2)+
  #scale_colour_gradient(low = "#4040a1", high = "firebrick")+
  scale_y_continuous(breaks = c(0.5,0.75,1))+
  scale_color_brewer(palette = "Paired")+
  # labs(x="Relative Expression (%)", 
  #      y="Relative Fitness", color="Relative Noise (%)")+
  labs(x="", 
       y="", color="")+
  facet_wrap(~Heritability,nrow = 1 )+ #scales = "free_y"
  theme_classic()+
  theme(axis.title = element_text(size = 17),
        legend.title = element_text(size=10,face="bold"),
        legend.text= element_text(size=14),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3,"cm"),
        #legend.position=c(0.025, 0.8),
        axis.text = element_text(size=16,face = "bold",color="black"),
        panel.spacing.x =unit(0.15, "lines") , 
        panel.spacing.y=unit(0.15,"lines"),
        strip.text = element_blank(),
        #strip.text = element_text(size=14,face = "bold"),
        strip.background = element_rect(
          color="transparent", fill="grey90"))


##### Figure 5
SumFitness_H0 <-SumFitness %>% 
  filter(Heritability==0)

#### Expression-fitness function
x<-seq(0.01,2,0.01)
cN <- 112 / max(dnorm(seq(0.01, 2, by = 0.01), mean = 1, sd = 0.4))
yN<-191 - cN * dnorm(x, mean = 1, sd = 0.4)

cL <- 90 / max(dlnorm(seq(0.01, 2, by = 0.01), mean = 0.5, sd = 0.7))
yL<-160 - cL * dlnorm(x, mean = 0.5, sd = 0.7)

# cMN <- 90 / max(0.85*(dnorm(seq(0.01, 2, by = 0.01), mean = 1, sd = 0.5))+
#                 0.15*(dnorm(seq(0.01, 2, by = 0.01), mean = 0.1, sd = 0.25)))
# 
# yMN<-160 - cMN * (0.85*(dnorm(seq(0.01, 2, by = 0.01), mean = 1, sd = 0.5))+
#   0.15*(dnorm(seq(0.01, 2, by = 0.01), mean = 0.1, sd = 0.25)))
#mixed normal
cMN <- 90 / max(0.6*(dnorm(seq(0.01, 2, by = 0.01), mean = 0.5, sd = 0.25))+
                  0.4*(dnorm(seq(0.01, 2, by = 0.01), mean = 0.005, sd = 0.5)))

yMN<-160 - cMN * (0.6*(dnorm(seq(0.01, 2, by = 0.01), mean = 0.5, sd = 0.25))+
                    0.4*(dnorm(seq(0.01, 2, by = 0.01), mean = 0.005, sd = 0.5)))

# mixed lognormal
cML <- 90 / max( (0.5*(dlnorm(seq(0.01, 2, by = 0.01), mean = 0.5, sd = 0.25))+
                    0.5*(dlnorm(seq(0.01, 2, by = 0.01), mean = 0.07, sd = 0.75)))/1)

yML<-160 - cMN * (0.5*(dlnorm(seq(0.01, 2, by = 0.01), mean = 0.5, sd = 0.25))+
                    0.5*(dlnorm(seq(0.01, 2, by = 0.01), mean = 0.07, sd = 0.75)))/1


Expression<-as.data.frame(cbind(Expression=rep(x, 4),
                                DT=c(yN, yL, yMN, yML),
                                Function=c(rep("Gaussian",200), rep("Lognormal",200),
                                           rep("Mixed-Gaussian",200), rep("Mixed-Lognormal",200)))) %>% 
  mutate(Expression=as.numeric(Expression)) %>% 
  mutate(DT=as.numeric(DT)) %>% 
  mutate(Function=factor(Function, 
                         levels=c("Gaussian", "Lognormal", "Mixed-Gaussian", "Mixed-Lognormal")))


###### Get inflection points
# 2. Fit a smooth spline to your data points
spline_fit <- smooth.spline(x, yML)

# 3. Predict the second derivative (deriv = 2)
second_derivative <- predict(spline_fit, x, deriv = 2)$y

# 4. Find the inflection points (where the second derivative crosses 0)
inflection_indices <- which(diff(sign(second_derivative)) != 0)
inflection_points  <- x[inflection_indices]


## Gaussian c(0.59 1.40), Lognormal c(0.37 1.65)
## Mixed Gaussian c(0.23 0.74) and Mixed Lognormal c(0.01 0.20 0.70 1.17 1.91)

labeldata<-as.data.frame(cbind(
  Function=c(rep("Gaussian",2), rep("Lognormal",2),
             rep("Mixed-Gaussian",2), rep("Mixed-Lognormal",5)),
  Expression=c(c(0.59 ,1.40),c(0.37, 1.65),c(0.23, 0.74),
               c(0.01, 0.20, 0.70, 1.17, 1.91)))
)  %>% 
  mutate(Expression=as.numeric(Expression)*100) %>% 
  mutate(Function=factor(Function, 
                         levels=c("Gaussian", "Lognormal", "Mixed-Gaussian", "Mixed-Lognormal")))


(Fig5_Func<-ggplot(Expression)+ 
    geom_line(aes(x=Expression*100, y=DT), size=1, color="black")+
    theme_classic()+labs(x="Relative Expression (%)", y="DT")+ 
    scale_x_continuous(limits = c(0,200))+
    scale_y_reverse()+
    facet_wrap(~Function,nrow=1, scales = "free_y" )+
    theme(axis.title = element_text(size = 16),
          legend.title = element_text(size=10,face="bold"),
          legend.background = element_rect(fill = "transparent"),
          legend.position=c(0.5, 0.9),
          axis.text = element_text(size=14,color="black"),
          legend.text= element_text(size=8),
          panel.spacing.x =unit(0.15, "lines") , 
          panel.spacing.y=unit(0.15,"lines"),
          panel.background = element_rect(fill = "transparent"),
          strip.text = element_text(size=14,face = "bold"),
          strip.background = element_rect(color="transparent", fill="grey90")
    )
  
)

(Fig5_Fit<-ggplot(SumFitness_H0, aes(y=Fitness, x=Mean))+
    geom_line(aes(y=Fitness, x=Mean, group=Noise, color=Noise),size=1)+
    geom_point(aes(y=Fitness, x=Mean, color=Noise), shape=21, size=2)+
    geom_linerange(aes(color=Noise,
                       ymin=Fitness-1.96*Fitness_se, ymax=Fitness+1.96*Fitness_se),show.legend = T)+
    geom_vline(data=labeldata, 
               aes(xintercept = Expression), color="grey50",lty=2)+
    #scale_colour_gradient(low = "#4040a1", high = "firebrick")+
    scale_color_brewer(palette = "Paired")+
    labs(x="Relative Expression (%)", 
         y="Relative Fitness", color="Relative Noise (%)")+
    facet_wrap(~Function,nrow=1, scales = "free_y" )+
    theme_classic()+
    theme(axis.title = element_text(size = 16),
          legend.title = element_text(size=10,face="bold"),
          legend.text= element_text(size=6),
          legend.background = element_rect(fill = "transparent"),
          legend.key.size = unit(0.3,"cm"),
          legend.position=c(0.7, 0.6),
          axis.text = element_text(size=14, color="black"),
          panel.spacing.x =unit(0.15, "lines") , 
          panel.spacing.y=unit(0.15,"lines"),
          strip.text = element_text(size=14,face = "bold"),
          strip.background = element_rect(
            color="transparent", fill="grey90"))
)

(Fig5<-ggdraw()+
    #draw_image("Plots/SimSchematic2.jpg",x=0,  y = 0.4, scale = 1)+
    draw_plot(Fig5_Fit,x = 0, y = 0,width = 1,height = 0.5)+ 
    draw_plot(Fig5_Func,x = 0, y = 0.5,width = 1,height = 0.5)+
    draw_label(x=0.05,y=0.975,label = "A)", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.29,y=0.975,label = "B)", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.54,y=0.975,label = "C)", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.78,y=0.975,label = "D)", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.05,y=0.48,label = "E)", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.29,y=0.48,label = "F)", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.54,y=0.48,label = "G)", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.78,y=0.48,label = "H)", color = "black", size = 16, fontface = "bold")
  
)


tiff("Plots/Fig5.tiff",width=12,height =8,units="in",res=300)
Fig5
dev.off()

png("Plots/Fig5.png",width=12,height = 8,units="in",res=300)
Fig5
dev.off()

#### Test low/high noise at different expression levels
#Gaussian
Fitness_NM<-process_sims_nosum("Data/Simultions/Python/DEAPSimOut/Final/Normal_Comp_twoD_G_c112_Fitvar1_1_Fitvar2_0.4.csv",
                               Relmean = 1, RelSD = 0.05) %>%
  mutate(Function="Gaussian")

#Close
t.test(Fitness_NM$Fitness[Fitness_NM$Mean==1 & Fitness_NM$Heritability==0 & Fitness_NM$SD==0.05],
       Fitness_NM$Fitness[Fitness_NM$Mean==1 & Fitness_NM$Heritability==0 & Fitness_NM$SD==0.8])

#far
t.test(Fitness_NM$Fitness[Fitness_NM$Mean==0 & Fitness_NM$Heritability==0 & Fitness_NM$SD==0.05],
       Fitness_NM$Fitness[Fitness_NM$Mean==0 & Fitness_NM$Heritability==0 & Fitness_NM$SD==0.8])

SumFitness_NM %>% 
  filter(Heritability==0) %>% 
  select(Mean,SD,Fitness) %>% 
  group_by(SD) %>%  
  summarise(max_val = max(Fitness, na.rm = TRUE),
            exp = Mean[which.max(Fitness)])
#Lognormal
Fitness_LN<-process_sims_nosum("Data/Simultions/Python/DEAPSimOut/Final/LogNormal_Comp_twoD_G_c90_Fitvar1_0.5_Fitvar2_0.7.csv",
                               Relmean = 1, RelSD = 0.05) %>%
  mutate(Function="Gaussian")

#Close
t.test(Fitness_LN$Fitness[Fitness_LN$Mean==1 & Fitness_LN$Heritability==0 & Fitness_LN$SD==0.05],
       Fitness_LN$Fitness[Fitness_LN$Mean==1 & Fitness_LN$Heritability==0 & Fitness_LN$SD==0.8])

#far
t.test(Fitness_LN$Fitness[Fitness_LN$Mean==0 & Fitness_LN$Heritability==0 & Fitness_LN$SD==0.05],
       Fitness_LN$Fitness[Fitness_LN$Mean==0 & Fitness_LN$Heritability==0 & Fitness_LN$SD==0.8])
t.test(Fitness_LN$Fitness[Fitness_LN$Mean==2 & Fitness_LN$Heritability==0 & Fitness_LN$SD==0.05],
       Fitness_LN$Fitness[Fitness_LN$Mean==2 & Fitness_LN$Heritability==0 & Fitness_LN$SD==0.8])

SumFitness_LN %>% 
  filter(Heritability==0) %>% 
  select(Mean,SD,Fitness) %>% 
  group_by(SD) %>%  
  summarise(max_val = max(Fitness, na.rm = TRUE),
            exp = Mean[which.max(Fitness)])

#Mixed-Normal
Fitness_MN<-process_sims_nosum("Data/Simultions/Python/DEAPSimOut/Final/MixedNormal_Comp_twoD_G_c90_Fitvar1_0.5_0.005_Fitvar2_0.25_0.5_w_0.6.csv",
                               Relmean = 1, RelSD = 0.05) %>%
  mutate(Function="Mixed-Gaussian")

#Close
t.test(Fitness_MN$Fitness[Fitness_MN$Mean==0.5 & Fitness_MN$Heritability==0 & Fitness_MN$SD==0.05],
       Fitness_MN$Fitness[Fitness_MN$Mean==0.5 & Fitness_MN$Heritability==0 & Fitness_MN$SD==0.8])

#far
t.test(Fitness_MN$Fitness[Fitness_MN$Mean==0 & Fitness_MN$Heritability==0 & Fitness_MN$SD==0.05],
       Fitness_MN$Fitness[Fitness_MN$Mean==0 & Fitness_MN$Heritability==0 & Fitness_MN$SD==0.8])
t.test(Fitness_MN$Fitness[Fitness_MN$Mean==2 & Fitness_MN$Heritability==0 & Fitness_MN$SD==0.05],
       Fitness_MN$Fitness[Fitness_MN$Mean==2 & Fitness_MN$Heritability==0 & Fitness_MN$SD==0.8])

#intermediate
t.test(Fitness_MN$Fitness[Fitness_MN$Mean==0 & Fitness_MN$Heritability==0 & Fitness_MN$SD==0.4],
       Fitness_MN$Fitness[Fitness_MN$Mean==0 & Fitness_MN$Heritability==0 & Fitness_MN$SD==0.8])
t.test(Fitness_MN$Fitness[Fitness_MN$Mean==0 & Fitness_MN$Heritability==0 & Fitness_MN$SD==0.4],
       Fitness_MN$Fitness[Fitness_MN$Mean==0 & Fitness_MN$Heritability==0 & Fitness_MN$SD==0.05])


#Mixed-lognormal
Fitness_MLN<-process_sims_nosum("Data/Simultions/Python/DEAPSimOut/Final/MixedLognormal_Comp_twoD_G_c90_Fitvar1_0.5_0.07_Fitvar2_0.25_0.75_w_0.5.csv",
                                Relmean = 1.5, RelSD = 0.05) %>%
  mutate(Function="Mixed-LogNormal")

#Close
t.test(Fitness_MLN$Fitness[Fitness_MLN$Mean==1.5 & Fitness_MLN$Heritability==0 & Fitness_MLN$SD==0.05],
       Fitness_MLN$Fitness[Fitness_MLN$Mean==1.5 & Fitness_MLN$Heritability==0 & Fitness_MLN$SD==0.8])

#far
t.test(Fitness_MLN$Fitness[Fitness_MLN$Mean==0 & Fitness_MLN$Heritability==0 & Fitness_MLN$SD==0.05],
       Fitness_MLN$Fitness[Fitness_MLN$Mean==0 & Fitness_MLN$Heritability==0 & Fitness_MLN$SD==0.8])
t.test(Fitness_MLN$Fitness[Fitness_MLN$Mean==2 & Fitness_MLN$Heritability==0 & Fitness_MLN$SD==0.05],
       Fitness_MLN$Fitness[Fitness_MLN$Mean==2 & Fitness_MLN$Heritability==0 & Fitness_MLN$SD==0.8])

#intermediate
t.test(Fitness_MLN$Fitness[Fitness_MLN$Mean==2 & Fitness_MLN$Heritability==0 & Fitness_MLN$SD==0.4],
       Fitness_MLN$Fitness[Fitness_MLN$Mean==2 & Fitness_MLN$Heritability==0 & Fitness_MLN$SD==0.8])
t.test(Fitness_MLN$Fitness[Fitness_MLN$Mean==2 & Fitness_MLN$Heritability==0 & Fitness_MLN$SD==0.4],
       Fitness_MLN$Fitness[Fitness_MLN$Mean==2 & Fitness_MLN$Heritability==0 & Fitness_MLN$SD==0.05])


(FigS5<-ggplot(SumFitness, aes(y=Fitness, x=Mean))+
    geom_line(aes(y=Fitness, x=Mean, group=Noise, color=Noise),size=1)+
    geom_point(aes(y=Fitness, x=Mean, color=Noise), shape=21, size=1)+
    geom_linerange(aes(color=Noise,
                       ymin=Fitness-1.96*Fitness_se, ymax=Fitness+1.96*Fitness_se),show.legend = T)+
    #geom_hline(yintercept = c(1,0.875,0.855), color="red3",lty=2)+
    #scale_colour_gradient(low = "#4040a1", high = "firebrick")+
    scale_color_brewer(palette = "Paired")+
    labs(x="Relative Expression (%)", 
         y="Relative Fitness", color="Relative Noise (%)")+
    facet_wrap(Function~Heritability, scales = "free_y" )+
    theme_classic()+
    theme(axis.title = element_text(size = 14),
          legend.title = element_text(size=10,face="bold"),
          legend.text= element_text(size=6),
          legend.background = element_rect(fill = "transparent"),
          legend.key.size = unit(0.3,"cm"),
          #legend.position=c(0.7, 0.6),
          axis.text = element_text(size=14),
          panel.spacing.x =unit(0.15, "lines") , 
          panel.spacing.y=unit(0.15,"lines"),
          strip.text = element_text(size=12,face = "bold"),
          strip.background = element_rect(
            color="transparent", fill="grey90"))
)

tiff("Plots/SupFig_S5.tiff",width=12,height =8,units="in",res=300)
FigS5
dev.off()

png("Plots/SupFig_S5.png",width=12,height = 8,units="in",res=300)
FigS5
dev.off()