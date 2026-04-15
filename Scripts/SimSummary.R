library(ggplot2)
library(tidyverse)

setwd("/home/taslima/data/WittkoppLab/Research_Projects/PlasticNoise/")

###### Figure 4
HSumFitness_NM<- process_sims("Data/Simultions/Python/DEAPSimOut/Normal_Constant_G_c112_Fitvar1_1_Fitvar2_0.4.csv",
                              Relmean = 1, RelSD = 0.05) %>% 
  mutate(Function="Gaussian") %>% 
  mutate(ExNoise = as.numeric(as.character(Noise))) %>% 
  filter(ExNoise <= 1600)

ggplot(HSumFitness_NM, aes(y=Fitness, x=Mean))+
  geom_line(aes(y=Fitness, x=Mean, group=Noise, color=Noise),size=1.0,
            show.legend = F)+
  geom_point(aes(y=Fitness, x=Mean, color=Noise), shape=21,size=2,
             show.legend = F)+
  geom_linerange(aes(color=Noise,
                     ymin=Fitness-1.96*Fitness_se, ymax=Fitness+1.96*Fitness_se),
                 show.legend = F)+
  #geom_hline(yintercept = c(1,0.875,0.855), color="red3",lty=2)+
  #scale_colour_gradient(low = "#4040a1", high = "firebrick")+
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

HFitness_NM<-process_sims_nosum("Data/Simultions/Python/DEAPSimOut/Normal_Constant_G_c112_Fitvar1_1_Fitvar2_0.4.csv") %>% 
  mutate(Function="Gaussian") %>% 
  mutate(Noise=as.numeric(as.character(Noise))) %>% 
  mutate(abs_distance=abs(Initial_ExpMean - 1))

t.test( HFitness_NM$Fitness[HFitness_NM$Initial_ExpMean==0 & HFitness_NM$Initial_ExpSD==0.05 & HFitness_NM$Heritability==0.00],
        HFitness_NM$Fitness[HFitness_NM$Initial_ExpMean==0 & HFitness_NM$Initial_ExpSD==1.4 & HFitness_NM$Heritability==0.00])

t.test( HFitness_NM$Fitness[HFitness_NM$Initial_ExpMean==2 & HFitness_NM$Initial_ExpSD==0.05 & HFitness_NM$Heritability==0.00],
        HFitness_NM$Fitness[HFitness_NM$Initial_ExpMean==2 & HFitness_NM$Initial_ExpSD==1.4 & HFitness_NM$Heritability==0.00])

t.test( HFitness_NM$Fitness[HFitness_NM$Initial_ExpMean==1 & HFitness_NM$Initial_ExpSD==0.05 & HFitness_NM$Heritability==0.00],
        HFitness_NM$Fitness[HFitness_NM$Initial_ExpMean==1 & HFitness_NM$Initial_ExpSD==1.4 & HFitness_NM$Heritability==0.00])


####### Figure 5
#Gaussian
SumFitness_NM<- process_sims("Data/Simultions/Python/DEAPSimOut/Normal_Constant_G_c112_Fitvar1_1_Fitvar2_0.4.csv",
                             Relmean = 1, RelSD = 0.05) %>%
  filter(Heritability==0) %>%
  mutate(Function="Gaussian")

#Lognormal
SumFitness_LN<-process_sims("Data/Simultions/Python/DEAPSimOut/LogNormal_Constant_G_c90_Fitvar1_0.5_Fitvar2_0.7.csv",
                            Relmean = 1, RelSD = 0.05) %>% 
  filter(Heritability==0) %>% 
  mutate(Function="Lognormal")

#Mixed Gaussian
SumFitness_MN<- process_sims("Data/Simultions/Python/DEAPSimOut/Normal_Constant_G_c90_Fitvar1_0.5_0.005_Fitvar2_0.25_0.5_w0.6.csv",
                             Relmean = 0.5, RelSD = 0.05) %>% 
  filter(Heritability==0) %>% 
  mutate(Function="Mixed-Gaussian")

#Mixed Lognormal
SumFitness_MLN<- process_sims("Data/Simultions/Python/DEAPSimOut/LogNormal_Constant_G_c90_Fitvar1_0.5_0.07_Fitvar2_0.25_0.75_w0.5.csv",
                              Relmean = 1.5, RelSD = 0.05) %>% 
  filter(Heritability==0) %>% 
  mutate(Function="Mixed-Lognormal")


SumFitness <-rbind(SumFitness_NM, SumFitness_LN, 
                   SumFitness_MN, SumFitness_MLN) %>% 
  mutate(Function=factor(Function, 
                         levels=c("Gaussian", "Lognormal", "Mixed-Gaussian", "Mixed-Lognormal"))) %>% 
  mutate(Mean=100*Mean) %>% 
  mutate(ExNoise = as.numeric(as.character(Noise))) %>% 
  filter(ExNoise <= 1600)

(Fig5_Fit<-ggplot(SumFitness, aes(y=Fitness, x=Mean))+
    geom_line(aes(y=Fitness, x=Mean, group=Noise, color=Noise),size=1)+
    geom_point(aes(y=Fitness, x=Mean, color=Noise), shape=21, size=2)+
    geom_linerange(aes(color=Noise,
                       ymin=Fitness-1.96*Fitness_se, ymax=Fitness+1.96*Fitness_se),show.legend = T)+
    #geom_hline(yintercept = c(1,0.875,0.855), color="red3",lty=2)+
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

#### Expression-fitness function
x<-seq(0.01,2,0.01)
#Gaussian
cN <- 112 / max(dnorm(seq(0.01, 2, by = 0.01), mean = 1, sd = 0.4))
yN<-191 - cN * dnorm(x, mean = 1, sd = 0.4)

#Lognormal
cL <- 90 / max(dlnorm(seq(0.01, 2, by = 0.01), mean = 0.5, sd = 0.7))
yL<-160 - cL * dlnorm(x, mean = 0.5, sd = 0.7)


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