library(ggplot2)
library(tidyverse)
library(broom)
library(lme4)
library(emmeans)
library(dplyr)

library(stringr)
library(readxl)
library(ggpubr)

library(cowplot)
library(ggh4x)

setwd("/home/taslima/data/WittkoppLab/Research_Projects/PlasticNoise/")

################## Joint analysis for expression-noise-fitness begins here #####################3
# Read expression data
EXP<-read.csv("Data/FabienOriginal.csv") %>% 
  as_tibble() %>% 
  mutate(MUTATION=if_else(STRAIN=="Y2671", "URA3DOUBLE", MUTATION)) %>% 
  dplyr::filter(!STRAIN=="Y2676") %>% #remove one extra URA3 strain
  mutate(ID=paste(MUTATION, ENVIRONMENT,sep = "_"))

# Raw Data from Fabien to get WT fitness
FITNESS_WT<- read.table("Data/ALL.DATA.FITNESS.txt",header=TRUE,as.is=TRUE) %>% 
  mutate(w.estimate.2=exp(s.estimate)) %>% 
  select(STRAIN,MUTATION, ENVIRONMENT,w.estimate.2) %>% 
  dplyr::filter(STRAIN=="Y1001") %>% 
  group_by(STRAIN,MUTATION, ENVIRONMENT) %>% 
  summarise(Fitness=median(w.estimate.2, na.rm = T))


FITNESS_WTDB<- read.table("Data/ALL.DATA.FITNESS.txt",header=TRUE,as.is=TRUE) %>% 
  mutate(w.estimate.2=exp(s.estimate)) %>% 
  select(STRAIN,MUTATION, ENVIRONMENT,w.estimate.2) %>% 
  dplyr::filter(STRAIN=="Y2682") %>% 
  group_by(STRAIN,MUTATION, ENVIRONMENT) %>% 
  summarise(Fitness=median(w.estimate.2, na.rm = T))

#This is RAW data from Fabien  
FITNESS<- read.table("Data/ALL.DATA.FITNESS.txt",header=TRUE,as.is=TRUE) %>% 
  mutate(YFP.CONSTRUCT=if_else(STRAIN=="Y2682", "DOUBLE", YFP.CONSTRUCT)) %>% 
  mutate(w.estimate.2= if_else(
    ENVIRONMENT=="GLUCOSE", exp(s.estimate)/FITNESS_WT$Fitness[FITNESS_WT$ENVIRONMENT=="GLUCOSE"],
    if_else(ENVIRONMENT=="GALACTOSE", exp(s.estimate)/FITNESS_WT$Fitness[FITNESS_WT$ENVIRONMENT=="GALACTOSE"], 
            if_else( ENVIRONMENT=="GLYCEROL", exp(s.estimate)/FITNESS_WT$Fitness[FITNESS_WT$ENVIRONMENT=="GLYCEROL"] ,
                     exp(s.estimate)/FITNESS_WT$Fitness[FITNESS_WT$ENVIRONMENT=="ETHANOL"]) ) )) %>% 
  mutate(w.estimate.2= if_else(YFP.CONSTRUCT=="SINGLE", w.estimate.2,
                               if_else(
                                 ENVIRONMENT=="GLUCOSE", exp(s.estimate)/FITNESS_WTDB$Fitness[FITNESS_WTDB$ENVIRONMENT=="GLUCOSE"],
                                 if_else(ENVIRONMENT=="GALACTOSE", exp(s.estimate)/FITNESS_WTDB$Fitness[FITNESS_WTDB$ENVIRONMENT=="GALACTOSE"],
                                         if_else( ENVIRONMENT=="GLYCEROL", exp(s.estimate)/FITNESS_WTDB$Fitness[FITNESS_WTDB$ENVIRONMENT=="GLYCEROL"] ,
                                                  exp(s.estimate)/FITNESS_WTDB$Fitness[FITNESS_WTDB$ENVIRONMENT=="ETHANOL"]) ) )
  ) ) %>%
  select(STRAIN,MUTATION, ENVIRONMENT,w.estimate.2) %>% 
  group_by(STRAIN,MUTATION, ENVIRONMENT) %>% 
  summarise(Fitness=median(w.estimate.2, na.rm = T),
            Fitness.mad=mad(w.estimate.2, na.rm = T),
            N=n(),
            Fitness.se=sd(w.estimate.2, na.rm = T)/sqrt(N)) %>% 
  mutate(Low.95=Fitness-1.96*Fitness.se,
         High.95=Fitness+.96*Fitness.se) %>% 
  select(STRAIN:Fitness,Low.95,High.95) %>% 
  mutate(ID=paste(MUTATION,ENVIRONMENT,sep = "_")) %>% 
  rename(STRAIN.Fitness=STRAIN,
         MUTATION.Fitness=MUTATION) %>% 
  relocate(ID,.before = STRAIN.Fitness )  %>% 
  mutate(MUTATION.Fitness=if_else(STRAIN.Fitness=="Y2679", "URA3DOUBLE", MUTATION.Fitness)) %>% 
  dplyr::filter(!STRAIN.Fitness %in% c("Y1189", "Y2679")) #1189: Barry;2679 2nd ura double

######## Main analysis begins
CombDat<-EXP %>% 
  # dplyr::select(MUTATION, YFP.MEDIAN.RELATIVE.MEAN,
  #                              YFP.FANO.RELATIVE.MEAN,ID) %>% 
  inner_join(.,FITNESS,
             by = c("ID" = "ID")) %>% 
  rename(ENVIRONMENT = ENVIRONMENT.x) %>% 
  mutate(Category=if_else(abs(1-Fitness) < 0.005, "Low-delta-Fi","High-delta-Fit")) %>% 
  drop_na(ENVIRONMENT) %>% 
  mutate(ENVIRONMENT=factor(ENVIRONMENT,
                            levels = c("GLUCOSE","GALACTOSE",
                                       "GLYCEROL","ETHANOL") ))


# ##################### Categorize low and high noise as Fabien classified
# ##-Fit non linear regression to Expression data separately for a given environment
# CombDat<-as_tibble(CombDat) %>% 
#   mutate(DeltaNoise = 0) 

ENVs<- levels(as.factor(CombDat$ENVIRONMENT))

for (i in ENVs) {
  dat <-as_tibble(CombDat) %>% 
    dplyr::filter(ENVIRONMENT==i)
  #dplyr::filter(YFP.MEDIAN.RELATIVE.MEAN <= 1.25) #Filter relative exp <=125%
  
  MODEL_FANO <- loess(dat$YFP.FANO.RELATIVE.MEAN ~dat$YFP.MEDIAN.RELATIVE.MEAN,degree=2,span=2/3)
  MODEL_SD <- loess(dat$YFP.SD.RELATIVE.MEAN ~dat$YFP.MEDIAN.RELATIVE.MEAN,degree=2,span=2/3)
  MODEL_CV <- loess(dat$YFP.CV.RELATIVE.MEAN ~dat$YFP.MEDIAN.RELATIVE.MEAN,degree=2,span=2/3)
  
  MODEL_Fit <- loess(dat$Fitness ~ dat$YFP.MEDIAN.RELATIVE.MEAN,degree=2,span=2/3)
  
  dat<- as_tibble(dat) %>% 
    mutate(DeltaFANO = YFP.FANO.RELATIVE.MEAN - predict(MODEL_FANO, dat$YFP.MEDIAN.RELATIVE.MEAN, se=TRUE)$fit ) %>% 
    mutate(DeltaSD = YFP.SD.RELATIVE.MEAN - predict(MODEL_SD, dat$YFP.MEDIAN.RELATIVE.MEAN, se=TRUE)$fit ) %>% 
    mutate(DeltaCV = YFP.CV.RELATIVE.MEAN - predict(MODEL_CV, dat$YFP.MEDIAN.RELATIVE.MEAN, se=TRUE)$fit ) %>% 
    mutate(DeltaFitness = Fitness - predict(MODEL_Fit, dat$YFP.MEDIAN.RELATIVE.MEAN, se=TRUE)$fit )
  
  assign(paste("CombDat", i,sep = "_"),dat)
  
  
  # Expression
  x.mid <- seq(min(dat$YFP.MEDIAN.RELATIVE.MEAN),max(dat$YFP.MEDIAN.RELATIVE.MEAN),by=0.001)
  y.mid <- predict(MODEL_SD, x.mid, se=TRUE)
  y.mid.cv <- predict(MODEL_CV, x.mid, se=TRUE)
  y.mid.fano <- predict(MODEL_FANO, x.mid, se=TRUE)
  
  #Fitness
  x.mid.fit <- seq(min(dat$YFP.MEDIAN.RELATIVE.MEAN),max(dat$YFP.MEDIAN.RELATIVE.MEAN),by=0.001)
  y.mid.fit <- predict(MODEL_Fit, x.mid, se=TRUE)
  
  assign(paste("Model", i,sep = "_"),
         data.frame(cbind(EXPRESSION=x.mid, SD=y.mid$fit, 
                          CV=y.mid.cv$fit,FANO=y.mid.fano$fit,
                          ENVIRONMENT=rep(i,length(x.mid)) )))
  
  assign(paste("ModelFIT", i,sep = "_"),
         data.frame(cbind(EXPRESSION=x.mid.fit, FITNESS=y.mid.fit$fit, ENVIRONMENT=rep(i,length(x.mid.fit)) )))
  
  #Get correlation stats
  LM_Cor<-cor.test(dat$YFP.MEDIAN.RELATIVE.MEAN,dat$YFP.SD.RELATIVE.SD, method = "pearson")
  assign(paste("LMModel", i,sep = "_"),
         data.frame(cbind(Pearson=as.numeric(LM_Cor$estimate), 
                          Pval=as.numeric(LM_Cor$p.value), ENVIRONMENT=rep(i,1) )))
  
  LM_Cor_CV<-cor.test(dat$YFP.MEDIAN.RELATIVE.MEAN,dat$YFP.CV.RELATIVE.SD, method = "pearson")
  assign(paste("LMModel_CV", i,sep = "_"),
         data.frame(cbind(Pearson=as.numeric(LM_Cor_CV$estimate), 
                          Pval=as.numeric(LM_Cor_CV$p.value), ENVIRONMENT=rep(i,1) )))
  
  LM_Cor_FANO<-cor.test(dat$YFP.MEDIAN.RELATIVE.MEAN,dat$YFP.FANO.RELATIVE.SD, method = "pearson")
  assign(paste("LMModel_FANO", i,sep = "_"),
         data.frame(cbind(Pearson=as.numeric(LM_Cor_FANO$estimate), 
                          Pval=as.numeric(LM_Cor_FANO$p.value), ENVIRONMENT=rep(i,1) )))
}

CombDat_NLR<-rbind(CombDat_GLUCOSE,CombDat_GALACTOSE,CombDat_GLYCEROL, CombDat_ETHANOL) %>% 
  mutate(Class= case_when(STRAIN %in% c("Y1002","Y2675", "Y1617") ~ "REF",
                          STRAIN %in% c("Y1898","Y1960","Y1897","Y1896",
                                        "Y1892","Y1901","Y1903","Y2455",
                                        "Y2458","Y2465","Y2460","Y2469",
                                        "Y2475","Y2478") ~ "TFBS",
                          STRAIN %in% c("Y2550","Y2651", "Y2656","Y2546",
                                        "Y2528","Y2566","Y2605") ~ "TATA",
                          STRAIN %in% c("Y2653","Y2812","Y2814","Y2820",
                                        "Y2816","Y2824","Y2826","Y2828",
                                        "Y2832","Y2834","Y2836","Y2840",
                                        "Y2842") ~ "TATA+TFBS",
                          STRAIN %in% c("Y2671","Y3019","Y3021","Y3022",
                                        "Y3010", "Y3028","Y3034","Y3040",
                                        "Y3006", "Y3011","Y3012","Y3025",
                                        "Y3015") ~ "TWO Copies"))

Concat_Model_SD <-rbind(Model_GLUCOSE,Model_GALACTOSE,Model_GLYCEROL,Model_ETHANOL) %>% 
  mutate(EXPRESSION=as.numeric(as.character(EXPRESSION))) %>% 
  mutate(SD=as.numeric(as.character(SD))) %>% 
  mutate(CV=as.numeric(as.character(CV))) %>% 
  mutate(FANO=as.numeric(as.character(FANO)))

Concat_Model_Fitness <-rbind(ModelFIT_GLUCOSE,ModelFIT_GALACTOSE,
                             ModelFIT_GLYCEROL,ModelFIT_ETHANOL) %>% 
  mutate(EXPRESSION=as.numeric(as.character(EXPRESSION))) %>% 
  mutate(FITNESS=as.numeric(as.character(FITNESS)))


Concat_Regression_SD<-rbind(LMModel_GLUCOSE, LMModel_GALACTOSE,
                            LMModel_GLYCEROL, LMModel_ETHANOL) %>% 
  mutate(Pearson=round(as.numeric(as.character(Pearson)),2)) %>% 
  mutate(Pval=round(as.numeric(as.character(Pval)),6))

Concat_Regression_CV<-rbind(LMModel_CV_GLUCOSE, LMModel_CV_GALACTOSE,
                            LMModel_CV_GLYCEROL, LMModel_CV_ETHANOL) %>% 
  mutate(Pearson=round(as.numeric(as.character(Pearson)),2)) %>% 
  mutate(Pval=round(as.numeric(as.character(Pval)),6))

Concat_Regression_FANO<-rbind(LMModel_FANO_GLUCOSE, LMModel_FANO_GALACTOSE,
                              LMModel_FANO_GLYCEROL, LMModel_FANO_ETHANOL) %>% 
  mutate(Pearson=round(as.numeric(as.character(Pearson)),2)) %>% 
  mutate(Pval=round(as.numeric(as.character(Pval)),6))


#### Test fitness difference between low and high noise group at a given env

dat<-CombDat_NLR %>% 
  mutate(N=as.numeric(as.character(N))) %>% 
  mutate(Category=if_else(DeltaSD < 0.01, "Low Noise","High Noise"))

# Fit the linear model (including interaction between Genotype and Environment)
model <- lm(Fitness ~ Category * ENVIRONMENT, data=dat)

# Use emmeans to get estimated marginal means
emm <- emmeans(model, ~ Category | ENVIRONMENT)

# Perform pairwise comparisons between environments for each genotype
pairwise_results <- pairs(emm, adjust = "bonferroni")

summary(pairwise_results)

#### Figure 2
(Fig2_A <-as_tibble(CombDat_NLR) %>% 
    mutate(N=as.numeric(as.character(N))) %>% 
    mutate(Category=if_else(DeltaSD < 0.01, "Low Noise","High Noise")) %>% 
    mutate(YFP.MEDIAN.RELATIVE.MEAN=100*YFP.MEDIAN.RELATIVE.MEAN,
           YFP.SD.RELATIVE.MEAN=100*YFP.SD.RELATIVE.MEAN,
           YFP.MEDIAN.RELATIVE.SD=100*YFP.MEDIAN.RELATIVE.SD,
           YFP.SD.RELATIVE.SD=100*YFP.SD.RELATIVE.SD
    ) %>% 
    ggplot()+
    geom_errorbar(aes(x=YFP.MEDIAN.RELATIVE.MEAN, y=YFP.SD.RELATIVE.MEAN,
                      xmin= YFP.MEDIAN.RELATIVE.MEAN-1.96*YFP.MEDIAN.RELATIVE.SD/sqrt(N),
                      xmax= YFP.MEDIAN.RELATIVE.MEAN+1.96*YFP.MEDIAN.RELATIVE.SD/sqrt(N),
    ),color="gray70") +
    geom_errorbar(aes(x=YFP.MEDIAN.RELATIVE.MEAN, y=YFP.SD.RELATIVE.MEAN,
                      ymin= YFP.SD.RELATIVE.MEAN-1.96*YFP.SD.RELATIVE.SD/sqrt(N),
                      ymax= YFP.SD.RELATIVE.MEAN+1.96*YFP.SD.RELATIVE.SD/sqrt(N),
    ),color="gray70") +
    geom_point(aes(x=YFP.MEDIAN.RELATIVE.MEAN,y=YFP.SD.RELATIVE.MEAN, 
                   fill=Category),shape=21 ,size=2.5,color="black")+
    geom_hline(yintercept = 100,linetype=2, color="grey50")+
    geom_vline(xintercept = 100,linetype=2, color="grey50")+
    geom_line(dat=Concat_Model_SD,aes(x=EXPRESSION*100, y=SD*100))+
    # geom_text(data=Concat_Regression_SD,aes(x=20, y=200,
    #                                         label=paste(paste0("r = ",Pearson), 
    #                                                     paste0("P-value = ", Pval),sep="\n" ) )) +
    # facet_wrap(~factor(ENVIRONMENT, levels = c("GLUCOSE","GALACTOSE",
    #                                            "GLYCEROL","ETHANOL")))+
    facet_wrap2(~factor(ENVIRONMENT, levels = c("GLUCOSE","GALACTOSE",
                                                "GLYCEROL","ETHANOL")), ncol=4,
                strip = strip_themed(
                  background_x = elem_list_rect(
                    fill = c("#56B4E9","#009E73","#D55E00","#0072B2"),color= rep("transparent",4))
                ) 
    )+
    scale_color_manual(values = c("#C76E00", "#A7DE70"))+ 
    scale_fill_manual(values = c("#C76E00", "#A7DE70"))+
    #scale_color_manual(values = c("#E6AB02", "#66A61E"))+ 
    #scale_fill_manual(values = c("#E6AB02", "#66A61E"))+
    labs(x="Median Expression relative to WT (%)", y="Expression Noise relative to WT (%)",
         color="", fill="")+
    theme_classic()+
    theme(axis.title = element_text(size = 14),
          legend.title = element_text(size=12,face="bold"),
          legend.background = element_rect(fill = "transparent"),
          legend.position=c(0.06, 0.8),
          axis.text = element_text(size=14),
          legend.text= element_text(size=12),
          panel.spacing.x =unit(0.15, "lines") , 
          panel.spacing.y=unit(0.15,"lines"),
          strip.text = element_text(size=12,face = "bold"),
          strip.background = element_rect(
            color="transparent", fill="grey90"))
)

(Fig2_B<-as_tibble(CombDat_NLR) %>% 
    mutate(N=as.numeric(as.character(N))) %>% 
    mutate(Category=if_else(DeltaSD < 0.01, "Low Noise","High Noise")) %>% 
    mutate(YFP.MEDIAN.RELATIVE.MEAN=100*YFP.MEDIAN.RELATIVE.MEAN,
           YFP.SD.RELATIVE.MEAN=100*YFP.SD.RELATIVE.MEAN,
           YFP.MEDIAN.RELATIVE.SD=100*YFP.MEDIAN.RELATIVE.SD,
           YFP.SD.RELATIVE.SD=100*YFP.SD.RELATIVE.SD
    ) %>% 
    ggplot()+
    geom_errorbar(aes(x=YFP.MEDIAN.RELATIVE.MEAN, y=Fitness,
                      xmin= YFP.MEDIAN.RELATIVE.MEAN-1.96*YFP.MEDIAN.RELATIVE.SD/sqrt(N),
                      xmax= YFP.MEDIAN.RELATIVE.MEAN+1.96*YFP.MEDIAN.RELATIVE.SD/sqrt(N),
    ),color="gray70", show.legend = F) +
    geom_errorbar(aes(x=YFP.MEDIAN.RELATIVE.MEAN, y=Fitness,
                      ymin= Low.95,
                      ymax= High.95,
    ),color="gray70", show.legend = F) +
    geom_point(aes(x=YFP.MEDIAN.RELATIVE.MEAN,y=Fitness, 
                   fill=Category),shape=21 ,size=2.5,color="black",
               show.legend = F)+
    geom_smooth(aes(x=YFP.MEDIAN.RELATIVE.MEAN, y=Fitness,color=Category,fill=Category),method = "loess", formula = y ~ x, 
                se = TRUE, size = 1, alpha = 0.2, show.legend = F)+
    geom_hline(yintercept = 1,linetype=2, color="grey50")+
    geom_vline(xintercept = 100,linetype=2, color="grey50")+
    # facet_wrap(~factor(ENVIRONMENT, levels = c("GLUCOSE","GALACTOSE",
    #                                            "GLYCEROL","ETHANOL")))+
    facet_wrap2(~factor(ENVIRONMENT, levels = c("GLUCOSE","GALACTOSE",
                                                "GLYCEROL","ETHANOL")), ncol=4,
                strip = strip_themed(
                  background_x = elem_list_rect(
                    fill = c("#56B4E9","#009E73","#D55E00","#0072B2"),color= rep("transparent",4))
                ) 
    )+
    #scale_color_manual(values = c("green3","orange"))+
    scale_y_continuous(limits = c(0.94,1.05))+
    #scale_color_manual(values = c("#E6AB02", "#66A61E"))+
    #scale_fill_manual(values = c("#E6AB02", "#66A61E"))+
    scale_color_manual(values = c("#C76E00", "#A7DE70"))+ 
    scale_fill_manual(values = c("#C76E00", "#A7DE70"))+
    labs(x="Median Expression relative to WT (%)", y="Relative Fitness")+
    theme_classic()+
    theme(axis.title = element_text(size = 14),
          legend.title = element_text(size=14,face="bold"),
          legend.background = element_rect(fill = "transparent"),
          legend.position=c(0.1, 0.85),
          axis.text = element_text(size=14),
          legend.text= element_text(size=14),
          panel.spacing.x =unit(0.15, "lines") , 
          panel.spacing.y=unit(0.15,"lines"),
          strip.text = element_text(size=12,face = "bold"),
          strip.background = element_rect(
            color="transparent", fill="grey90"))
  
)

(Fig2 <-
    ggdraw() +
    draw_plot(Fig2_B,x = 0, y = 0,width = 1,height = 0.5)+ 
    draw_plot(Fig2_A,x = 0, y = 0.5,width = 1,height = 0.5)
  #draw_label(x=0.02,y=0.985,label = "A)", color = "black", size = 16, fontface = "bold")+
  #draw_label(x=0.02,y=0.487,label = "B)", color = "black", size = 16, fontface = "bold")
)

tiff("Plots/Fig2_part2.tiff",width=12,height =8,units="in",res=300)
Fig2
dev.off()

# Figure 3
# Define the optimum as reduction of 0.5% fitness relative to WT
# Since it is a non-linear curve first we listed 10 data points at a given enviornment 
# and chose the least expression one as an optimum threshold
max_fitness<-c(max(Concat_Model_Fitness$FITNESS[Concat_Model_Fitness$ENVIRONMENT=="GLUCOSE"]),
               max(Concat_Model_Fitness$FITNESS[Concat_Model_Fitness$ENVIRONMENT=="GALACTOSE"]),
               max(Concat_Model_Fitness$FITNESS[Concat_Model_Fitness$ENVIRONMENT=="GLYCEROL"]),
               max(Concat_Model_Fitness$FITNESS[Concat_Model_Fitness$ENVIRONMENT=="ETHANOL"])
)

max_exp <-c(1.288223,1.0356148,1.455407,1.714618)
Optimum<-rbind(Concat_Model_Fitness[Concat_Model_Fitness$ENVIRONMENT=="GLUCOSE",][order(abs(Concat_Model_Fitness$FITNESS[Concat_Model_Fitness$ENVIRONMENT=="GLUCOSE"]
                                                                                            - (max_fitness[1] -max_fitness[1]*0.005)))[c(1:10)],],
               Concat_Model_Fitness[Concat_Model_Fitness$ENVIRONMENT=="GALACTOSE",][order(abs(Concat_Model_Fitness$FITNESS[Concat_Model_Fitness$ENVIRONMENT=="GALACTOSE"]
                                                                                              -(max_fitness[2] -max_fitness[2]*0.005)))[c(1:10)],],
               Concat_Model_Fitness[Concat_Model_Fitness$ENVIRONMENT=="GLYCEROL",][order(abs(Concat_Model_Fitness$FITNESS[Concat_Model_Fitness$ENVIRONMENT=="GLYCEROL"]
                                                                                             - (max_fitness[3] -max_fitness[3]*0.005)))[c(1:10)],],
               Concat_Model_Fitness[Concat_Model_Fitness$ENVIRONMENT=="ETHANOL",][order(abs(Concat_Model_Fitness$FITNESS[Concat_Model_Fitness$ENVIRONMENT=="ETHANOL"]
                                                                                            - (max_fitness[4] -max_fitness[4]*0.005)))[c(1:10)],]
)

#This with summarized fitness from Mo
#Optimum<-Optimum[c(1,14,22,31),]
Optimum<-Optimum[c(9,18,23,39),] #get the minimum expression
Optimum$Threshold<-(max_fitness- max_fitness*0.005)
Optimum$MaxFitness<-max_fitness
Optimum$MaxExp<-max_exp

(Fig3A<-as_tibble(CombDat_NLR) %>% 
    mutate(N=as.numeric(as.character(N))) %>% 
    mutate(Category=if_else(DeltaSD < 0.01, "Low Noise","High Noise")) %>% 
    mutate(YFP.MEDIAN.RELATIVE.MEAN=100*YFP.MEDIAN.RELATIVE.MEAN,
           YFP.SD.RELATIVE.MEAN=100*YFP.SD.RELATIVE.MEAN,
           YFP.MEDIAN.RELATIVE.SD=100*YFP.MEDIAN.RELATIVE.SD,
           YFP.SD.RELATIVE.SD=100*YFP.SD.RELATIVE.SD
    ) %>% 
    ggplot()+
    geom_errorbar(aes(x=YFP.MEDIAN.RELATIVE.MEAN, y=Fitness,
                      xmin= YFP.MEDIAN.RELATIVE.MEAN-1.96*YFP.MEDIAN.RELATIVE.SD/sqrt(N),
                      xmax= YFP.MEDIAN.RELATIVE.MEAN+1.96*YFP.MEDIAN.RELATIVE.SD/sqrt(N),
    ),color="gray70") +
    geom_errorbar(aes(x=YFP.MEDIAN.RELATIVE.MEAN, y=Fitness,
                      ymin= Low.95,
                      ymax= High.95,
    ),color="gray70") +
    geom_point(aes(x=YFP.MEDIAN.RELATIVE.MEAN,y=Fitness, 
                   fill=Category),shape=21 ,size=2.5,color="black")+
    geom_line(dat=Concat_Model_Fitness,aes(x=EXPRESSION*100, y=FITNESS))+
    #geom_hline(yintercept = 1,linetype=2, color="grey50")+
    #geom_hline(yintercept = 0.995,linetype=2, color="grey50")+
    geom_hline(data=Optimum, aes(yintercept = Threshold ),linetype=2, color="grey50")+
    geom_hline(data=Optimum, aes(yintercept = MaxFitness ),linetype=2, color="grey10")+
    geom_vline(data=Optimum, aes(xintercept = MaxExp*100 ),linetype=2, color="grey50")+
    #geom_vline(xintercept = 100,linetype=2, color="grey50")+
    geom_vline(data=Optimum,aes(xintercept = EXPRESSION*100),linetype=1, color="firebrick")+
    facet_wrap2(~factor(ENVIRONMENT, levels = c("GLUCOSE","GALACTOSE",
                                                "GLYCEROL","ETHANOL")), ncol=4,
                strip = strip_themed(
                  background_x = elem_list_rect(
                    fill = c("#56B4E9","#009E73","#D55E00","#0072B2"),color= rep("transparent",4))
                ) 
    )+
    scale_y_continuous(limits = c(0.94,1.05))+
    # scale_color_manual(values = c("#E6AB02", "#66A61E"))+
    # scale_fill_manual(values = c("#E6AB02", "#66A61E"))+
    scale_color_manual(values = c("#C76E00", "#A7DE70"))+ 
    scale_fill_manual(values = c("#C76E00", "#A7DE70"))+
    labs(x="Median Expression relative to WT (%)", y="Relative Fitness")+
    theme_classic()+
    theme(axis.title = element_text(size = 16,color="black"),
          #legend.title = element_text(size=14,face="bold"),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.position=c(0.05, 0.80),
          axis.text = element_text(size=16,color="black"),
          legend.text= element_text(size=14),
          panel.spacing.x =unit(0.15, "lines") , 
          panel.spacing.y=unit(0.15,"lines"),
          strip.text = element_text(size=13,face = "bold"),
          strip.background = element_rect(
            color="transparent", fill="grey90"))
)

#Figure 3
#####3D-E Plot close and far from optimum class 
CombDat_NLR<-CombDat_NLR %>% 
  mutate( OPT_Group = case_when(
    ENVIRONMENT=="GLUCOSE" & 
      YFP.MEDIAN.RELATIVE.MEAN < Optimum$EXPRESSION[(Optimum$ENVIRONMENT=="GLUCOSE")] ~ "Far from optimum",
    ENVIRONMENT=="GALACTOSE" & 
      YFP.MEDIAN.RELATIVE.MEAN < Optimum$EXPRESSION[(Optimum$ENVIRONMENT=="GALACTOSE")] ~ "Far from optimum",
    ENVIRONMENT=="GLYCEROL" & 
      YFP.MEDIAN.RELATIVE.MEAN < Optimum$EXPRESSION[(Optimum$ENVIRONMENT=="GLYCEROL")] ~ "Far from optimum",
    ENVIRONMENT=="ETHANOL" & 
      YFP.MEDIAN.RELATIVE.MEAN < Optimum$EXPRESSION[(Optimum$ENVIRONMENT=="ETHANOL")] ~ "Far from optimum",
    .default = "Close to optimum")
  )

(Fig3B<-as_tibble(CombDat_NLR) %>% 
    mutate(N=as.numeric(as.character(N))) %>% 
    mutate(Category=if_else(DeltaSD < 0.01, "Delta-SD<1%","Delta-SD>1%")) %>% 
    mutate(YFP.MEDIAN.RELATIVE.MEAN=100*YFP.MEDIAN.RELATIVE.MEAN,
           YFP.SD.RELATIVE.MEAN=100*YFP.SD.RELATIVE.MEAN,
           YFP.MEDIAN.RELATIVE.SD=100*YFP.MEDIAN.RELATIVE.SD,
           YFP.SD.RELATIVE.SD=100*YFP.SD.RELATIVE.SD
    ) %>% 
    ggplot()+
    geom_point(aes(x=DeltaSD,y=DeltaFitness,fill=YFP.MEDIAN.RELATIVE.MEAN),
               shape=21 ,size=2.5,color="black", show.legend = T)+
    stat_cor(aes(x=DeltaSD,y=DeltaFitness),label.y = 0.02,
             cor.coef.name = c("r"),size=5)+ 
    geom_smooth(aes(x=DeltaSD,y=DeltaFitness),method='lm')+
    geom_hline(yintercept = 0,linetype=2, color="grey50")+
    geom_vline(xintercept = 0,linetype=2, color="grey50")+
    scale_fill_gradient2()+
    facet_wrap2(OPT_Group~factor(ENVIRONMENT, levels = c("GLUCOSE","GALACTOSE",
                                                         "GLYCEROL","ETHANOL")), ncol=4,
                scales = "free_x",
                strip = strip_themed(
                  background_x = elem_list_rect(
                    fill = c("#56B4E9","#009E73","#D55E00","#0072B2"),color= rep("transparent",4))
                ) 
    )+
    # facet_wrap(OPT_Group~factor(ENVIRONMENT, levels = c("GLUCOSE","GALACTOSE",
    #                                                     "GLYCEROL","ETHANOL")),
    #            nrow = 2,scales = "free_x")+
    labs(x="Delta Noise (Standard deviation)", y="Delta Fitness", fill="Median Expression")+
    theme_classic()+
    theme(axis.title = element_text(size = 16, color="black"),
          legend.title = element_text(size=12,face="bold"),
          legend.background = element_rect(fill = "transparent"),
          legend.position=c(0.05, 0.1),
          axis.text = element_text(size=16,color="black"),
          legend.text= element_text(size=12),
          panel.spacing.x =unit(0.15, "lines") , 
          panel.spacing.y=unit(0.15,"lines"),
          strip.text = element_text(size=13,face = "bold"),
          strip.background = element_rect(
            color="transparent", fill="grey90"))
)

# (Fig3 <-
#     ggdraw() +
#     draw_plot(Fig3B,x = 0, y = 0,width = 1,height = 0.7)+ 
#     draw_plot(Fig3A,x = 0, y = 0.7,width = 1,height = 0.3)+
#     draw_label(x=0.04,y=0.985,label = "A)", color = "black", size = 16, fontface = "bold")+
#     draw_label(x=0.04,y=0.7,label = "B)", color = "black", size = 16, fontface = "bold")+
#     draw_label(x=0.04,y=0.35,label = "C)", color = "black", size = 16, fontface = "bold")
# )

(Fig3 <-
    ggdraw() +
    draw_plot(Fig3B,x = 0, y = 0,width = 1,height = 0.7)+ 
    draw_plot(Fig3A,x = 0, y = 0.7,width = 1,height = 0.3)+
    draw_label(x=0.07,y=0.98,label = "A", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.31,y=0.98,label = "B", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.54,y=0.98,label = "C", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.78,y=0.98,label = "D", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.07,y=0.67,label = "E", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.31,y=0.67,label = "F", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.54,y=0.67,label = "G", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.78,y=0.67,label = "H", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.07,y=0.33,label = "I", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.31,y=0.33,label = "J", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.54,y=0.33,label = "K", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.78,y=0.33,label = "L", color = "black", size = 16, fontface = "bold")
  #draw_label(x=0.04,y=0.7,label = "B)", color = "black", size = 16, fontface = "bold")+
  #draw_label(x=0.04,y=0.35,label = "C)", color = "black", size = 16, fontface = "bold")
)

tiff("Plots/Fig3.tiff",width=16,height =12,units="in",res=300)
Fig3
dev.off()

png("Plots/Fig3.png",width=16,height =12,units="in",res=300)
Fig3
dev.off()

###### Supplementary Figure2 noise as FANO to support Figure 3
(Sup_Fig3A<-as_tibble(CombDat_NLR) %>% 
    mutate(N=as.numeric(as.character(N))) %>% 
    mutate(Category=if_else(DeltaFANO < 0.01, "Low Noise","High Noise")) %>% 
    mutate(YFP.MEDIAN.RELATIVE.MEAN=100*YFP.MEDIAN.RELATIVE.MEAN,
           YFP.FANO.RELATIVE.MEAN=100*YFP.FANO.RELATIVE.MEAN,
           YFP.MEDIAN.RELATIVE.SD=100*YFP.MEDIAN.RELATIVE.SD,
           YFP.FANO.RELATIVE.SD=100*YFP.FANO.RELATIVE.SD
    ) %>% 
    ggplot()+
    geom_errorbar(aes(x=YFP.MEDIAN.RELATIVE.MEAN, y=Fitness,
                      xmin= YFP.MEDIAN.RELATIVE.MEAN-1.96*YFP.MEDIAN.RELATIVE.SD/sqrt(N),
                      xmax= YFP.MEDIAN.RELATIVE.MEAN+1.96*YFP.MEDIAN.RELATIVE.SD/sqrt(N),
    ),color="gray70") +
    geom_errorbar(aes(x=YFP.MEDIAN.RELATIVE.MEAN, y=Fitness,
                      ymin= Low.95,
                      ymax= High.95,
    ),color="gray70") +
    geom_point(aes(x=YFP.MEDIAN.RELATIVE.MEAN,y=Fitness, 
                   fill=Category),shape=21 ,size=2.5,color="black")+
    geom_line(dat=Concat_Model_Fitness,aes(x=EXPRESSION*100, y=FITNESS))+
    geom_hline(data=Optimum, aes(yintercept = Threshold ),linetype=2, color="grey50")+
    geom_hline(data=Optimum, aes(yintercept = MaxFitness ),linetype=2, color="grey10")+
    geom_vline(data=Optimum, aes(xintercept = MaxExp*100 ),linetype=2, color="grey50")+
    #geom_vline(xintercept = 100,linetype=2, color="grey50")+
    geom_vline(data=Optimum,aes(xintercept = EXPRESSION*100),linetype=1, color="firebrick")+    facet_wrap2(~factor(ENVIRONMENT, levels = c("GLUCOSE","GALACTOSE",
                                                                                                                                            "GLYCEROL","ETHANOL")), ncol=4,
                                                                                                            strip = strip_themed(
                                                                                                              background_x = elem_list_rect(
                                                                                                                fill = c("#56B4E9","#009E73","#D55E00","#0072B2"),color= rep("transparent",4))
                                                                                                            ) 
    )+
    scale_y_continuous(limits = c(0.94,1.05))+
    # scale_color_manual(values = c("#E6AB02", "#66A61E"))+
    # scale_fill_manual(values = c("#E6AB02", "#66A61E"))+
    scale_color_manual(values = c("#C76E00", "#A7DE70"))+ 
    scale_fill_manual(values = c("#C76E00", "#A7DE70"))+
    labs(x="Median Expression relative to WT (%)", y="Relative Fitness")+
    theme_classic()+
    theme(axis.title = element_text(size = 16,color="black"),
          #legend.title = element_text(size=14,face="bold"),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.position=c(0.05, 0.8),
          axis.text = element_text(size=16,color="black"),
          legend.text= element_text(size=14),
          panel.spacing.x =unit(0.15, "lines") , 
          panel.spacing.y=unit(0.15,"lines"),
          strip.text = element_text(size=13,face = "bold"),
          strip.background = element_rect(
            color="transparent", fill="grey90"))
)

(Sup_Fig3B<-as_tibble(CombDat_NLR) %>% 
    mutate(N=as.numeric(as.character(N))) %>% 
    mutate(Category=if_else(DeltaFANO < 0.01, "Delta-SD<1%","Delta-SD>1%")) %>% 
    mutate(YFP.MEDIAN.RELATIVE.MEAN=100*YFP.MEDIAN.RELATIVE.MEAN,
           YFP.SD.RELATIVE.MEAN=100*YFP.SD.RELATIVE.MEAN,
           YFP.MEDIAN.RELATIVE.SD=100*YFP.MEDIAN.RELATIVE.SD,
           YFP.SD.RELATIVE.SD=100*YFP.SD.RELATIVE.SD
    ) %>% 
    ggplot()+
    geom_point(aes(x=DeltaFANO,y=DeltaFitness,fill=YFP.MEDIAN.RELATIVE.MEAN),
               shape=21 ,size=2.5,color="black", show.legend = T)+
    stat_cor(aes(x=DeltaFANO,y=DeltaFitness),label.y = 0.02, 
             cor.coef.name = c("r"),size=5)+ 
    geom_smooth(aes(x=DeltaFANO,y=DeltaFitness),method='lm')+
    geom_hline(yintercept = 0,linetype=2, color="grey50")+
    geom_vline(xintercept = 0,linetype=2, color="grey50")+
    scale_fill_gradient2()+
    facet_wrap2(OPT_Group~factor(ENVIRONMENT, levels = c("GLUCOSE","GALACTOSE",
                                                         "GLYCEROL","ETHANOL")), ncol=4,
                scales = "free_x",
                strip = strip_themed(
                  background_x = elem_list_rect(
                    fill = c("#56B4E9","#009E73","#D55E00","#0072B2"),color= rep("transparent",4))
                ) 
    )+
    # facet_wrap(OPT_Group~factor(ENVIRONMENT, levels = c("GLUCOSE","GALACTOSE",
    #                                                     "GLYCEROL","ETHANOL")),
    #            nrow = 2,scales = "free_x")+
    labs(x="Delta Noise (Fano Factor)", y="Delta Fitness", fill="Median Expression")+
    theme_classic()+
    theme(axis.title = element_text(size = 16, color="black"),
          legend.title = element_text(size=12,face="bold"),
          legend.background = element_rect(fill = "transparent"),
          legend.position=c(0.05, 0.1),
          axis.text = element_text(size=16,color="black"),
          legend.text= element_text(size=12),
          panel.spacing.x =unit(0.15, "lines") , 
          panel.spacing.y=unit(0.15,"lines"),
          strip.text = element_text(size=13,face = "bold"),
          strip.background = element_rect(
            color="transparent", fill="grey90"))
)

(SupFig_S3 <-
    ggdraw() +
    draw_plot(Sup_Fig3B,x = 0, y = 0,width = 1,height = 0.7)+ 
    draw_plot(Sup_Fig3A,x = 0, y = 0.7,width = 1,height = 0.3)+
    draw_label(x=0.07,y=0.98,label = "A", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.31,y=0.98,label = "B", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.54,y=0.98,label = "C", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.78,y=0.98,label = "D", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.07,y=0.67,label = "E", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.31,y=0.67,label = "F", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.54,y=0.67,label = "G", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.78,y=0.67,label = "H", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.07,y=0.33,label = "I", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.31,y=0.33,label = "J", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.54,y=0.33,label = "K", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.78,y=0.33,label = "L", color = "black", size = 16, fontface = "bold")
)

tiff("Plots/SupFig_S3.tiff",width=16,height =12,units="in",res=300)
SupFig_S3
dev.off()

png("Plots/SupFig_S3.png",width=16,height =12,units="in",res=300)
SupFig_S3
dev.off()

###### Supplementary Figure4 noise as CV to support Figure 3
(Sup_Fig4A<-as_tibble(CombDat_NLR) %>% 
    mutate(N=as.numeric(as.character(N))) %>% 
    mutate(Category=if_else(DeltaCV < 0.01, "Low Noise","High Noise")) %>% 
    mutate(YFP.MEDIAN.RELATIVE.MEAN=100*YFP.MEDIAN.RELATIVE.MEAN,
           YFP.FANO.RELATIVE.MEAN=100*YFP.FANO.RELATIVE.MEAN,
           YFP.MEDIAN.RELATIVE.SD=100*YFP.MEDIAN.RELATIVE.SD,
           YFP.FANO.RELATIVE.SD=100*YFP.FANO.RELATIVE.SD
    ) %>% 
    ggplot()+
    geom_errorbar(aes(x=YFP.MEDIAN.RELATIVE.MEAN, y=Fitness,
                      xmin= YFP.MEDIAN.RELATIVE.MEAN-1.96*YFP.MEDIAN.RELATIVE.SD/sqrt(N),
                      xmax= YFP.MEDIAN.RELATIVE.MEAN+1.96*YFP.MEDIAN.RELATIVE.SD/sqrt(N),
    ),color="gray70") +
    geom_errorbar(aes(x=YFP.MEDIAN.RELATIVE.MEAN, y=Fitness,
                      ymin= Low.95,
                      ymax= High.95,
    ),color="gray70") +
    geom_point(aes(x=YFP.MEDIAN.RELATIVE.MEAN,y=Fitness, 
                   fill=Category),shape=21 ,size=2.5,color="black")+
    geom_line(dat=Concat_Model_Fitness,aes(x=EXPRESSION*100, y=FITNESS))+
    geom_hline(data=Optimum, aes(yintercept = Threshold ),linetype=2, color="grey50")+
    geom_hline(data=Optimum, aes(yintercept = MaxFitness ),linetype=2, color="grey10")+
    geom_vline(data=Optimum, aes(xintercept = MaxExp*100 ),linetype=2, color="grey50")+
    #geom_vline(xintercept = 100,linetype=2, color="grey50")+
    geom_vline(data=Optimum,aes(xintercept = EXPRESSION*100),linetype=1, color="firebrick")+    facet_wrap2(~factor(ENVIRONMENT, levels = c("GLUCOSE","GALACTOSE",
                                                                                                                                            "GLYCEROL","ETHANOL")), ncol=4,
                                                                                                            strip = strip_themed(
                                                                                                              background_x = elem_list_rect(
                                                                                                                fill = c("#56B4E9","#009E73","#D55E00","#0072B2"),color= rep("transparent",4))
                                                                                                            ) 
    )+
    scale_y_continuous(limits = c(0.94,1.05))+
    # scale_color_manual(values = c("#E6AB02", "#66A61E"))+
    # scale_fill_manual(values = c("#E6AB02", "#66A61E"))+
    scale_color_manual(values = c("#C76E00", "#A7DE70"))+ 
    scale_fill_manual(values = c("#C76E00", "#A7DE70"))+
    labs(x="Median Expression relative to WT (%)", y="Relative Fitness")+
    theme_classic()+
    theme(axis.title = element_text(size = 16,color="black"),
          #legend.title = element_text(size=14,face="bold"),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.position=c(0.05, 0.8),
          axis.text = element_text(size=16,color="black"),
          legend.text= element_text(size=14),
          panel.spacing.x =unit(0.15, "lines") , 
          panel.spacing.y=unit(0.15,"lines"),
          strip.text = element_text(size=13,face = "bold"),
          strip.background = element_rect(
            color="transparent", fill="grey90"))
)

(Sup_Fig4B<-as_tibble(CombDat_NLR) %>% 
    mutate(N=as.numeric(as.character(N))) %>% 
    mutate(Category=if_else(DeltaCV < 0.01, "Delta-SD<1%","Delta-SD>1%")) %>% 
    mutate(YFP.MEDIAN.RELATIVE.MEAN=100*YFP.MEDIAN.RELATIVE.MEAN,
           YFP.SD.RELATIVE.MEAN=100*YFP.SD.RELATIVE.MEAN,
           YFP.MEDIAN.RELATIVE.SD=100*YFP.MEDIAN.RELATIVE.SD,
           YFP.SD.RELATIVE.SD=100*YFP.SD.RELATIVE.SD
    ) %>% 
    ggplot()+
    geom_point(aes(x=DeltaCV,y=DeltaFitness,fill=YFP.MEDIAN.RELATIVE.MEAN),
               shape=21 ,size=2.5,color="black", show.legend = T)+
    stat_cor(aes(x=DeltaCV,y=DeltaFitness),label.y = 0.02, 
             cor.coef.name = c("r"),size=5)+ 
    geom_smooth(aes(x=DeltaCV,y=DeltaFitness),method='lm')+
    geom_hline(yintercept = 0,linetype=2, color="grey50")+
    geom_vline(xintercept = 0,linetype=2, color="grey50")+
    scale_fill_gradient2()+
    facet_wrap2(OPT_Group~factor(ENVIRONMENT, levels = c("GLUCOSE","GALACTOSE",
                                                         "GLYCEROL","ETHANOL")), ncol=4,
                scales = "free_x",
                strip = strip_themed(
                  background_x = elem_list_rect(
                    fill = c("#56B4E9","#009E73","#D55E00","#0072B2"),color= rep("transparent",4))
                ) 
    )+
    # facet_wrap(OPT_Group~factor(ENVIRONMENT, levels = c("GLUCOSE","GALACTOSE",
    #                                                     "GLYCEROL","ETHANOL")),
    #            nrow = 2,scales = "free_x")+
    labs(x="Delta Noise (Coefficient of Variation)", y="Delta Fitness", fill="Median Expression")+
    theme_classic()+
    theme(axis.title = element_text(size = 16, color="black"),
          legend.title = element_text(size=12,face="bold"),
          legend.background = element_rect(fill = "transparent"),
          legend.position=c(0.05, 0.1),
          axis.text = element_text(size=16,color="black"),
          legend.text= element_text(size=12),
          panel.spacing.x =unit(0.15, "lines") , 
          panel.spacing.y=unit(0.15,"lines"),
          strip.text = element_text(size=13,face = "bold"),
          strip.background = element_rect(
            color="transparent", fill="grey90"))
)

(SupFig_S4 <-
    ggdraw() +
    draw_plot(Sup_Fig4B,x = 0, y = 0,width = 1,height = 0.7)+ 
    draw_plot(Sup_Fig4A,x = 0, y = 0.7,width = 1,height = 0.3)+
    draw_label(x=0.07,y=0.98,label = "A", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.31,y=0.98,label = "B", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.54,y=0.98,label = "C", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.78,y=0.98,label = "D", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.07,y=0.67,label = "E", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.31,y=0.67,label = "F", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.54,y=0.67,label = "G", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.78,y=0.67,label = "H", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.07,y=0.33,label = "I", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.31,y=0.33,label = "J", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.54,y=0.33,label = "K", color = "black", size = 16, fontface = "bold")+
    draw_label(x=0.78,y=0.33,label = "L", color = "black", size = 16, fontface = "bold")
)

tiff("Plots/SupFig_S4.tiff",width=16,height =12,units="in",res=300)
SupFig_S4
dev.off()

png("Plots/SupFig_S4.png",width=16,height =12,units="in",res=300)
SupFig_S4
dev.off()

## Write final Fitness estimates
CombDat_NLR %>% 
  select(ID, STRAIN,MUTATION,N,STRAIN.Fitness,
         ENVIRONMENT,YFP.MEDIAN.RELATIVE.MEAN,YFP.MEDIAN.RELATIVE.SD,
         YFP.SD.RELATIVE.MEAN,YFP.SD.RELATIVE.SD,
         YFP.FANO.RELATIVE.MEAN, YFP.FANO.RELATIVE.SD,
         YFP.CV.RELATIVE.MEAN,YFP.CV.RELATIVE.SD,
         Fitness,Low.95,High.95,
         Category,DeltaSD,DeltaFANO,DeltaCV, DeltaFitness, Class) %>% 
  write_csv("Results/Supplementary_Table10.csv")