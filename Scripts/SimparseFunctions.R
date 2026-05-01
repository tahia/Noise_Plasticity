library(ggplot2)
library(tidyverse)

#### Function to process sims by DEAP
process_sims_comp<-function(file=file,Relmean=Relmean, RelSD=RelSD) {
  
  colNames<-c("Initial_ExpMean","Initial_ExpSD","Fitfun","FIT_var1","FIT_var2",
              "Heritability","iteration","time","mean_A","SD_A",
              "mean_B","SD_B","freq_A","freq_B","PopSize","mean_dt_A","mean_dt_B")
  
  Fitness<-read.csv(file,
                    header = FALSE, col.names = colNames  ) %>% 
    mutate(Time= 1+ (time/200)) %>% 
    mutate(logR= log(freq_A/freq_B)) %>% 
    select(Initial_ExpMean,Initial_ExpSD,Fitfun,
           Heritability,iteration,time,freq_A,freq_B,logR,Time) %>% 
    group_by(Initial_ExpMean,Initial_ExpSD,Heritability, iteration) %>% 
    nest() %>% 
    mutate(model = map(data, ~ lm(logR ~ Time, data = .x))) %>%
    mutate(tidied_model = map(model, tidy)) %>%
    unnest(tidied_model) %>%
    filter(term == "Time") %>% 
    mutate(Fitness=exp(estimate)) %>% 
    select(Initial_ExpMean,Initial_ExpSD,Heritability, iteration,Fitness) 
  
  #We already know from the fitness function the optimum
  MaxFit<-mean(Fitness$Fitness[which(Fitness$Initial_ExpMean==Relmean & 
                                       Fitness$Initial_ExpSD==RelSD)])
  
  SumFitness<-Fitness %>% 
    mutate(Fitness=Fitness/MaxFit) %>% 
    group_by(Initial_ExpMean,Initial_ExpSD,Heritability) %>% 
    summarise_all(funs(mean,sd, se =sd(.)/sqrt(n()))) %>% 
    rename(Mean=Initial_ExpMean,
           SD=Initial_ExpSD, 
           Fitness=Fitness_mean) %>%
    mutate(Noise=as.factor(100*SD/0.05) ) 
  
  return(SumFitness)
  
}