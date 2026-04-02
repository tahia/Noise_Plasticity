#### Function to process sims by DEAP
process_sims_tidy<-function(file=file) {
  
  colNames<-c("Initial_ExpMean","Initial_ExpSD","Fitfun","FIT_var1","FIT_var2",
              "Heritability","iteration","time","mean_doubling_time","SD","Fano","CV")
  
  Dist1<-read.csv(file,
                  header = FALSE, col.names = colNames  ) %>% 
    mutate(Time= 1+ (time/480)) %>% 
    group_by(Initial_ExpMean,Initial_ExpSD,Heritability, iteration) %>% 
    nest() %>% 
    mutate(model = map(data, ~ lm(mean_doubling_time ~ Time, data = .x))) %>%
    mutate(tidied_model = map(model, tidy)) %>%
    unnest(tidied_model) %>%
    filter(term == "Time") %>% 
    select(Initial_ExpMean,Initial_ExpSD,Heritability, iteration,estimate,std.error)
  
  WTREF <-Dist1 %>% 
    dplyr::filter(Initial_ExpMean==1 & Initial_ExpSD==0.1) %>% 
    group_by(Heritability) %>% 
    summarise_all(funs(mean,sd, se =sd(.)/sqrt(n())))
  
  WTREF_s<- -(WTREF$estimate_mean[WTREF$Heritability==0])
  
  SumFitness<-Dist1 %>% 
    mutate(s= -estimate) %>% 
    mutate(Fitness = 1+s - WTREF_s) %>% 
    group_by(Initial_ExpMean,Initial_ExpSD,Heritability) %>% 
    summarise_all(funs(mean,sd, se =sd(.)/sqrt(n()))) %>% 
    rename(Mean=Initial_ExpMean,
           SD=Initial_ExpSD, 
           Fitness=Fitness_mean) %>%
    mutate(Noise=as.factor(100*SD/0.1) ) 
  
  return(SumFitness)
  
}

process_sims<-function(file=file) {
  colNames<-c("Initial_ExpMean","Initial_ExpSD","Fitfun","FIT_var1","FIT_var2",
              "Heritability","iteration","time","mean_doubling_time","SD","Fano","CV")
  
  Dist1<-read.csv(file,
                  header = FALSE, col.names = colNames  )
  
  Dist1_T1<-as_tibble(Dist1) %>% 
    dplyr::filter(time==0 ) %>% 
    select(Initial_ExpMean,Initial_ExpSD,Heritability,iteration,mean_doubling_time) %>% 
    rename(DT0=mean_doubling_time)
  
  Dist1_T1_4<-as_tibble(Dist1) %>% 
    dplyr::filter(time==1440 ) %>% 
    select(Initial_ExpMean,Initial_ExpSD,Heritability,iteration,mean_doubling_time) %>% 
    rename(DT4=mean_doubling_time) %>% 
    mutate(DT0=Dist1_T1$DT0)
  
  #WTREF: Mean 1 SD 0.1
  WTREF <-Dist1_T1_4 %>% 
    dplyr::filter(Initial_ExpMean==1 & Initial_ExpSD==0.1) %>% 
    group_by(Heritability) %>% 
    summarise_all(funs(mean,sd, se =sd(.)/sqrt(n())))
  
  WTREF_DT0<- WTREF$DT0_mean[WTREF$Heritability==0]
  WTREF_DT4<- WTREF$DT4_mean[WTREF$Heritability==0]
  
  #Fitness
  SumFitness<-Dist1_T1_4 %>% 
    #mutate(Fitness = 1-(log(DT4/WTREF_DT4 )/4)) %>% 
    mutate(s= ((WTREF_DT4/DT4) -1)/3) %>% 
    #mutate(Fitness = WTREF_DT4/DT4 ) %>% 
    mutate(Fitness = 1+s ) %>% 
    group_by(Initial_ExpMean,Initial_ExpSD,Heritability) %>% 
    summarise_all(funs(mean,sd, se =sd(.)/sqrt(n()))) %>% 
    rename(Mean=Initial_ExpMean,
           SD=Initial_ExpSD, 
           Fitness=Fitness_mean) %>%
    mutate(Noise=as.factor(100*SD/0.1) ) 
  
  return(SumFitness)
}

process_sims_nosum<-function(file=file) {
  colNames<-c("Initial_ExpMean","Initial_ExpSD","Fitfun","FIT_var1","FIT_var2",
              "Heritability","iteration","time","mean_doubling_time","SD","Fano","CV")
  
  Dist1<-read.csv(file,
                  header = FALSE, col.names = colNames  )
  
  Dist1_T1<-as_tibble(Dist1) %>% 
    dplyr::filter(time==0 ) %>% 
    select(Initial_ExpMean,Initial_ExpSD,Heritability,iteration,mean_doubling_time) %>% 
    rename(DT0=mean_doubling_time)
  
  Dist1_T1_4<-as_tibble(Dist1) %>% 
    dplyr::filter(time==1440 ) %>% 
    select(Initial_ExpMean,Initial_ExpSD,Heritability,iteration,mean_doubling_time) %>% 
    rename(DT4=mean_doubling_time) %>% 
    mutate(DT0=Dist1_T1$DT0)
  
  #WTREF: Mean 1 SD 0.1
  WTREF <-Dist1_T1_4 %>% 
    dplyr::filter(Initial_ExpMean==1 & Initial_ExpSD==0.1) %>% 
    group_by(Heritability) %>% 
    summarise_all(funs(mean,sd, se =sd(.)/sqrt(n())))
  
  WTREF_DT0<- WTREF$DT0_mean[WTREF$Heritability==0]
  WTREF_DT4<- WTREF$DT4_mean[WTREF$Heritability==0]
  
  #Fitness
  Fitness<-Dist1_T1_4 %>% 
    #mutate(Fitness = 1-(log(DT4/WTREF_DT4 )/4)) %>% 
    mutate(s= ((WTREF_DT4/DT4) -1)/3) %>% 
    #mutate(Fitness = WTREF_DT4/DT4 ) %>% 
    mutate(Fitness = 1+s ) %>% 
    mutate(Noise=as.factor(100*Initial_ExpSD/0.1) ) 
  
  return(Fitness)
}