for (i in seq(0,2,0.1)) {
  #for (j in c(0.05,0.1,0.15, 0.2, 0.3, 0.4,0.5, 0.6,0.8, 1, 1.2, 1.4)) {
  for (j in c(0.05,0.1,0.15, 0.2, 0.3, 0.4,0.5, 0.6,0.8)) {
    for (k in seq(0,1,0.25)) 
      {
      
      cat("python PopsimExpFitDeapConstantNoise.py  --EXPR_MEAN_A ", i, " --EXPR_SD_A ", j,
         "--EXPR_MEAN_B 0.5 --EXPR_SD_B 0.05 ",
          " --FIT_var1 0.5 --FIT_var1_pair 0.005 --FIT_var2 0.25 --FIT_var2_pair 0.5 --weight 0.6 --h " , k, 
          " --pop_size 10000 --total_time 600 --dt 200 --iterations 5 --fitness_function mixednorm ", 
          " --output_dir /home/path_to_output_dir/ ",
          "--outfile ", paste("SIM_CompD_mean_" ,sprintf("%.2f",i),"_sd_", sprintf("%.2f",j),
                                 "_her_",sprintf("%.2f",k),"_mixednorm.csv",sep = ""),
          "\n") 
    }
  }
}

