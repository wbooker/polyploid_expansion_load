avg_mat <- NULL
track_mat <- NULL
deme_plot <- 2000
sim <- 5
setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/test_12_14_2022/auto_DFE_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-2501/")
dirs <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/test_12_14_2022/auto_DFE_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-2501/", full.names = TRUE)
file <- "_meanFitness"
#for(dir in dirs[c(seq(2,5),8,seq(11,21),seq(24,32))]){
#for(dir in dirs[c(6,7)]){
dir <- dirs
  for(i in c(1,2,3,4,5)){
    dat <- read.csv(paste(c(dir,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
    dat <- dat[1:(deme_plot/10),1:850]
    track_mat_temp <- (dat != 0)*1
    if(i == 1){
      avg_mat <- dat
      track_mat <- track_mat_temp
    }
    else{
      avg_mat <- avg_mat + dat
      track_mat <- track_mat + track_mat_temp
    }
  }

  avgs_mat <- avg_mat / track_mat

  library(fields)
  x = as.integer(row.names(dat))
  y = seq(1:850)
  z = avgs_mat
  png(file=paste(c(dir,"/","all",file,".png"),sep="", collapse=""), width = 600, height = 900)
  #image.plot(x,y,as.matrix(z), zlim=c(0,200), xlab="Generation", ylab="Deme")
  #image.plot(x,y,as.matrix(z), zlim=c(0,50), xlab="Generation", ylab="Deme")
  #image.plot(x,y,as.matrix(z), zlim=c(1e-07, 0.0001), xlab="Generation", ylab="Deme")
  #image.plot(x,y,as.matrix(z), zlim=c(0, 0.2), xlab="Generation", ylab="Deme")
  #image.plot(x,y,as.matrix(z), zlim=c(0,75), xlab="Generation", ylab="Deme")
  image.plot(x,y,as.matrix(z), zlim=c(0.9,1.05), xlab="Generation", ylab="Deme")
  #image.plot(x,y,as.matrix(z), zlim=c(0, 0.05), xlab="Generation", ylab="Deme")
  #image.plot(x,y,as.matrix(z), zlim=c(0.000000001, 0.0001), xlab="Generation", ylab="Deme")
  #image.plot(x,y,as.matrix(z), zlim=c(0, 75), xlab="Generation", ylab="Deme")

  dev.off()
#}



