avg_mat <- NULL
track_mat <- NULL
neut_avgs <- NULL
deme_plot <- 2000
sim <- 5
demes <- seq(1:1000)
gen <- 180
z = NULL
setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/2d_model/diploidization_2d-dom_meiotic_dipLambda-20_recessive_2dsize-10_K-50_m-0.05_r-0.693147_u_r-5.0e-08_u_d-5.0e-09_rho-5.0e-08_bs-0.0_ds--0.001472_g-999999_start-2501")
#dirs <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/allo_DFE_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501", full.names = TRUE)
file <- "_meanFitnessScaled"

i <- 1
dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
#pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
dat_1 <- dat[(gen/10)+1,]
dat_1 <- unname(dat_1)
dat_2d <- matrix(nrow=400,ncol=10) 
#dat_2d <- data.frame(unname(dat_1[1,1:5]))
counter <- 1
for(i in seq(1,2000,10)){
  dat_2d[counter,1:10] <- as.numeric(dat_1[1,i:(i+9)])
  counter = counter+1
}


#pops <- pops[1:(deme_plot/10),demes]
#neut <- read.csv(paste(c(i,"/",i,"_meanNeutralMutations",".csv"),sep="", collapse=""), row.names = 1)
#neut <- neut[1:(deme_plot/10),demes]
#track_mat_temp <- (pops != 0)*1
#if(i == 1){
#  avg_mat <- dat
#  track_mat <- track_mat_temp
#  #neut_avgs <- neut
#}
#else{
  #avg_mat <- avg_mat + dat
 # track_mat <- track_mat + track_mat_temp
 # #neut_avgs <- neut_avgs + neut
#}

avgs_mat <- avg_mat / track_mat
#neut_mat <- neut_avgs / track_mat
#diff_mat <- avgs_mat - neut_mat

library(fields)
x = seq(1:10)
y = seq(1:50)
z = dat_2d[1:50,]
png(file=paste(c("2d_",gen,file,".png"),sep="", collapse=""), width = 700, height = 1100)
#image.plot(x,y,as.matrix(z), zlim=c(0,200), xlab="Generation", ylab="Deme")
#image.plot(x,y,as.matrix(z), zlim=c(0,50), xlab="Generation", ylab="Deme")
#image.plot(x,y,as.matrix(z), zlim=c(1e-07, 0.0001), xlab="Generation", ylab="Deme")
#image.plot(x,y,as.matrix(z), zlim=c(0, 0.3), xlab="Generation", ylab="Deme")
#image.plot(x,y,as.matrix(z), zlim=c(0,75), xlab="Generation", ylab="Deme")
image.plot(x,y,t(as.matrix(z)), zlim=c(0.8,1.05), xlab="Deme Position X", ylab="Deme Position Y")
#image.plot(x,y,as.matrix(z), zlim=c(0, 0.05), xlab="Generation", ylab="Deme")
#image.plot(x,y,as.matrix(z), zlim=c(0.000000001, 0.0001), xlab="Generation", ylab="Deme")
#image.plot(x,y,as.matrix(z), zlim=c(-5, 5), xlab="Generation", ylab="Deme")
#image.plot(x,y,t(as.matrix(z)), zlim=c(0.5,1), xlab="Deme Position X", ylab="Deme Position Y")


dev.off()
#}

