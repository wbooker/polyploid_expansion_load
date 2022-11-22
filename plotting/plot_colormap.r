avg_mat <- NULL
track_mat <- NULL
deme_plot <- 500
sim <- 1
setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/output_new/2/dip_recessive_K-50_m-0.05_r-0.693147_u-2.5e-11_rho-2.5e-11_bs-0.0_ds--0.005_g-999999999")
for(i in c(sim)){
  dat <- read.csv(paste(c(i,"_fixedMutations.csv"),sep="", collapse=""), row.names = 1)
  dat <- dat[1:deme_plot,1:850]
  track_mat_temp <- (dat != 0)*1
  if(i == sim){
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
png(file=paste(c(sim,"_fixedMutations.png"),sep="", collapse=""), width = (deme_plot+100), height = 900)
image.plot(x,y,as.matrix(z), zlim=c(0,100), xlab="Generation", ylab="Deme")
#image.plot(x,y,as.matrix(z), zlim=c(0.1,30), xlab="Generation", ylab="Deme")
#image.plot(x,y,as.matrix(z), zlim=c(0, 0.015), xlab="Generation", ylab="Deme")
#image.plot(x,y,as.matrix(z), zlim=c(0.69,1.01), xlab="Generation", ylab="Deme")
#image.plot(x,y,as.matrix(z), zlim=c(0, 0.0075), xlab="Generation", ylab="Deme")
#image.plot(x,y,as.matrix(z), zlim=c(0, 75), xlab="Generation", ylab="Deme")

dev.off()




