avg_mat <- NULL
track_mat <- NULL
library(fields)

for(i in 1:5){
  dat <- read.csv(paste(c("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/allo_K-25_m-0.05_r-0.693147_u-2.5e-11_rho-2.5e-11_bs-0.0_ds--0.005_g-999999999/",i,"_fitnesses.csv"),sep="", collapse=""), row.names = 1, col.names = 1, header = 1)
  dat <- dat[1:60,1:850]
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
png(file="C:/Datamentor/R-tutorial/saving_plot2.png")
image.plot(x,y,as.matrix(z), zlim=c(min(z, na.rm = T),max(z, na.rm = T)), xlab="Generation", ylab="Deme")
dev.off()




