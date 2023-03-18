avg_mat <- NULL
track_mat <- NULL
deme_plot <- 3000
sim <- 5
setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/auto_additive_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501")
dirs <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/test_12_15_2022/auto_DFE_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-2501/", full.names = TRUE)
file <- "_log"
#for(dir in dirs[c(seq(2,5),8,seq(11,21),seq(24,32))]){
#for(dir in dirs[c(6,7)]){
dir <- dirs
  for(i in seq(1,50)){
    dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
    dat <- dat[1:(deme_plot/10),]
    track_mat_temp <- (dat$current_deme != 0)*1
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
  x_auto = (as.integer(row.names(dat))-2500)
  y_auto = avgs_mat$current_deme

avg_mat <- NULL
track_mat <- NULL
deme_plot <- 3000
sim <- 5
setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/allo_additive_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501")
dirs <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/test_12_15_2022/allo_DFE_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-2501/", full.names = TRUE)
file <- "_log"
#for(dir in dirs[c(seq(2,5),8,seq(11,21),seq(24,32))]){
#for(dir in dirs[c(6,7)]){
dir <- dirs
  for(i in seq(1,50)){
    dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
    dat <- dat[1:(deme_plot/10),]
    track_mat_temp <- (dat$current_deme != 0)*1
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
  x_allo = (as.integer(row.names(dat))-2500)
  y_allo = avgs_mat$current_deme


avg_mat <- NULL
track_mat <- NULL
deme_plot <- 3000
sim <- 5
setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/diploid_additive_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501")
dirs <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/test_12_15_2022/diploid_DFE_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-2501/", full.names = TRUE)
file <- "_log"
#for(dir in dirs[c(seq(2,5),8,seq(11,21),seq(24,32))]){
#for(dir in dirs[c(6,7)]){
dir <- dirs
  for(i in seq(1,50)){
    dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
    dat <- dat[1:(deme_plot/10),]
    track_mat_temp <- (dat$current_deme != 0)*1
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
  x_diploid = (as.integer(row.names(dat))-2500)
  y_diploid = avgs_mat$current_deme

  png(file=paste(c("demes_by_gen",".png"),sep="", collapse=""), width = 600, height = 900)
  plot((x_allo),y_allo, type="l", col="red", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
  #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
  lines(x_diploid,y_diploid,col="black",lwd=5)
  lines(x_auto,y_auto,col="blue",lwd=5)
  legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
       col=c("Black", "Blue", "Red"), lty=1, cex=0.8)
  #legend(1, 450, legend=c("Diploid", "Auto"),
 #     col=c("Black", "Blue"), lty=1, cex=0.8)

  dev.off()


############## Fitness distribution Figs ###################
allo <- read.table("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/allo_additive_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501/core_fits_2501.txt")
auto <- read.table("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/auto_additive_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501/core_fits_2501.txt")
dip <- read.table("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/diploid_additive_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501/core_fits_2501.txt")
dat <- data.frame(values = c(allo[,1],auto[,1],dip[,1]), group=c(rep("allo", 50),rep("auto", 50),rep("dip", 49)))
png(file=paste(c("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/additive_noBeneficial",".png"),sep="", collapse=""), width = 600, height = 900)
boxplot(values ~ group, dat)  
dev.off()

allo <- read.table("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/allo_additive_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0045_ds--0.0045_g-999999_start-2501/core_fits_2501.txt")
auto <- read.table("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/auto_additive_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0045_ds--0.0045_g-999999_start-2501/core_fits_2501.txt")
dip <- read.table("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/diploid_additive_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0045_ds--0.0045_g-999999_start-2501/core_fits_2501.txt")
dat <- data.frame(values = c(allo[,1],auto[,1],dip[,1]), group=c(rep("allo", 50),rep("auto", 50),rep("dip", 50)))
png(file=paste(c("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/additive_withBeneficial",".png"),sep="", collapse=""), width = 600, height = 900)
boxplot(values ~ group, dat)  
dev.off()

allo <- read.table("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/allo_DFE_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501/core_fits_2501.txt")
auto <- read.table("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/auto_DFE_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501/core_fits_2501.txt")
dip <- read.table("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/diploid_DFE_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501/core_fits_2501.txt")
dat <- data.frame(values = c(allo[,1],auto[,1],dip[,1]), group=c(rep("allo", 50),rep("auto", 50),rep("dip", 50)))
png(file=paste(c("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/DFE_core",".png"),sep="", collapse=""), width = 600, height = 900)
boxplot(values ~ group, dat)  
dev.off()

auto <- read.table("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/auto_recessive_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501/core_fits_2501.txt")
dip <- read.table("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/diploid_recessive_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501/core_fits_2501.txt")
dat <- data.frame(values = c(auto[,1],dip[,1]), group=c(rep("auto", 50),rep("dip", 50)))
png(file=paste(c("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/recessive_noBeneficial_core",".png"),sep="", collapse=""), width = 600, height = 900)
boxplot(values ~ group, dat)  
dev.off()

auto <- read.table("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/auto_recessive_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0045_ds--0.0045_g-999999_start-2501/core_fits_2501.txt")
dip <- read.table("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/diploid_recessive_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0045_ds--0.0045_g-999999_start-2501/core_fits_2501.txt")
dat <- data.frame(values = c(auto[,1],dip[,1]), group=c(rep("auto", 50),rep("dip", 50)))
png(file=paste(c("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/recessive_Beneficial_core",".png"),sep="", collapse=""), width = 600, height = 900)
boxplot(values ~ group, dat)  
dev.off()