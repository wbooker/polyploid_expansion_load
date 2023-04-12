library(fields)
avg_mat <- NULL
track_mat <- NULL
neut_avgs <- NULL
deme_plot <- 2500
sim <- 5
demes <- seq(1:1000)
z = NULL
mod_dirs <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/new_fit_calc_modelTest_WeekendRun", full.names = TRUE, recursive = FALSE)
#setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/test_diploidization_3_17_2023/diploidization-dom_dipLambda-1_remDipMuts-1_additive_K-50_m-0.05_r-0.693147_u_r-5.0e-08_u_d-5.0e-10_rho-5.0e-08_bs-0.0_ds--0.001472_g-999999_start-2501")
#dirs <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/allo_DFE_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501", full.names = TRUE)
#file <- "_fixedMutations"

##Note here just to remember to get presence for averaging from popSize, not from if value is 0 or not, since some measurements like heterozygosity will be 0 and that will mess up the averaging as it is currently

stats <- c("_meanFitnessScaled")
starts <- c(0)
ends <- c(1)
for(stat in 1:length(stats)){
  file <- stats[stat]
  #mod <- "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/diploid_additive_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501"
  for(mod in mod_dirs){
    setwd(mod)
    for(i in c(1:27,29:50)){
      print(c(file,i))
      dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
      dat <- dat[1:(deme_plot/10),demes]
      pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
      pops <- pops[1:(deme_plot/10),demes]
      #neut <- read.csv(paste(c(i,"/",i,"_meanNeutralMutations",".csv"),sep="", collapse=""), row.names = 1)
      #neut <- neut[1:(deme_plot/10),demes]
      track_mat_temp <- (pops != 0)*1
      if(i == 1){
        avg_mat <- dat
        track_mat <- track_mat_temp
        #neut_avgs <- neut
      }
      else{
        avg_mat <- avg_mat + dat
        track_mat <- track_mat + track_mat_temp
    #    #neut_avgs <- neut_avgs + neut
      }

      single_mat <- dat / track_mat_temp
    #neut_mat <- neut_avgs / track_mat
    #diff_mat <- avgs_mat - neut_mat

    x = as.integer(row.names(dat))
    y = demes
    z = single_mat
    png(file=paste(c(i,"/",i,file,".png"),sep="", collapse=""), width = 700, height = 1100)
    #image.plot(x,y,as.matrix(z), zlim=c(starts[i],ends[i]), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(1,70), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(1e-07, 0.001), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0, 0.01), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0,75), xlab="Generation", ylab="Deme")
    image.plot(x,y,as.matrix(z), zlim=c(0.6,1.05), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0, 0.05), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0.000000001, 0.0001), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(-5, 5), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0, 1), xlab="Generation", ylab="Deme")


    dev.off()
    #}
    }
    avgs_mat <- avg_mat / track_mat
    x = as.integer(row.names(dat))
    y = demes
    z = avgs_mat
    png(file=paste(c("all",file,".png"),sep="", collapse=""), width = 700, height = 1100)
    #image.plot(x,y,as.matrix(z), zlim=c(starts[i],ends[i]), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(1,70), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(1e-07, 0.001), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0, 0.01), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0,75), xlab="Generation", ylab="Deme")
    image.plot(x,y,as.matrix(z), zlim=c(0.6,1.05), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0, 0.05), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0.000000001, 0.0001), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(-5, 5), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0, 1), xlab="Generation", ylab="Deme")


    dev.off()
  }
}



#stats <- c("_dipIndex", "_meanFitnessScaled", "_heterozygosity","_exp_load_prop_1","_varFitness","_maxPhs","_unscaledFitness","_popSize","_meanMutationsPerInd","_meanMutationsPerInd", "_fixedMutations")
#starts <- c(0,0.7,1e-07,0,1e-07,0,0.6,1,0,0,0)
#ends <- c(1,1.05,0.0001,0.35,1e-03,0.01,1.05,70,500,300,300)