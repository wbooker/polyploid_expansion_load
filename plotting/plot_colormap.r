library(fields)
library(scales)
avg_mat <- NULL
track_mat <- NULL
neut_avgs <- NULL
deme_plot <- 1000
sim <- 100
demes <- seq(1:600)
nrep <- 15
z = NULL
#plot_colors = viridis(64,option = "C")
#mod_dirs <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/pair_fit_dip_FIXED/test_recomb/no_s", full.names = TRUE, recursive = FALSE)
#setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/test_diploidization_3_17_2023/diploidization-dom_dipLambda-1_remDipMuts-1_additive_K-50_m-0.05_r-0.693147_u_r-5.0e-08_u_d-5.0e-10_rho-5.0e-08_bs-0.0_ds--0.001472_g-999999_start-2501")
#dirs <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/full_run_12_22_2022/allo_DFE_K-50_m-0.05_r-0.693147_u-2.5e-08_rho-2.5e-08_bs-0.0_ds--0.0045_g-999999_start-2501", full.names = TRUE)
#file <- "_fixedMutations"

##Note here just to remember to get presence for averaging from popSize, not from if value is 0 or not, since some measurements like heterozygosity will be 0 and that will mess up the averaging as it is currently

stats <- c("_dipIndex")
starts <- c(0)
ends <- c(1)
for(stat in 1:length(stats)){
  file <- stats[stat]
mod_dir <- list.dirs(path = "/work/users/w/w/wwbooker/polyploid_expansion_load/output/test_dip", full.names = TRUE, recursive = FALSE)

for(mod in mod_dir[c(4)]){
  #mod <- "/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/basic/auto_additive_fixed_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.0_ds--0.005_g-999999_start-5001"
  #for(mod in mod_dirs){
    setwd(mod)
   init_count = 1
  cur_rep <- 1
#    for(i in c(1:11,13:18,20,23:25)){
  for(i in c(1:15)){
      if(cur_rep > nrep){
        next
      }
      dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
      dat <- dat[1:(deme_plot/10),demes]
      pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
      pops <- pops[1:(deme_plot/10),demes]
      if(length(na.omit(pops[,1])) < 100){
        next
      }
      print(c(file,i))

      #neut <- read.csv(paste(c(i,"/",i,"_meanNeutralMutations",".csv"),sep="", collapse=""), row.names = 1)
      #neut <- neut[1:(deme_plot/10),demes]
      track_mat_temp <- (pops != 0)*1
      if(init_count == 1){
        init_count = 0
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
    #z[z>0.3] = 0.3
    png(file=paste(c(i,"/",i,i,file,".png"),sep="", collapse=""))
    #image.plot(x,y,as.matrix(z), zlim=c(starts[i],ends[i]), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(1,70), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(1e-07, 0.001), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0, 0.01), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0,75), xlab="Generation", ylab="Deme")
    image.plot(x,y,as.matrix(z), zlim=c(0,1), xlab="Generation", ylab="Deme", las = 1)
    #image.plot(x,y,as.matrix(z), zlim=c(0, 0.05), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0.000000001, 0.0001), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(-5, 5), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0, 1), xlab="Generation", ylab="Deme", col = plot_colors)


    dev.off()
    pdf(file=paste(c(i,"/",i,i,file,".pdf"),sep="", collapse=""))
    image.plot(x,y,as.matrix(z), zlim=c(0,1), xlab="Generation", ylab="Deme", las = 1)
    dev.off()
    cur_rep <- cur_rep + 1
  
    #}
    }
    t <- track_mat
    #t[which(t < 75)] <- 0
    avgs_mat <- avg_mat / t
    #x = as.integer(row.names(dat))
    y = demes[1:200]
    z = avgs_mat[,1:200]
    #z[z>0.3] = 0.3
    png(file=paste(c("all",file,".png"),sep="", collapse=""))
    #image.plot(x,y,as.matrix(z), zlim=c(starts[i],ends[i]), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(1,70), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(1e-07, 0.001), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0, 0.01), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0,75), xlab="Generation", ylab="Deme")
    image.plot(x,y,as.matrix(z), zlim=c(0,1), xlab="Generation", ylab="Deme", las = 1)
    #image.plot(x,y,as.matrix(z), zlim=c(0, 0.05), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0.000000001, 0.0001), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(-5, 5), xlab="Generation", ylab="Deme")
    #image.plot(x,y,as.matrix(z), zlim=c(0, 1), xlab="Generation", ylab="Deme", col = plot_colors)


    dev.off()
    pdf(file=paste(c("all",file,".pdf"),sep="", collapse=""))
    image.plot(x,y,as.matrix(z), zlim=c(0,1), xlab="Generation", ylab="Deme", las = 1)
    dev.off()

  }
}
#}


#stats <- c("_dipIndex", "_meanFitnessScaled", "_heterozygosity","_exp_load_prop_1","_varFitness","_maxPhs","_unscaledFitness","_popSize","_meanMutationsPerInd","_meanMutationsPerInd", "_fixedMutations")
#starts <- c(0,0.7,1e-07,0,1e-07,0,0.6,1,0,0,0)
#ends <- c(1,1.05,0.0001,0.35,1e-03,0.01,1.05,70,500,300,300)