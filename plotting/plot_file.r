######### Plot edge value #######

nrep <- 5
deme_plot <- 2000

file_list <- c("_meanFitnessScaled")
#file_list <- c("_fixedMutations","_meanBenMutationsPerGenome","_meanBenMutationsPerInd","_unscaledFitness","_maxPhs","_heterozygosity","_fixedMutations","_meanDelMutationsPerGenome","_meanDelMutationsPerInd","_meanFitnessScaled")


for(file in file_list){
    print(file)
    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/new_fit_calc_modelTest_WeekendRun/auto_DFE_K-100_m-0.05_r-0.693147_u_del-2.5e-11_u_ben-2.5e-12_rho-2.5e-11_bs-0.0_ds--0.005_g-999999999_start-2501")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,i] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
    }

    auto_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    auto_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/new_fit_calc_modelTest_WeekendRun/allo_DFE_K-100_m-0.05_r-0.693147_u_del-2.5e-11_u_ben-2.5e-12_rho-2.5e-11_bs-0.0_ds--0.005_g-999999999_start-2501")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,i] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
    }

    allo_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    allo_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)


    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/new_fit_calc_modelTest_WeekendRun/diploid_DFE_K-100_m-0.05_r-0.693147_u_del-2.5e-11_u_ben-2.5e-12_rho-2.5e-11_bs-0.0_ds--0.005_g-999999999_start-2501")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]
        edge_deme <- NULL
        for(gen in 1:length(pops[,1])){
            edge_mat[gen,i] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
    }

    diploid_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    diploid_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    x <- as.integer(row.names(dat))
    
    png(file=paste(c("all",file,"_byDeme.png"),sep="", collapse=""), width = 1000, height = 500)
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = 'blue')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = 'blue')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x,diploid_mean,col="black",lwd=5)
    lines(x, diploid_mean+diploid_sd, lty = 'dashed', col = 'black')
    lines(x, diploid_mean-diploid_sd, lty = 'dashed', col = 'black')
    lines(x,allo_mean,col="red",lwd=5)
    lines(x, allo_mean+allo_sd, lty = 'dashed', col = 'red')
    lines(x, allo_mean-allo_sd, lty = 'dashed', col = 'red')
    legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
        col=c("Black", "Blue", "Red"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

}



nrep <- 10
deme_plot <- 2000

file_list <- c("_meanFitnessScaled")
#file_list <- c("_fixedMutations","_meanBenMutationsPerGenome","_meanBenMutationsPerInd","_unscaledFitness","_maxPhs","_heterozygosity","_fixedMutations","_meanDelMutationsPerGenome","_meanDelMutationsPerInd","_meanFitnessScaled")


for(file in file_list){
    print(file)
    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/PPToL_run/Basic/diploid_recessive_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-2501")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,i] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
    }

    auto_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    auto_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/PPToL_run/Diploidization/diploidization-dom_dipLambda-1_remDipMuts-0_recessive_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-2.5e-10_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-2501")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,i] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
    }

    diploidized_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    diploidized_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    x <- as.integer(row.names(dat))
    
    png(file=paste(c("Diploidized_all_",file,"_byDeme.png"),sep="", collapse=""), width = 1000, height = 500)
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploidized_mean-diploidized_sd)),max(c(auto_mean+auto_sd,diploidized_mean+diploidized_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = 'blue')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = 'blue')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x,diploidized_mean,col="black",lwd=5)
    lines(x, diploidized_mean+diploidized_sd, lty = 'dashed', col = 'orange')
    lines(x, diploidized_mean-diploidized_sd, lty = 'dashed', col = 'orange')
    legend(1, 450, legend=c("No diploidization", "Diploidization"),
        col=c("Blue", "Orange"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

}