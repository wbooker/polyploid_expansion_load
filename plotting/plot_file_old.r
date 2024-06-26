######### Plot edge value #######

nrep <- 25
deme_plot <- 500
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
file_list <- c("_meanFitnessScaled")
edge_lab <- "Value"
#file_list <- c("_fixedMutations","_unscaledFitness","_heterozygosity","_fixedMutations","_meanDelMutationsPerGenome","_meanDelMutationsPerInd","_meanFitnessScaled")


for(file in file_list){
    print(file)
    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/dipAllo/allo_DFE_chunk-10000_K-100_m-0.05_r-0.693147_u_del-5.0e-08_u_ben-5.0e-09_rho-5.0e-08_bs-0.0_ds--0.005_g-999999_start-2501")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in c(1:4, 6:10)){
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


    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/basic/diploid_DFE_K-100_m-0.05_r-0.693147_u_del-5.0e-08_u_ben-5.0e-09_rho-5.0e-08_bs-0.0_ds--0.001472_g-999999_start-2501")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in c(1:4, 6:10)){
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

    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/dipAllo/allo_DFE_chunk-1000_K-100_m-0.05_r-0.693147_u_del-5.0e-08_u_ben-5.0e-09_rho-5.0e-08_bs-0.0_ds--0.005_g-999999_start-2501")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in c(1:4, 6:10)){
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

    x <- as.integer(row.names(dat))
    
    png(file=paste(c("all",file,"_byDeme_edge2.png"),sep="", collapse=""), width = 1000, height = 500)
    plot((x),auto_mean, type="l", col="#A1B5D8", lwd=5, xlab="Generation", ylab=paste(c("Edge", edge_lab), collapse = " "), ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = '#A1B5D8')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = '#A1B5D8')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x,diploid_mean,col="black",lwd=5)
    lines(x, diploid_mean+diploid_sd, lty = 'dashed', col = 'black')
    lines(x, diploid_mean-diploid_sd, lty = 'dashed', col = 'black')
    lines(x,allo_mean,col="#DCA92E",lwd=5)
    lines(x, allo_mean+allo_sd, lty = 'dashed', col = '#DCA92E')
    lines(x, allo_mean-allo_sd, lty = 'dashed', col = '#DCA92E')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

    pdf(file=paste(c("all",file,"_byDeme_edge2.pdf"),sep="", collapse=""))
    plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=paste(c("Edge", edge_lab), collapse = " "), ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = '#2B6F59')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = '#2B6F59')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x,diploid_mean,col="black",lwd=5)
    lines(x, diploid_mean+diploid_sd, lty = 'dashed', col = 'black')
    lines(x, diploid_mean-diploid_sd, lty = 'dashed', col = 'black')
    lines(x,allo_mean,col="#A1B5D8",lwd=5)
    lines(x, allo_mean+allo_sd, lty = 'dashed', col = '#A1B5D8')
    lines(x, allo_mean-allo_sd, lty = 'dashed', col = '#A1B5D8')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "2B6F59", "DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

}

####### Core value ###########

for(file in file_list){
    print(file)
    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/basic/auto_DFE_K-100_m-0.05_r-0.693147_u_del-5.0e-08_u_ben-5.0e-09_rho-5.0e-08_bs-0.0_ds--0.001472_g-999999_start-2501")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,i] <- dat[gen,1]
        }
        print(i)
    }

    auto_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    auto_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/basic/allo_DFE_K-100_m-0.05_r-0.693147_u_del-5.0e-08_u_ben-5.0e-09_rho-5.0e-08_bs-0.0_ds--0.001472_g-999999_start-2501")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,i] <- dat[gen,1]
        }
        print(i)
    }

    allo_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    allo_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)


    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/basic/diploid_DFE_K-100_m-0.05_r-0.693147_u_del-5.0e-08_u_ben-5.0e-09_rho-5.0e-08_bs-0.0_ds--0.001472_g-999999_start-2501")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]
        edge_deme <- NULL
        for(gen in 1:length(pops[,1])){
            edge_mat[gen,i] <- dat[gen,1]
        }
        print(i)
    }

    diploid_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    diploid_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    x <- as.integer(row.names(dat))
    
    png(file=paste(c("all",file,"_byDeme_core.png"),sep="", collapse=""), width = 1000, height = 500)
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=paste(c("Core", edge_lab), collapse = " "), ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = '#2B6F59')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = '#2B6F59')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x,diploid_mean,col="black",lwd=5)
    lines(x, diploid_mean+diploid_sd, lty = 'dashed', col = 'black')
    lines(x, diploid_mean-diploid_sd, lty = 'dashed', col = 'black')
    lines(x,allo_mean,col="#DCA92E",lwd=5)
    lines(x, allo_mean+allo_sd, lty = 'dashed', col = '#DCA92E')
    lines(x, allo_mean-allo_sd, lty = 'dashed', col = '#DCA92E')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

    pdf(file=paste(c("all",file,"_byDeme_core.pdf"),sep="", collapse=""))
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=paste(c("Core", edge_lab), collapse = " "), ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = '#2B6F59')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = '#2B6F59')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x,diploid_mean,col="black",lwd=5)
    lines(x, diploid_mean+diploid_sd, lty = 'dashed', col = 'black')
    lines(x, diploid_mean-diploid_sd, lty = 'dashed', col = 'black')
    lines(x,allo_mean,col="#DCA92E",lwd=5)
    lines(x, allo_mean+allo_sd, lty = 'dashed', col = '#DCA92E')
    lines(x, allo_mean-allo_sd, lty = 'dashed', col = '#DCA92E')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "2B6F59", "DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

}

######### Mid value ##########

for(file in file_list){
    print(file)
    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/basic/auto_additive_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-2501")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,i] <- dat[gen,max(which(pops[gen,] > 0)) %/% 2]
        }
        print(i)
    }

    auto_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    auto_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/basic/allo_additive_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-2501")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,i] <- dat[gen,max(which(pops[gen,] > 0)) %/% 2]
        }
        print(i)
    }

    allo_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    allo_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)


    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/basic/diploid_additive_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-2501")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]
        edge_deme <- NULL
        for(gen in 1:length(pops[,1])){
            edge_mat[gen,i] <- dat[gen,max(which(pops[gen,] > 0)) %/% 2]
        }
        print(i)
    }

    diploid_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    diploid_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    x <- as.integer(row.names(dat))
    
    png(file=paste(c("all",file,"_byDeme_mid.png"),sep="", collapse=""), width = 1000, height = 500)
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=paste(c("Mid", edge_lab), collapse = " "), ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = '#2B6F59')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = '#2B6F59')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x,diploid_mean,col="black",lwd=5)
    lines(x, diploid_mean+diploid_sd, lty = 'dashed', col = 'black')
    lines(x, diploid_mean-diploid_sd, lty = 'dashed', col = 'black')
    lines(x,allo_mean,col="#DCA92E",lwd=5)
    lines(x, allo_mean+allo_sd, lty = 'dashed', col = '#DCA92E')
    lines(x, allo_mean-allo_sd, lty = 'dashed', col = '#DCA92E')
    legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
        col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

    pdf(file=paste(c("all",file,"_byDeme_mid.pdf"),sep="", collapse=""))
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=paste(c("Mid", edge_lab), collapse = " "), ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = '#2B6F59')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = '#2B6F59')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x,diploid_mean,col="black",lwd=5)
    lines(x, diploid_mean+diploid_sd, lty = 'dashed', col = 'black')
    lines(x, diploid_mean-diploid_sd, lty = 'dashed', col = 'black')
    lines(x,allo_mean,col="#DCA92E",lwd=5)
    lines(x, allo_mean+allo_sd, lty = 'dashed', col = '#DCA92E')
    lines(x, allo_mean-allo_sd, lty = 'dashed', col = '#DCA92E')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "2B6F59", "DCA92E"), lty=1, cex=0.8)
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


################ diff model prop disomic #################

nrep <- 25
deme_plot <- 1000
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
file_list <- c("_dipIndex")
edge_lab <- "Dip Index"
mod_dir <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/pair_fit_dip_FIXED/test_recomb/recessive", full.names = TRUE, recursive = FALSE)
#file_list <- c("_fixedMutations","_meanBenMutationsPerGenome","_meanBenMutationsPerInd","_unscaledFitness","_maxPhs","_heterozygosity","_fixedMutations","_meanDelMutationsPerGenome","_meanDelMutationsPerInd","_meanFitnessScaled")
for(mod in mod_dir){
    #mod <- "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/pair_fit_dip_FIXED/no_s/diploidization-diff_dipLambda-100_remDipMuts-1_recessive_K-100_m-0.05_r-0.693147_u_del-1.0e-10_u_ben-2.5e-09_u_dip-0.002_rho-1.0e-10_bs-0_ds-0_g-999999_start-2501"
    file <- file_list[1]
    if(!grepl("diploidization-", mod, fixed=TRUE)){
        next
    }
    rep_count <- 0

    setwd(mod)
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]
        if(length(na.omit(pops[,1])) < 100){
            next
        }
        rep_count <- rep_count + 1
        for(gen in 1:length(pops[,1])){
            edge_mat[gen,i] <- dat[gen,(max(which(pops[gen,] > 0)))]
            core_mat[gen,i] <- dat[gen,1]

        }
        print(i)
    }

    edge_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    edge_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    core_mean <- apply(core_mat, 1, mean, na.rm = TRUE)
    core_sd <- apply(core_mat, 1, sd, na.rm = TRUE)

    x <- as.integer(row.names(dat))
    
    png(file=paste(c("all",file,"_byDeme_edge.png"),sep="", collapse=""))
    par( mar=c(5,5,2,2)) 
    plot((x),edge_mean, type="l", col="gray", lwd=5, xlab="", ylab="", ylim = c(0,1), las=1, cex.axis=1.5, cex.lab=1.5)
    mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    mtext(side=2, text="Diploidization Index", line=3.5, cex = 1.5) 
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, edge_mean+edge_sd, lty = 'dashed', col = 'gray')
    lines(x, edge_mean-edge_sd, lty = 'dashed', col = 'gray')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x, core_mean,col="black",lwd=5)
    lines(x, core_mean+core_sd, lty = 'dashed', col = 'black')
    lines(x, core_mean-core_sd, lty = 'dashed', col = 'black')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

    pdf(file=paste(c("all",file,"_byDeme_edge.pdf"),sep="", collapse=""))
    par( mar=c(5,5,2,2)) 
    plot((x),edge_mean, type="l", col="gray", lwd=5, xlab="", ylab="", ylim = c(0,1), las=1, cex.axis=1.5, cex.lab=1.5)
    mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    mtext(side=2, text="Diploidization Index", line=3.5, cex = 1.5) 
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, edge_mean+edge_sd, lty = 'dashed', col = 'gray')
    lines(x, edge_mean-edge_sd, lty = 'dashed', col = 'gray')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x, core_mean,col="black",lwd=5)
    lines(x, core_mean+core_sd, lty = 'dashed', col = 'black')
    lines(x, core_mean-core_sd, lty = 'dashed', col = 'black')
    #    col=c("Black", "2B6F59", "DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

#}

    prop_mat <- matrix(nrow = 100, ncol = rep_count)
    equal_mat <- matrix(nrow = 100, ncol = rep_count)
    prop_1_mat_edge <- matrix(nrow = 100, ncol = rep_count)
    prop_1_mat_core <- matrix(nrow = 100, ncol = rep_count)

    for(i in 1:100){
        for(j in 1:rep_count){
            if(edge_mat[i,j]>=core_mat[i,j]){
                prop_mat[i,j] <- 1
            }
            else{
                prop_mat[i,j] <- 0
            }
            if(edge_mat[i,j]==core_mat[i,j]){
                equal_mat[i,j] <- 1
            }
            if(edge_mat[i,j] >= 0.9){
                prop_1_mat_edge[i,j] = 1
            }
            if(core_mat[i,j] >= 0.9){
                prop_1_mat_core[i,j] = 1
            }
        }
    }
    prop <- apply(prop_mat, 1, sum, na.rm = TRUE)
    equal <- apply(equal_mat, 1, sum, na.rm = TRUE)
    prop_greater <- (prop - equal) / (rep(rep_count, length(equal)) - equal)
    prop_1_edge <- apply(prop_1_mat_edge, 1, sum, na.rm = TRUE) / rep_count
    prop_1_core <- apply(prop_1_mat_core, 1, sum, na.rm = TRUE) / rep_count
    prop_less <- ((rep(rep_count, length(equal)) -prop - equal) / (rep(rep_count, length(equal)) - equal))
    prop_equal <- equal / rep_count

    edge_mean <- apply(edge_mat, 1, var, na.rm = TRUE)
    edge_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    core_mean <- apply(core_mat, 1, var, na.rm = TRUE)
    core_sd <- apply(core_mat, 1, sd, na.rm = TRUE)

    x <- as.integer(row.names(dat))
    
    png(file=paste(c("all",file,"_prop_disomic_byDeme.png"),sep="", collapse=""))
    par( mar=c(5,5,2,2)) 
    plot((x),prop_1_edge, type="l", col="gray", lwd=5, xlab="", ylab="", ylim = c(0,1), las=1, cex.axis=1.5, cex.lab=1.5)
    mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    mtext(side=2, text="Proportion of Simulations Disomic", line=3.5, cex = 1.5)
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #lines(x, edge_mean+edge_sd, lty = 'dashed', col = 'gray')
    #lines(x, edge_mean-edge_sd, lty = 'dashed', col = 'gray')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x, prop_1_core,col="black",lwd=5)
    #lines(x,prop_equal,col='blue',lwd=5)
    #lines(x, core_mean+core_sd, lty = 'dashed', col = 'black')
    #lines(x, core_mean-core_sd, lty = 'dashed', col = 'black')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

    pdf(file=paste(c("all",file,"_prop_disomic_byDeme.pdf"),sep="", collapse=""))
    par( mar=c(5,5,2,2)) 
    plot((x),prop_1_edge, type="l", col="gray", lwd=5, xlab="", ylab="", ylim = c(0,1), las=1, cex.axis=1.5, cex.lab=1.5)
    mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    mtext(side=2, text="Proportion of Simulations Disomic", line=3.5, cex = 1.5)
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #lines(x, edge_mean+edge_sd, lty = 'dashed', col = 'gray')
    #lines(x, edge_mean-edge_sd, lty = 'dashed', col = 'gray')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x, prop_1_core,col="black",lwd=5)
    #lines(x,prop_equal,col='blue',lwd=5)
    #lines(x, core_mean+core_sd, lty = 'dashed', col = 'black')
    #lines(x, core_mean-core_sd, lty = 'dashed', col = 'black')
    #    col=c("Black", "2B6F59", "DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

}

################ dom model prop disomic and dip index #################
nrep <- 25
deme_plot <- 1000
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
file_list <- c("_dipIndex")
edge_lab <- "Dip Index"
mod_dir <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/dip", full.names = TRUE, recursive = FALSE)
#file_list <- c("_fixedMutations","_meanBenMutationsPerGenome","_meanBenMutationsPerInd","_unscaledFitness","_maxPhs","_heterozygosity","_fixedMutations","_meanDelMutationsPerGenome","_meanDelMutationsPerInd","_meanFitnessScaled")
for(mod in mod_dir){
   if(!grepl("diploidization-dom_dipLambda-5", mod, fixed=TRUE)){
        if(!grepl("diploidization-dom_dipLambda-10", mod, fixed=TRUE)){
            next
        }
    }
    nrep <- 25
    #mod <- "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/dip/diploidization-dom_dipLambda-1_remDipMuts-1_DFE_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-2.5e-10_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-2501"
    if(length(list.dirs(path=mod)) >= 100){
        nrep <- 100
    }
    file <- file_list[1]
    setwd(mod)
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]
        if(length(na.omit(pops[,1])) < 100){
            next
        }
        for(gen in 1:length(pops[,1])){
            edge_mat[gen,i] <- dat[gen,(max(which(pops[gen,] > 0)))]
            core_mat[gen,i] <- dat[gen,1]

        }
        print(i)
    }

    edge_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    edge_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    core_mean <- apply(core_mat, 1, mean, na.rm = TRUE)
    core_sd <- apply(core_mat, 1, sd, na.rm = TRUE)

    x <- as.integer(row.names(dat))
    png(file=paste(c("all",file,"_byDeme_edge.png"),sep="", collapse=""))
    par( mar=c(5,5,2,2)) 
    plot((x),edge_mean, type="l", col="gray", lwd=5, xlab="", ylab="", ylim = c(0,1), las=1, cex.axis=1.5, cex.lab=1.5)
    mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    mtext(side=2, text="Diploidization Index", line=3.5, cex = 1.5)
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, edge_mean+edge_sd, lty = 'dashed', col = 'gray')
    lines(x, edge_mean-edge_sd, lty = 'dashed', col = 'gray')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x, core_mean,col="black",lwd=5)
    lines(x, core_mean+core_sd, lty = 'dashed', col = 'black')
    lines(x, core_mean-core_sd, lty = 'dashed', col = 'black')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

    pdf(file=paste(c("all",file,"_byDeme_edge.pdf"),sep="", collapse=""))
    par( mar=c(5,5,2,2)) 
    plot((x),edge_mean, type="l", col="gray", lwd=5, xlab="", ylab="", ylim = c(0,1), las=1, cex.axis=1.5, cex.lab=1.5)
    mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    mtext(side=2, text="Diploidization Index", line=3.5, cex = 1.5)    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, edge_mean+edge_sd, lty = 'dashed', col = 'gray')
    lines(x, edge_mean-edge_sd, lty = 'dashed', col = 'gray')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x, core_mean,col="black",lwd=5)
    lines(x, core_mean+core_sd, lty = 'dashed', col = 'black')
    lines(x, core_mean-core_sd, lty = 'dashed', col = 'black')
    #    col=c("Black", "2B6F59", "DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()
}

nrep <- 25
deme_plot <- 1000
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
file_list <- c("_dipIndex")
edge_lab <- "Dip Index"
mod_dir <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/dip", full.names = TRUE, recursive = FALSE)
#file_list <- c("_fixedMutations","_meanBenMutationsPerGenome","_meanBenMutationsPerInd","_unscaledFitness","_maxPhs","_heterozygosity","_fixedMutations","_meanDelMutationsPerGenome","_meanDelMutationsPerInd","_meanFitnessScaled")
for(mod in mod_dir){
    if(!grepl("diploidization-dom_dipLambda-5", mod, fixed=TRUE)){
        if(!grepl("diploidization-dom_dipLambda-10", mod, fixed=TRUE)){
            next
        }
    }
    nrep <- 25
    if(length(list.dirs(path=mod)) >= 100){
        nrep <- 100
    }
file <- c("_dipIndex")
edge_lab <- "Dip Index"
rep_count <- 0
    #mod <- "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/dip/diploidization-dom_dipLambda-1_remDipMuts-1_DFE_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-2.5e-10_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-2501"
    setwd(mod)
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        if(length(na.omit(pops[,1])) < 100){
            next
        }
        rep_count <- rep_count + 1

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,rep_count] <- dat[gen,max(which(pops[gen,] > 0))]
            core_mat[gen,rep_count] <- dat[gen,1]

        }
        print(i)
    }

    prop_mat <- matrix(nrow = 100, ncol = rep_count)
    equal_mat <- matrix(nrow = 100, ncol = rep_count)
    prop_1_mat_edge <- matrix(nrow = 100, ncol = rep_count)
    prop_1_mat_core <- matrix(nrow = 100, ncol = rep_count)

    for(i in 1:100){
        for(j in 1:rep_count){
            if(edge_mat[i,j]>=core_mat[i,j]){
                prop_mat[i,j] <- 1
            }
            else{
                prop_mat[i,j] <- 0
            }
            if(edge_mat[i,j]==core_mat[i,j]){
                equal_mat[i,j] <- 1
            }
            if(edge_mat[i,j] == 1.0){
                prop_1_mat_edge[i,j] = 1
            }
            if(core_mat[i,j] == 1.0){
                prop_1_mat_core[i,j] = 1
            }
        }
    }
    prop <- apply(prop_mat, 1, sum, na.rm = TRUE)
    equal <- apply(equal_mat, 1, sum, na.rm = TRUE)
    prop_greater <- (prop - equal) / (rep(rep_count, length(equal)) - equal)
    prop_1_edge <- apply(prop_1_mat_edge, 1, sum, na.rm = TRUE) / rep_count
    prop_1_core <- apply(prop_1_mat_core, 1, sum, na.rm = TRUE) / rep_count
    prop_less <- ((rep(rep_count, length(equal)) -prop - equal) / (rep(rep_count, length(equal)) - equal))
    prop_equal <- equal / rep_count

    edge_mean <- apply(edge_mat, 1, var, na.rm = TRUE)
    edge_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    core_mean <- apply(core_mat, 1, var, na.rm = TRUE)
    core_sd <- apply(core_mat, 1, sd, na.rm = TRUE)

    x <- as.integer(row.names(dat))
    op <- par(mar = c(5,7,4,2) + 0.1)

    png(file=paste(c("all",file,"_prop_1_byDeme.png"),sep="", collapse=""))
    par( mar=c(5,5,2,2)) 
    #plot((x),prop_1_edge, type="l", col="gray", lwd=5, xlab="Generation", ylab="Proportion of Simulations Disomic", ylim = c(0,1), las=1, cex.axis=1.5, cex.lab=1.5)
    plot((x),prop_1_edge, type="l", col="gray", lwd=5, xlab="", ylab="", ylim = c(0,1), las=1, cex.axis=1.5, cex.lab=1.5)
    mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    mtext(side=2, text="Proportion of Simulations Disomic", line=3.5, cex = 1.5)
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #lines(x, edge_mean+edge_sd, lty = 'dashed', col = 'gray')
    #lines(x, edge_mean-edge_sd, lty = 'dashed', col = 'gray')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x, prop_1_core,col="black",lwd=5)
    #lines(x,prop_equal,col='blue',lwd=5)
    #lines(x, core_mean+core_sd, lty = 'dashed', col = 'black')
    #lines(x, core_mean-core_sd, lty = 'dashed', col = 'black')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

    pdf(file=paste(c("all",file,"_prop_1_byDeme.pdf"),sep="", collapse=""))
    par( mar=c(5,5,2,2)) 
    plot((x),prop_1_edge, type="l", col="gray", lwd=5, xlab="", ylab="", ylim = c(0,1), las=1, cex.axis=1.5, cex.lab=1.5)
    mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    mtext(side=2, text="Proportion of Simulations Disomic", line=3.5, cex = 1.5)
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #lines(x, edge_mean+edge_sd, lty = 'dashed', col = 'gray')
    #lines(x, edge_mean-edge_sd, lty = 'dashed', col = 'gray')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x, prop_1_core,col="black",lwd=5)
    #lines(x,prop_equal,col='blue',lwd=5)
    #lines(x, core_mean+core_sd, lty = 'dashed', col = 'black')
    #lines(x, core_mean-core_sd, lty = 'dashed', col = 'black')
    #    col=c("Black", "2B6F59", "DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()
}


#################### SFS ###########################

nrep <- 25
deme_plot <- 1000
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
file_list <- c("_coreSiteFrequencies")
edge_lab <- "Core SFS"
mod_dir <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/basic", full.names = TRUE, recursive = FALSE)
#file_list <- c("_fixedMutations","_meanBenMutationsPerGenome","_meanBenMutationsPerInd","_unscaledFitness","_maxPhs","_heterozygosity","_fixedMutations","_meanDelMutationsPerGenome","_meanDelMutationsPerInd","_meanFitnessScaled")
#for(mod in mod_dir){
file <- c("_coreSiteFrequencies")
edge_lab <- "Core SFS"
rep_count <- 0
    mod <- "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/basic/diploid_recessive_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-2501"
    setwd(mod)
    start_mat <- NULL
    end_mat <- NULL
    for(i in 1:nrep){
        curr_file <- paste(c(i,"/",i,file,".csv"),sep="", collapse="")
        n_muts <- count.fields(curr_file, sep = ',')
        n_cols <- max(count.fields(curr_file, sep = ','))
        dat <- read.table(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), header = FALSE, sep = ",", row.names = 1, col.names = paste0("V",seq_len(n_cols)), fill = TRUE)
        tmp_start <- as.vector(dat[1,1:n_muts[1]-1])
        names(tmp_start) <- NULL
        tmp_start <- unlist(tmp_start)

        tmp_end <- as.vector(dat[length(dat[,1]),1:n_muts[length(n_muts)]-1])
        names(tmp_end) <- NULL
        tmp_end <- unlist(tmp_end)

        start_mat <- c(start_mat,tmp_start)
        end_mat <- c(end_mat,tmp_end)
    }
    
    write.table(as.vector(start_mat), "core_start_sfs.csv", col.names=FALSE, row.names=FALSE, sep = ",")
    write.table(as.vector(end_mat), "core_end_sfs.csv", col.names=FALSE, row.names=FALSE, sep = ",")


    png(file=paste(c("all",file,"_prop_1_byDeme.png"),sep="", collapse=""), width = 1000, height = 500)
    plot((x),prop_1_edge, type="l", col="gray", lwd=5, xlab="Generation", ylab="Prop Dip Index > 0.9", ylim = c(0,1))
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #lines(x, edge_mean+edge_sd, lty = 'dashed', col = 'gray')
    #lines(x, edge_mean-edge_sd, lty = 'dashed', col = 'gray')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x, prop_1_core,col="black",lwd=5)
    #lines(x,prop_equal,col='blue',lwd=5)
    #lines(x, core_mean+core_sd, lty = 'dashed', col = 'black')
    #lines(x, core_mean-core_sd, lty = 'dashed', col = 'black')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

    pdf(file=paste(c("all",file,"_prop_greater_byDeme.pdf"),sep="", collapse=""))
    plot((x),prop_1_edge, type="l", col="gray", lwd=5, xlab="Generation", ylab="Prop Dip Index > 0.9", ylim = c(0,1))
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #lines(x, edge_mean+edge_sd, lty = 'dashed', col = 'gray')
    #lines(x, edge_mean-edge_sd, lty = 'dashed', col = 'gray')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x, prop_1_core,col="black",lwd=5)
    #lines(x,prop_equal,col='blue',lwd=5)
    #lines(x, core_mean+core_sd, lty = 'dashed', col = 'black')
    #lines(x, core_mean-core_sd, lty = 'dashed', col = 'black')
    #    col=c("Black", "2B6F59", "DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()
