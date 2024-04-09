
######### Plot edge value #######

nrep <- 100
deme_plot <- 2010
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
#file_list <- c("_unscaledFitness")
#edge_lab <- c("Edge Mean Fitness (unscaled)")
file_list <- c("_meanFitnessScaled", "_meanDelMutationsPerGenome", "_meanDelMutationsPerInd", "_meanBenMutationsPerGenome", "_meanBenMutationsPerInd", "_fixedMutations")
edge_lab <- c("Edge Mean Fitness","Edge Mean Deleterious Mutations Per Genome", "Edge Mean Deleterious Mutations Per Individual", "Edge Mean Beneficial Mutations Per Genome", "Edge Mean Beneficial Mutations Per Individual", "Edge Fixed Deleterious Mutations")


for(my_file in 1:length(file_list)){
    file <- file_list[my_file]
    print(file)
    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/basic/auto_additive_fixed_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.005_ds--0.005_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
        n_ <- n_ + 1
    }

    auto_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    auto_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)


    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/basic/diploid_additive_fixed_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.005_ds--0.005_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
        n_ <- n_ + 1
    }

    diploid_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    diploid_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/basic/allo_additive_fixed_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.005_ds--0.005_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
        n_ <- n_ + 1
    }

    allo_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    allo_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    allo_mean <- allo_mean[!is.na(allo_mean)]
    auto_mean <- auto_mean[!is.na(auto_mean)]
    diploid_mean <- diploid_mean[!is.na(diploid_mean)]
    allo_sd <- allo_sd[!is.na(allo_sd)]
    auto_sd <- auto_sd[!is.na(auto_sd)]
    diploid_sd <- diploid_sd[!is.na(diploid_sd)]
    x <- seq(0,10*(length(auto_sd)-1),10)

    
    png(file=paste(c("all",file,"_byDeme_edge2.png"),sep="", collapse=""), width = 1000, height = 1000)
    plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab[my_file], ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
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

    pdf(file=paste(c("all",file,"_byDeme_edge2.pdf"),sep="", collapse=""))
    plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab[my_file], ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
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

##############################################################
nrep <- 100
deme_plot <- 2010
for(my_file in 1:length(file_list)){
    file <- file_list[my_file]
    print(file)
    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/basic/auto_DFE_gamma_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.001472_ds--0.001472_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
        n_ <- n_ + 1
    }

    auto_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    auto_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)


    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/basic/diploid_DFE_gamma_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.001472_ds--0.001472_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
        n_ <- n_ + 1
    }

    diploid_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    diploid_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/basic/allo_DFE_gamma_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.001472_ds--0.001472_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
        n_ <- n_ + 1
    }

    allo_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    allo_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    allo_mean <- allo_mean[!is.na(allo_mean)]
    auto_mean <- auto_mean[!is.na(auto_mean)]
    diploid_mean <- diploid_mean[!is.na(diploid_mean)]
    allo_sd <- allo_sd[!is.na(allo_sd)]
    auto_sd <- auto_sd[!is.na(auto_sd)]
    diploid_sd <- diploid_sd[!is.na(diploid_sd)]
    x <- seq(0,10*(length(auto_sd)-1),10)

    
    png(file=paste(c("all",file,"_byDeme_edge2.png"),sep="", collapse=""), width = 1000, height = 1000)
    plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab[my_file], ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
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

    pdf(file=paste(c("all",file,"_byDeme_edge2.pdf"),sep="", collapse=""))
    plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab[my_file], ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = '#2B6F59')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = '#2B6F59')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x,diploid_mean,col="black",lwd=5)
    lines(x, diploid_mean+diploid_sd, lty = 'dashed', col = 'black')
    lines(x, diploid_mean-diploid_sd, lty = 'dashed', col = 'black')
    lines(x,allo_mean,col="#DCA92E",lwd=5)
    lines(x, allo_mean+allo_sd, lty = 'dashed', col = '#DCA92E')
    lines(x, allo_mean-allo_sd, lty = 'dashed', col = '#DCA92E')
    #    col=c("Black", "2B6F59", "DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

}

nrep <- 100
deme_plot <- 2010
file_list <- c("_meanFitnessScaled", "_meanDelMutationsPerGenome", "_meanDelMutationsPerInd", "_meanBenMutationsPerGenome", "_meanBenMutationsPerInd", "_fixedMutations", "_unscaledFitness")
edge_lab <- c("Edge Mean Fitness","Edge Mean Deleterious Mutations Per Genome", "Edge Mean Deleterious Mutations Per Individual", "Edge Mean Beneficial Mutations Per Genome", "Edge Mean Beneficial Mutations Per Individual", "Edge Fixed Deleterious Mutations", "Edge Mean Fitness (Unscaled)")

for(my_file in 1:length(file_list)){
    file <- file_list[my_file]
    print(file)
    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/basic/auto_additive_gamma_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.001472_ds--0.001472_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        if(length(dat[,1]) < deme_plot/10 - 1){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
        n_ <- n_ + 1
    }

    auto_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    auto_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)


    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/basic/diploid_additive_gamma_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.001472_ds--0.001472_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        if(length(dat[,1]) < deme_plot/10 - 1){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
        n_ <- n_ + 1
    }

    diploid_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    diploid_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/basic/allo_additive_gamma_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.001472_ds--0.001472_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        if(length(dat[,1]) < deme_plot/10 - 1){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
        n_ <- n_ + 1
    }

    allo_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    allo_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    #x <- as.integer(row.names(dat))[1:201]
    #x <- x[!is.na(x)]
    allo_mean <- allo_mean[!is.na(allo_mean)]
    auto_mean <- auto_mean[!is.na(auto_mean)]
    diploid_mean <- diploid_mean[!is.na(diploid_mean)]
    allo_sd <- allo_sd[!is.na(allo_sd)]
    auto_sd <- auto_sd[!is.na(auto_sd)]
    diploid_sd <- diploid_sd[!is.na(diploid_sd)]
    x <- seq(0,10*(length(auto_sd)-1),10)

    
    png(file=paste(c("all",file,"_byDeme_edge2.png"),sep="", collapse=""), width = 1000, height = 1000)
    plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab[my_file], ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
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

    pdf(file=paste(c("all",file,"_byDeme_edge2.pdf"),sep="", collapse=""))
    plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab[my_file], ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
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

nrep <- 100
deme_plot <- 2010
file_list <- c("_meanBenMutationsPerInd", "_fixedMutations")
edge_lab <- c("Edge Mean Beneficial Mutations Per Individual", "Edge Fixed Deleterious Mutations")


for(my_file in 1:length(file_list)){
    file <- file_list[my_file]
    print(file)
    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/basic/auto_bd_dr_gamma_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.001472_ds--0.001472_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        if(length(dat[,1]) < deme_plot/10 - 1){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
        n_ <- n_ + 1
    }

    auto_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    auto_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)


    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/basic/diploid_bd_dr_gamma_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.001472_ds--0.001472_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        if(length(dat[,1]) < deme_plot/10 - 1){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
        n_ <- n_ + 1
    }

    diploid_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    diploid_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/test_burnin_2/auto_duplex_gamma_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.001472_ds--0.001472_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        if(length(dat[,1]) < deme_plot/10 - 1){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,max(which(pops[gen,] > 0))]
        }
        print(i)
        n_ <- n_ + 1
    }

    allo_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    allo_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    #x <- as.integer(row.names(dat))[1:201]
    #x <- x[!is.na(x)]
    allo_mean <- allo_mean[!is.na(allo_mean)]
    auto_mean <- auto_mean[!is.na(auto_mean)]
    diploid_mean <- diploid_mean[!is.na(diploid_mean)]
    allo_sd <- allo_sd[!is.na(allo_sd)]
    auto_sd <- auto_sd[!is.na(auto_sd)]
    allo_mean <- allo_mean[1:length(auto_mean)]
    allo_sd <- allo_sd[1:length(auto_mean)]
    diploid_sd <- diploid_sd[!is.na(diploid_sd)]
    x <- seq(0,10*(length(auto_sd)-1),10)

    
    png(file=paste(c("all",file,"_byDeme_edge2.png"),sep="", collapse=""), width = 1000, height = 1000)
    plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab[my_file], ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
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
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)   

    dev.off()

    pdf(file=paste(c("all",file,"_byDeme_edge2.pdf"),sep="", collapse=""))
    plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab[my_file], ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
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



################### cores ###############
####### Core value ###########
nrep <- 100
deme_plot <- 2010
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

file_list <- c("_meanFitnessScaled")
edge_lab <- c("Mean Fitness")

for(file in file_list){
    print(file)
    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/basic/auto_bd_dr_gamma_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.001472_ds--0.001472_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        if(length(dat[,1]) < deme_plot/10 - 1){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,1]
        }
        print(i)
        n_ <- n_ + 1
    }

    auto_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    auto_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/basic/diploid_bd_dr_gamma_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.001472_ds--0.001472_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        if(length(dat[,1]) < deme_plot/10 - 1){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,1]
        }
        print(i)
        n_ <- n_ + 1
    }

    diploid_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    diploid_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/test_burnin_2/auto_duplex_gamma_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.001472_ds--0.001472_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:103)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        if(length(dat[,1]) < deme_plot/10 - 1){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,1]
        }
        print(i)
        n_ <- n_ + 1
    }

    allo_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    allo_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)




    allo_mean <- allo_mean[!is.na(allo_mean)]
    auto_mean <- auto_mean[!is.na(auto_mean)]
    diploid_mean <- diploid_mean[!is.na(diploid_mean)]
    allo_sd <- allo_sd[!is.na(allo_sd)]
    auto_sd <- auto_sd[!is.na(auto_sd)]
    diploid_sd <- diploid_sd[!is.na(diploid_sd)]
    x <- seq(0,10*(length(auto_sd)-1),10)
    
    png(file=paste(c("all",file,"_byDeme_core.png"),sep="", collapse=""), width = 1000, height = 500)
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=paste(c("Core", edge_lab), collapse = " "), ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
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
    lines(x,allo_mean,col="#A1B5D8",lwd=5)
    lines(x, allo_mean+allo_sd, lty = 'dashed', col = '#A1B5D8')
    lines(x, allo_mean-allo_sd, lty = 'dashed', col = '#A1B5D8')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "2B6F59", "DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

}

################### pi burnin ###############
nrep <- 100
deme_plot <- 15010
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

file <- c("_heterozygosity")
edge_lab <- c("Nucleotide Diversity")

    print(file)
    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/test_burnin_2/diploid_DFE_gamma_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.001472_ds--0.001472_g-999999_start-15000")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    n_ <- 1
    for(i in c(1:20)){
        if(n_ > nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        if(length(dat[,1]) < deme_plot/10 - 1){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,n_] <- dat[gen,1]
        }
        print(i)
        n_ <- n_ + 1
    }

    auto_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    auto_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)
    auto_mean <- auto_mean[!is.na(auto_mean)]
    auto_sd <- auto_sd[!is.na(auto_sd)]

    x <- seq(0,10*(length(auto_sd)-1),10)
    
    png(file=paste(c("burnin",file,"_byDeme_core.png"),sep="", collapse=""), width = 1000, height = 500)
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    plot((x),auto_mean, type="l", col="black", lwd=5, xlab="Generation", ylab=paste(c("Core", edge_lab), collapse = " "), ylim = c(min(c(auto_mean-auto_sd)),max(c(auto_mean+auto_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = 'black')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = 'black')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

    pdf(file=paste(c("burnin",file,"_byDeme_core.pdf"),sep="", collapse=""))
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    plot((x),auto_mean, type="l", col="black", lwd=5, xlab="Generation", ylab=paste(c("Core", edge_lab), collapse = " "), ylim = c(min(c(auto_mean-auto_sd)),max(c(auto_mean+auto_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = 'black')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = 'black')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

################ diff model prop disomic and mean dip by s #################

nrep <- 100
deme_plot <- 1000
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
file_list <- c("_dipIndex")
edge_lab <- "Dip Index"
#mod_dir <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/pair_fit_dip_FIXED", full.names = TRUE, recursive = FALSE)
#file_list <- c("_fixedMutations","_meanBenMutationsPerGenome","_meanBenMutationsPerInd","_unscaledFitness","_maxPhs","_heterozygosity","_fixedMutations","_meanDelMutationsPerGenome","_meanDelMutationsPerInd","_meanFitnessScaled")
#for(mod in mod_dir){
    mod <- "/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/dip/diploidization-dom_dipLambda-2_remDipMuts-1_bd_drexp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-5.0e-10_rho-1.0e-08_bs-0.01_ds--0.01_g-999999_start-5001"
    file <- file_list[1]
    rep_count <- 0
    #mod <- "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/pair_fit_dip_FIXED/diploidization-diff_dipLambda-100_remDipMuts-1_recessive_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-2.5e-08_bs-0.0_ds--0.01_g-999999_start-2501"
    setwd(mod)
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:120){
        if(rep_count == nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        if(length(na.omit(pops[,1])) < 100){
            next
        }
        rep_count <- rep_count + 1

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,rep_count] <- dat[gen,(max(which(pops[gen,] > 0)))]
            core_mat[gen,rep_count] <- dat[gen,1]

        }
        print(i)
    }
    print(rep_count)
    edge_mean_1e2 <- apply(edge_mat, 1, mean, na.rm = TRUE)
    edge_sd_1e2 <- apply(edge_mat, 1, sd, na.rm = TRUE)

    core_mean_1e2 <- apply(core_mat, 1, mean, na.rm = TRUE)
    core_sd_1e2 <- apply(core_mat, 1, sd, na.rm = TRUE)

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
    prop_1_edge_1e2 <- apply(prop_1_mat_edge, 1, sum, na.rm = TRUE) / rep_count
    prop_1_core_1e2 <- apply(prop_1_mat_core, 1, sum, na.rm = TRUE) / rep_count

    mod <- "/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/dip/diploidization-dom_dipLambda-2_remDipMuts-1_bd_drexp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-5.0e-10_rho-1.0e-08_bs-0.005_ds--0.005_g-999999_start-5001"
    file <- file_list[1]
    rep_count <- 0
    #mod <- "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/pair_fit_dip_FIXED/diploidization-diff_dipLambda-100_remDipMuts-1_recessive_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-2.5e-08_bs-0.0_ds--0.01_g-999999_start-2501"
    setwd(mod)
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:120){
        if(rep_count == nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        if(length(na.omit(pops[,1])) < 100){
            next
        }
        rep_count <- rep_count + 1

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,rep_count] <- dat[gen,(max(which(pops[gen,] > 0)))]
            core_mat[gen,rep_count] <- dat[gen,1]

        }
        print(i)
    }
    print(rep_count)
    edge_mean_5e3 <- apply(edge_mat, 1, mean, na.rm = TRUE)
    edge_sd_5e3 <- apply(edge_mat, 1, sd, na.rm = TRUE)

    core_mean_5e3 <- apply(core_mat, 1, mean, na.rm = TRUE)
    core_sd_5e3 <- apply(core_mat, 1, sd, na.rm = TRUE)

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
    prop_1_edge_5e3 <- apply(prop_1_mat_edge, 1, sum, na.rm = TRUE) / rep_count
    prop_1_core_5e3 <- apply(prop_1_mat_core, 1, sum, na.rm = TRUE) / rep_count

    mod <- "/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/dip/diploidization-dom_dipLambda-2_remDipMuts-1_bd_drexp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-5.0e-10_rho-1.0e-08_bs-0.001_ds--0.001_g-999999_start-5001"
    file <- file_list[1]
    rep_count <- 0
    #mod <- "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/pair_fit_dip_FIXED/diploidization-diff_dipLambda-100_remDipMuts-1_recessive_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-2.5e-08_bs-0.0_ds--0.01_g-999999_start-2501"
    setwd(mod)
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in c(1:120)){
        if(rep_count == nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        if(length(na.omit(pops[,1])) < 100){
            next
        }
        rep_count <- rep_count + 1

        for(gen in 1:length(pops[,1])){
            edge_mat[gen,rep_count] <- dat[gen,(max(which(pops[gen,] > 0)))]
            core_mat[gen,rep_count] <- dat[gen,1]

        }
        print(i)
    }
    print(rep_count)
    edge_mean_1e3 <- apply(edge_mat, 1, mean, na.rm = TRUE)
    edge_sd_1e3 <- apply(edge_mat, 1, sd, na.rm = TRUE)

    core_mean_1e3 <- apply(core_mat, 1, mean, na.rm = TRUE)
    core_sd_1e3 <- apply(core_mat, 1, sd, na.rm = TRUE)
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
    prop_1_edge_1e3 <- apply(prop_1_mat_edge, 1, sum, na.rm = TRUE) / rep_count
    prop_1_core_1e3 <- apply(prop_1_mat_core, 1, sum, na.rm = TRUE) / rep_count

    x <- as.integer(row.names(dat))
    
    png(file=paste(c("all",file,"_prop_disomic_byDeme_by_s_n100.png"),sep="", collapse=""), width = 1000, height = 500)
    plot((x),prop_1_edge_1e3, type="l", col="gray", lwd=5, xlab="Generation", ylab="Proportion of Simulations Disomic", ylim = c(0,1))
    lines(x, prop_1_edge_5e3,col="#87CCB2",lwd=5)
    lines(x, prop_1_edge_1e2,col="#2B6F59",lwd=5)
    lines(x, prop_1_core_1e3,col="gray", lwd=5, lty = 'twodash')
    lines(x, prop_1_core_5e3,col="#87CCB2", lwd=5, lty = 'twodash')
    lines(x, prop_1_core_1e2,col="#2B6F59", lwd=5, lty = 'twodash')

    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #lines(x, edge_mean+edge_sd, lty = 'dashed', col = 'gray')
    #lines(x, edge_mean-edge_sd, lty = 'dashed', col = 'gray')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    #lines(x,prop_equal,col='blue',lwd=5)
    #lines(x, core_mean+core_sd, lty = 'dashed', col = 'black')
    #lines(x, core_mean-core_sd, lty = 'dashed', col = 'black')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

    pdf(file=paste(c("all",file,"_prop_disomic_byDeme_by_s_n100.pdf"),sep="", collapse=""))
    plot((x),prop_1_edge_1e3, type="l", col="gray", lwd=5, xlab="Generation", ylab="Proportion of Simulations Disomic", ylim = c(0,1))
    lines(x, prop_1_edge_5e3,col="#87CCB2",lwd=5)
    lines(x, prop_1_edge_1e2,col="#2B6F59",lwd=5)
    lines(x, prop_1_core_1e3,col="gray", lwd=5, lty = 'twodash')
    lines(x, prop_1_core_5e3,col="#87CCB2", lwd=5, lty = 'twodash')
    lines(x, prop_1_core_1e2,col="#2B6F59", lwd=5, lty = 'twodash')

    dev.off()

    x <- seq(0,990,10)
    
    png(file=paste(c("all",file,"_mean_by_s_n100.png"),sep="", collapse=""), width = 1000, height = 500)
    plot((x),edge_mean_1e3, type="l", col="gray", lwd=5, xlab="Generation", ylab="Mean Diploidization Index", ylim = c(0,1))
    lines(x, edge_mean_5e3,col="#87CCB2",lwd=5)
    lines(x, edge_mean_1e2,col="#2B6F59",lwd=5)
    lines(x, core_mean_1e3,col="gray", lwd=5, lty = 'twodash')
    lines(x, core_mean_5e3,col="#87CCB2", lwd=5, lty = 'twodash')
    lines(x, core_mean_1e2,col="#2B6F59", lwd=5, lty = 'twodash')
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #lines(x, edge_mean+edge_sd, lty = 'dashed', col = 'gray')
    #lines(x, edge_mean-edge_sd, lty = 'dashed', col = 'gray')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    #lines(x,prop_equal,col='blue',lwd=5)
    #lines(x, core_mean+core_sd, lty = 'dashed', col = 'black')
    #lines(x, core_mean-core_sd, lty = 'dashed', col = 'black')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

    pdf(file=paste(c("all",file,"_mean_by_s_n100.pdf"),sep="", collapse=""))
    plot((x),edge_mean_1e3, type="l", col="gray", lwd=5, xlab="Generation", ylab="Mean Diploidization Index", ylim = c(0,1))
    lines(x, edge_mean_5e3,col="#87CCB2",lwd=5)
    lines(x, edge_mean_1e2,col="#2B6F59",lwd=5)
    lines(x, core_mean_1e3,col="gray", lwd=5, lty = 'twodash')
    lines(x, core_mean_5e3,col="#87CCB2", lwd=5, lty = 'twodash')
    lines(x, core_mean_1e2,col="#2B6F59", lwd=5, lty = 'twodash')


    dev.off()
 

#### Prop disomic by rho ####

nrep <- 100
deme_plot <- 1000
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
file_list <- c("_dipIndex")
edge_lab <- "Dip Index"
#mod_dir <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/pair_fit_dip_FIXED", full.names = TRUE, recursive = FALSE)
#file_list <- c("_fixedMutations","_meanBenMutationsPerGenome","_meanBenMutationsPerInd","_unscaledFitness","_maxPhs","_heterozygosity","_fixedMutations","_meanDelMutationsPerGenome","_meanDelMutationsPerInd","_meanFitnessScaled")
#for(mod in mod_dir){
    mod <- "/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/dip/diploidization-diff_dipLambda-100_remDipMuts-1_meioticfFitness-1_pe_inflection-85_pe_slope-1_bd_dr_exp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-1.0e-06_bs-0.01_ds--0.01_g-999999_start-5001"
    file <- file_list[1]
    rep_count <- 0
    #mod <- "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/pair_fit_dip_FIXED/diploidization-diff_dipLambda-100_remDipMuts-1_recessive_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-2.5e-08_bs-0.0_ds--0.01_g-999999_start-2501"
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
            edge_mat[gen,rep_count] <- dat[gen,(max(which(pops[gen,] > 0)))]
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
            if(edge_mat[i,j] >= 0.9){
                prop_1_mat_edge[i,j] = 1
            }
            if(core_mat[i,j] >= 0.9){
                prop_1_mat_core[i,j] = 1
            }
        }
    }
    prop <- apply(prop_mat, 1, sum, na.rm = TRUE)
    prop_1_edge_1e6 <- apply(prop_1_mat_edge, 1, sum, na.rm = TRUE) / rep_count
    prop_1_core_1e6 <- apply(prop_1_mat_core, 1, sum, na.rm = TRUE) / rep_count

    mod <- "/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/dip/diploidization-diff_dipLambda-100_remDipMuts-1_meioticfFitness-1_pe_inflection-85_pe_slope-1_bd_dr_exp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-1.0e-07_bs-0.01_ds--0.01_g-999999_start-5001"
    file <- file_list[1]
    rep_count <- 0
    #mod <- "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/pair_fit_dip_FIXED/diploidization-diff_dipLambda-100_remDipMuts-1_recessive_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-2.5e-08_bs-0.0_ds--0.01_g-999999_start-2501"
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
            edge_mat[gen,rep_count] <- dat[gen,(max(which(pops[gen,] > 0)))]
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
            if(edge_mat[i,j] >= 0.9){
                prop_1_mat_edge[i,j] = 1
            }
            if(core_mat[i,j] >= 0.9){
                prop_1_mat_core[i,j] = 1
            }
        }
    }
    prop <- apply(prop_mat, 1, sum, na.rm = TRUE)
    prop_1_edge_1e7 <- apply(prop_1_mat_edge, 1, sum, na.rm = TRUE) / rep_count
    prop_1_core_1e7 <- apply(prop_1_mat_core, 1, sum, na.rm = TRUE) / rep_count

    mod <- "/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/dip/diploidization-diff_dipLambda-100_remDipMuts-1_meioticfFitness-1_pe_inflection-85_pe_slope-1_bd_dr_exp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-1.0e-08_bs-0.01_ds--0.01_g-999999_start-5001"
    file <- file_list[1]
    rep_count <- 0
    #mod <- "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/pair_fit_dip_FIXED/diploidization-diff_dipLambda-100_remDipMuts-1_recessive_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-2.5e-08_bs-0.0_ds--0.01_g-999999_start-2501"
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
            edge_mat[gen,rep_count] <- dat[gen,(max(which(pops[gen,] > 0)))]
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
            if(edge_mat[i,j] >= 0.9){
                prop_1_mat_edge[i,j] = 1
            }
            if(core_mat[i,j] >= 0.9){
                prop_1_mat_core[i,j] = 1
            }
        }
    }
    prop <- apply(prop_mat, 1, sum, na.rm = TRUE)
    prop_1_edge_1e8 <- apply(prop_1_mat_edge, 1, sum, na.rm = TRUE) / rep_count
    prop_1_core_1e8 <- apply(prop_1_mat_core, 1, sum, na.rm = TRUE) / rep_count

    x <- as.integer(row.names(dat))
    
    png(file=paste(c("all",file,"_prop_disomic_byDeme_by_rho.png"),sep="", collapse=""), width = 1000, height = 500)
    plot((x),prop_1_edge_1e6, type="l", col="gray", lwd=5, xlab="Generation", ylab="Proportion of Simulations Disomic", ylim = c(0,1))
    lines(x, prop_1_edge_1e7,col="#87CCB2",lwd=5)
    lines(x, prop_1_edge_1e8,col="#2B6F59",lwd=5)
    lines(x, prop_1_core_1e6,col="gray", lwd=5, lty = 'twodash')
    lines(x, prop_1_core_1e7,col="#87CCB2", lwd=5, lty = 'twodash')
    lines(x, prop_1_core_1e8,col="#2B6F59", lwd=5, lty = 'twodash')

    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #lines(x, edge_mean+edge_sd, lty = 'dashed', col = 'gray')
    #lines(x, edge_mean-edge_sd, lty = 'dashed', col = 'gray')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    #lines(x,prop_equal,col='blue',lwd=5)
    #lines(x, core_mean+core_sd, lty = 'dashed', col = 'black')
    #lines(x, core_mean-core_sd, lty = 'dashed', col = 'black')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

    pdf(file=paste(c("all",file,"_prop_disomic_byDeme_by_rho.pdf"),sep="", collapse=""))
    plot((x),prop_1_edge_1e6, type="l", col="gray", lwd=3, xlab="Generation", ylab="Proportion of Simulations Disomic", ylim = c(0,1))
    lines(x, prop_1_edge_1e7,col="#87CCB2",lwd=3)
    lines(x, prop_1_edge_1e8,col="#2B6F59",lwd=3)
    lines(x, prop_1_core_1e6,col="gray", lwd=3, lty = 'twodash')
    lines(x, prop_1_core_1e7,col="#87CCB2", lwd=3, lty = 'twodash')
    lines(x, prop_1_core_1e8,col="#2B6F59", lwd=3, lty = 'twodash')

    dev.off()


#### avg dip by s ####

nrep <- 50
deme_plot <- 1000
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
file_list <- c("_dipIndex")
edge_lab <- "Dip Index"

    mod <- "/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/dip/diploidization-diff_dipLambda-100_remDipMuts-1_meioticfFitness-1_pe_inflection-70_pe_slope-1_bd_dr_exp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-1.0e-07_bs-0.01_ds--0.01_g-999999_start-5001"
    file <- file_list[1]
    setwd(mod)
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    rep_count <- 0
    for(i in 1:100){
        if(rep_count == nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
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
        rep_count <- rep_count + 1
        print(i)
    }

    edge_mean_5e3 <- apply(edge_mat, 1, mean, na.rm = TRUE)
    edge_sd_5e3 <- apply(edge_mat, 1, sd, na.rm = TRUE)

    core_mean_5e3 <- apply(core_mat, 1, mean, na.rm = TRUE)
    core_sd_5e3 <- apply(core_mat, 1, sd, na.rm = TRUE)

    mod <- "/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/dip/diploidization-diff_dipLambda-100_remDipMuts-1_meioticfFitness-1_pe_inflection-70_pe_slope-1_bd_dr_exp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-1.0e-07_bs-0.005_ds--0.005_g-999999_start-5001"
    file <- file_list[1]
    setwd(mod)
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    rep_count <- 0
    for(i in 1:100){
        if(rep_count == nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
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
        rep_count <- rep_count + 1
        print(i)
    }

    edge_mean_1e2 <- apply(edge_mat, 1, mean, na.rm = TRUE)
    edge_sd_1e2 <- apply(edge_mat, 1, sd, na.rm = TRUE)

    core_mean_1e2 <- apply(core_mat, 1, mean, na.rm = TRUE)
    core_sd_1e2 <- apply(core_mat, 1, sd, na.rm = TRUE)

    mod <- "/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/dip/diploidization-diff_dipLambda-100_remDipMuts-1_meioticfFitness-1_pe_inflection-70_pe_slope-1_bd_dr_exp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-1.0e-07_bs-0.001_ds--0.001_g-999999_start-5001"
    file <- file_list[1]
    setwd(mod)
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    rep_count <- 0
    for(i in 1:100){
        if(rep_count == nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
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
        rep_count <- rep_count + 1
        print(i)
    }

    edge_mean_1e3 <- apply(edge_mat, 1, mean, na.rm = TRUE)
    edge_sd_1e3 <- apply(edge_mat, 1, sd, na.rm = TRUE)

    core_mean_1e3 <- apply(core_mat, 1, mean, na.rm = TRUE)
    core_sd_1e3 <- apply(core_mat, 1, sd, na.rm = TRUE)

    x <- seq(0,990,10)
    
    png(file=paste(c("all",file,"_mean_by_s.png"),sep="", collapse=""), width = 1000, height = 500)
    plot((x),edge_mean_1e3, type="l", col="gray", lwd=5, xlab="Generation", ylab="Mean Diploidization Index", ylim = c(0,1))
    lines(x, edge_mean_5e3,col="#87CCB2",lwd=5)
    lines(x, edge_mean_1e2,col="#2B6F59",lwd=5)
    lines(x, core_mean_1e3,col="gray", lwd=5, lty = 'twodash')
    lines(x, core_mean_5e3,col="#87CCB2", lwd=5, lty = 'twodash')
    lines(x, core_mean_1e2,col="#2B6F59", lwd=5, lty = 'twodash')
    #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    #lines(x, edge_mean+edge_sd, lty = 'dashed', col = 'gray')
    #lines(x, edge_mean-edge_sd, lty = 'dashed', col = 'gray')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    #lines(x,prop_equal,col='blue',lwd=5)
    #lines(x, core_mean+core_sd, lty = 'dashed', col = 'black')
    #lines(x, core_mean-core_sd, lty = 'dashed', col = 'black')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

    pdf(file=paste(c("all",file,"_mean_by_s.pdf"),sep="", collapse=""))
    plot((x),edge_mean_1e3, type="l", col="gray", lwd=5, xlab="Generation", ylab="Mean Diploidization Index", ylim = c(0,1))
    lines(x, edge_mean_5e3,col="#87CCB2",lwd=5)
    lines(x, edge_mean_1e2,col="#2B6F59",lwd=5)
    lines(x, core_mean_1e3,col="gray", lwd=5, lty = 'twodash')
    lines(x, core_mean_5e3,col="#87CCB2", lwd=5, lty = 'twodash')
    lines(x, core_mean_1e2,col="#2B6F59", lwd=5, lty = 'twodash')


    dev.off()
 


####### plot edge non vs. dip fitness ########    
nrep <- 49
deme_plot <- 1000
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
file_list <- c("_heterozygosity")
edge_lab <- "Value"
#file_list <- c("_fixedMutations","_unscaledFitness","_heterozygosity","_fixedMutations","_meanDelMutationsPerGenome","_meanDelMutationsPerInd","_meanFitnessScaled")


for(file in file_list){
    print(file)
    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/no_PAFC/auto_bd_dr_exp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.005_ds--0.005_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    rep_count <- 0
    for(i in c(1:100)){
        if(rep_count == nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
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
        }
        print(i)
    }

    auto_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    auto_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)


    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/dip/diploidization-diff_dipLambda-100_remDipMuts-1_meioticfFitness-1_pe_inflection-85_pe_slope-1_bd_dr_exp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-1.0e-06_bs-0.005_ds--0.005_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    rep_count <- 0
    for(i in c(1:100)){
        if(rep_count == nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
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
        }
        print(i)
    }

    diploidized_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    diploidized_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/no_PAFC/diploidization-diff_dipLambda-100_remDipMuts-1_meioticfFitness-0_pe_inflection-85_pe_slope-1_bd_dr_exp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-1.0e-06_bs-0.005_ds--0.005_g-999999_start-5001")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    rep_count <- 0
    for(i in c(1:100)){
        if(rep_count == nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
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
        }
        print(i)
    }
    noPAFC_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    noPAFC_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    x <- as.integer(row.names(dat))
    
    png(file=paste(c("dip_not_fitness",file,"_byDeme_edge.png"),sep="", collapse=""))
    plot((x),auto_mean, type="l", col="black", lwd=5, xlab="Generation", ylab="Heterozygosity", las=1, cex.axis=1.5, cex.lab=1.5, ylim = c(min(c(auto_mean-auto_sd,diploidized_mean-diploidized_sd,noPAFC_mean-noPAFC_sd)),max(c(auto_mean+auto_sd,diploidized_mean+diploidized_sd,noPAFC_mean+noPAFC_sd))))
    par( mar=c(5,5,2,2), omi=c(0,0,0,0)) 
    #mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    #mtext(side=2, text="Heterozygosity", line=3.5, cex = 1.5)#plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = 'black')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = 'black')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x,diploidized_mean,col="#A1B5D8",lwd=5)
    lines(x, diploidized_mean+diploidized_sd, lty = 'dashed', col = '#A1B5D8')
    lines(x, diploidized_mean-diploidized_sd, lty = 'dashed', col = '#A1B5D8')
    lines(x,noPAFC_mean,col="#2B6F59",lwd=5)
    lines(x, noPAFC_mean+noPAFC_sd, lty = 'dashed', col = '#2B6F59')
    lines(x, noPAFC_mean-noPAFC_sd, lty = 'dashed', col = '#2B6F59')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()


    pdf(file=paste(c("dip_not_fitness",file,"_byDeme_edge.pdf"),sep="", collapse=""))
    par( mar=c(5,5,2,2)) 
    plot((x),auto_mean, type="l", col="black", lwd=5, xlab="", ylab="", las=1, cex.axis=1.5, cex.lab=1.5, ylim = c(min(c(auto_mean-auto_sd,diploidized_mean-diploidized_sd,noPAFC_mean-noPAFC_sd)),max(c(auto_mean+auto_sd,diploidized_mean+diploidized_sd,noPAFC_mean+noPAFC_sd))))
    mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    mtext(side=2, text="Heterozygosity", line=3.5, cex = 1.5)#plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = 'black')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = 'black')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x,diploidized_mean,col="#A1B5D8",lwd=5)
    lines(x, diploidized_mean+diploidized_sd, lty = 'dashed', col = '#A1B5D8')
    lines(x, diploidized_mean-diploidized_sd, lty = 'dashed', col = '#A1B5D8')
    lines(x,noPAFC_mean,col="#2B6F59",lwd=5)
    lines(x, noPAFC_mean+noPAFC_sd, lty = 'dashed', col = '#2B6F59')
    lines(x, noPAFC_mean-noPAFC_sd, lty = 'dashed', col = '#2B6F59')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "2B6F59", "DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

}


 ####### plot core non vs. dip fitness ########    
nrep <- 50
deme_plot <- 1000
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
file_list <- c("_heterozygosity")
edge_lab <- "Value"
for(file in file_list){
    print(file)
    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/no_PAFC/auto_bd_dr_exp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.005_ds--0.005_g-999999_start-5001")
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    rep_count <- 0
    for(i in c(1:100)){
        if(rep_count == nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]
        if(length(na.omit(pops[,1])) < 100){
            next
        }
        rep_count <- rep_count + 1
        for(gen in 1:length(pops[,1])){
            core_mat[gen,rep_count] <- dat[gen,1]
        }
        print(i)
    }

    auto_mean <- apply(core_mat, 1, mean, na.rm = TRUE)
    auto_sd <- apply(core_mat, 1, sd, na.rm = TRUE)


    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/dip/diploidization-diff_dipLambda-100_remDipMuts-1_meioticfFitness-1_pe_inflection-85_pe_slope-1_bd_dr_exp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-1.0e-06_bs-0.005_ds--0.005_g-999999_start-5001")
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    rep_count <- 0
    for(i in c(1:100)){
        if(rep_count == nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]
        if(length(na.omit(pops[,1])) < 100){
            next
        }
        rep_count <- rep_count + 1
        for(gen in 1:length(pops[,1])){
            core_mat[gen,rep_count] <- dat[gen,1]
        }
        print(i)
    }

    diploidized_mean <- apply(core_mat, 1, mean, na.rm = TRUE)
    diploidized_sd <- apply(core_mat, 1, sd, na.rm = TRUE)

    setwd("/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/no_PAFC/diploidization-diff_dipLambda-100_remDipMuts-1_meioticfFitness-0_pe_inflection-85_pe_slope-1_bd_dr_exp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-1.0e-06_bs-0.005_ds--0.005_g-999999_start-5001")
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    rep_count <- 0
    for(i in c(1:100)){
        if(rep_count == nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]
        if(length(na.omit(pops[,1])) < 100){
            next
        }
        rep_count <- rep_count + 1
        for(gen in 1:length(pops[,1])){
            core_mat[gen,rep_count] <- dat[gen,1]
        }
        print(i)
    }

    noPAFC_mean <- apply(core_mat, 1, mean, na.rm = TRUE)
    noPAFC_sd <- apply(core_mat, 1, sd, na.rm = TRUE)

    x <- as.integer(row.names(dat))
      
    png(file=paste(c("dip_not",file,"_byDeme_core.png"),sep="", collapse=""))
    plot((x),auto_mean, type="l", col="black", lwd=5, xlab="Generation", ylab="Nucleotide Diversity", las=1, cex.axis=1.5, cex.lab=1.5, ylim = c(min(c(auto_mean-auto_sd,diploidized_mean-diploidized_sd,noPAFC_mean-noPAFC_sd)),max(c(auto_mean+auto_sd,diploidized_mean+diploidized_sd,noPAFC_mean+noPAFC_sd))))
    par( mar=c(5,5,2,2), omi=c(0,0,0,0)) 
    #mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    #mtext(side=2, text="Heterozygosity", line=3.5, cex = 1.5)#plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = 'black')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = 'black')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x,diploidized_mean,col="#A1B5D8",lwd=5)
    lines(x, diploidized_mean+diploidized_sd, lty = 'dashed', col = '#A1B5D8')
    lines(x, diploidized_mean-diploidized_sd, lty = 'dashed', col = '#A1B5D8')
    lines(x,noPAFC_mean,col="#2B6F59",lwd=5)
    lines(x, noPAFC_mean+noPAFC_sd, lty = 'dashed', col = '#2B6F59')
    lines(x, noPAFC_mean-noPAFC_sd, lty = 'dashed', col = '#2B6F59')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()


    pdf(file=paste(c("dip_not",file,"_byDeme_core.pdf"),sep="", collapse=""))
    par( mar=c(5,5,2,2)) 
    plot((x),auto_mean, type="l", col="black", lwd=5, xlab="", ylab="", las=1, cex.axis=1.5, cex.lab=1.5, ylim = c(min(c(auto_mean-auto_sd,diploidized_mean-diploidized_sd,noPAFC_mean-noPAFC_sd)),max(c(auto_mean+auto_sd,diploidized_mean+diploidized_sd,noPAFC_mean+noPAFC_sd))))
    mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    mtext(side=2, text="Nucleotide Diversity", line=3.5, cex = 1.5)#plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
    lines(x, auto_mean+auto_sd, lty = 'dashed', col = 'black')
    lines(x, auto_mean-auto_sd, lty = 'dashed', col = 'black')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x,diploidized_mean,col="#A1B5D8",lwd=5)
    lines(x, diploidized_mean+diploidized_sd, lty = 'dashed', col = '#A1B5D8')
    lines(x, diploidized_mean-diploidized_sd, lty = 'dashed', col = '#A1B5D8')
    lines(x,noPAFC_mean,col="#2B6F59",lwd=5)
    lines(x, noPAFC_mean+noPAFC_sd, lty = 'dashed', col = '#2B6F59')
    lines(x, noPAFC_mean-noPAFC_sd, lty = 'dashed', col = '#2B6F59')
    #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
    #    col=c("Black", "2B6F59", "DCA92E"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

}

########## PAFC no PAFC fig w/ transparent lines etc


 
addTrans <- function(color,trans)
{
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.

  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))

  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}


nrep <- 50
deme_plot <- 1000
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
file_list <- c("_pairEfficiency")
edge_lab <- "Pairing Efficiency"
#mod_dir <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/pair_fit_dip_FIXED", full.names = TRUE, recursive = FALSE)
#file_list <- c("_fixedMutations","_meanBenMutationsPerGenome","_meanBenMutationsPerInd","_unscaledFitness","_maxPhs","_heterozygosity","_fixedMutations","_meanDelMutationsPerGenome","_meanDelMutationsPerInd","_meanFitnessScaled")
#for(mod in mod_dir){
    mod <- "/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/no_PAFC/diploidization-diff_dipLambda-100_remDipMuts-1_meioticfFitness-0_pe_inflection-85_pe_slope-1_bd_dr_exp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-1.0e-06_bs-0.005_ds--0.005_g-999999_start-5001"

    file <- file_list[1]
    setwd(mod)
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    edge_mat_di <- matrix(nrow=deme_plot/10,ncol = nrep)
    core_mat_di <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        dat2 <- read.csv(paste(c(i,"/",i,"_dipIndex",".csv"),sep="", collapse=""), row.names = 1)
        dat2 <- dat2[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]
        if(length(na.omit(pops[,1])) < 100){
            next
        }
        for(gen in 1:length(pops[,1])){
            edge_mat[gen,i] <- dat[gen,(max(which(pops[gen,] > 0)))]
            core_mat[gen,i] <- dat[gen,1]
            edge_mat_di[gen,i] <- dat2[gen,(max(which(pops[gen,] > 0)))]
            core_mat_di[gen,i] <- dat2[gen,1]

        }
        print(i)
    }

    edge_mean <- apply(edge_mat, 1, mean, na.rm = TRUE)
    edge_sd <- apply(edge_mat, 1, sd, na.rm = TRUE)

    core_mean <- apply(core_mat, 1, mean, na.rm = TRUE)
    core_sd <- apply(core_mat, 1, sd, na.rm = TRUE)

    x <- as.integer(row.names(dat))

   png(file=paste(c("dip_not",file,"edge_byDeme.png"),sep="", collapse=""))
   par( mar=c(5,5,2, 2))     #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    if(max(edge_mat_di[,1]) >= 0.9){
        plot(x,edge_mat[,1], type="l", col=addTrans("#A1B5D8",70), xlab="", ylab="", ylim = c(0,1), lwd = 3.0, las=1, cex.lab=2, cex.axis=1.5, cex.sub=2)
    } else {
        plot(x,edge_mat[,1], type="l", col=addTrans("#DCA92E",70), xlab="", ylab="", ylim = c(0,1), lwd = 3.0, las=1, cex.lab=2, cex.axis=1.5, cex.sub=2)
    }
    mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    mtext(side=2, text="Pairing Efficiency", line=3.5, cex = 1.5)
    for(gen in 2:length(edge_mat[1,])){
        if(max(edge_mat_di[,gen]) >= 0.90){
            lines(x,edge_mat[,gen], col = addTrans("#A1B5D8",70), lwd = 3.0)
        } else {
            lines(x,edge_mat[,gen], col = addTrans("#DCA92E",70), lwd = 3.0)
        }
    }
    #for(gen in 1:length(core_mat[1,])){
    #        lines(x,core_mat[,gen], col = addTrans("black",50))
    #}
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)
    dev.off() 

    pdf(file=paste(c("dip_not",file,"edge_byDeme.pdf"),sep="", collapse=""))
    par( mar=c(5,5,2, 2))     #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    if(max(edge_mat_di[,1]) >= 0.9){
        plot(x,edge_mat[,1], type="l", col=addTrans("#A1B5D8",70), xlab="", ylab="", ylim = c(0,1), lwd = 3.00, las=1, cex.lab=2, cex.axis=1.5, cex.sub=2)
    } else {
        plot(x,edge_mat[,1], type="l", col=addTrans("#DCA92E",70), xlab="", ylab="", ylim = c(0,1), lwd = 3.00, las=1, cex.lab=2, cex.axis=1.5, cex.sub=2)
    }
    mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    mtext(side=2, text="Pairing Efficiency", line=3.5, cex = 1.5)
    for(gen in 2:length(edge_mat[1,])){
        if(max(edge_mat_di[,gen]) >= 0.9){
            lines(x,edge_mat[,gen], col = addTrans("#A1B5D8",70), lwd = 3.0)
        } else {
            lines(x,edge_mat[,gen], col = addTrans("#DCA92E",70), lwd = 3.0)
        }
    }
    #for(gen in 1:length(core_mat[1,])){
    #        lines(x,core_mat[,gen], col = addTrans("black",50))
    #}
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)
    dev.off()

    png(file=paste(c("dip_not",file,"core_byDeme.png"),sep="", collapse=""))
    par( mar=c(5,5,2, 2))     #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
#plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    if(max(core_mat_di[,1]) >= 0.9){
        plot(x,core_mat[,1], type="l", col=addTrans("#A1B5D8",70), xlab="", ylab="", ylim = c(0,1), lwd = 3.00, las=1, cex.lab=2, cex.axis=1.5, cex.sub=2)
    } else {
        plot(x,core_mat[,1], type="l", col=addTrans("#DCA92E",70), xlab="", ylab="", ylim = c(0,1), lwd = 3.00, las=1, cex.lab=2, cex.axis=1.5, cex.sub=2)
    }
    mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    mtext(side=2, text="Pairing Efficiency", line=3.5, cex = 1.5)
    for(gen in 2:length(core_mat[1,])){
        if(max(core_mat_di[,gen]) >= 0.9){
            lines(x,core_mat[,gen], col = addTrans("#A1B5D8",70), lwd = 3.0)
        } else {
            lines(x,core_mat[,gen], col = addTrans("#DCA92E",70), lwd = 3.0)
        }
    }
    #for(gen in 1:length(core_mat[1,])){
    #        lines(x,core_mat[,gen], col = addTrans("black",50))
    #}
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)
    dev.off()

    pdf(file=paste(c("dip_not",file,"core_byDeme.pdf"),sep="", collapse=""))
    par( mar=c(5,5,2, 2))     #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
      if(max(core_mat_di[,1]) >= 0.9){
        plot(x,core_mat[,1], type="l", col=addTrans("#A1B5D8",70), xlab="", ylab="", ylim = c(0,1), lwd = 3.00, las=1, cex.lab=2, cex.axis=1.5, cex.sub=2)
    } else {
        plot(x,core_mat[,1], type="l", col=addTrans("#DCA92E",70), xlab="", ylab="", ylim = c(0,1), lwd = 3.00, las=1, cex.lab=2, cex.axis=1.5, cex.sub=2)
    }
    mtext(side=1, text="Generation", line=3.5, cex = 1.5)
    mtext(side=2, text="Pairing Efficiency", line=3.5, cex = 1.5)
    for(gen in 2:length(core_mat[1,])){
        if(max(core_mat_di[,gen]) >= 0.9){
            lines(x,core_mat[,gen], col = addTrans("#A1B5D8",70), lwd = 3.0)
        } else {
            lines(x,core_mat[,gen], col = addTrans("#DCA92E",70), lwd = 3.0)
        }
    }
    #for(gen in 1:length(core_mat[1,])){
    #        lines(x,core_mat[,gen], col = addTrans("black",50))
    #}
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)
    dev.off()
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))









    png(file=paste(c("all",file,"core_byDeme.png"),sep="", collapse=""), width = 1000, height = 500)
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    plot(x,core_mat[,1], type="l", col=addTrans("black",50), xlab="Generation", ylab=edge_lab, ylim = c(0,1))
    for(gen in 2:length(core_mat[1,])){
            lines(x,core_mat[,gen], col = addTrans("black",50))
    }
    #for(gen in 1:length(core_mat[1,])){
    #        lines(x,core_mat[,gen], col = addTrans("black",50))
    #}
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)
    dev.off()

    pdf(file=paste(c("all",file,"core_byDeme.pdf"),sep="", collapse=""))
        #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    plot(x,core_mat[,1], type="l", col=addTrans("black",50), xlab="Generation", ylab=edge_lab, ylim = c(0,1))
    for(gen in 2:length(core_mat[1,])){
            lines(x,core_mat[,gen], col = addTrans("black",50))
    }
    #for(gen in 1:length(core_mat[1,])){
    #        lines(x,core_mat[,gen], col = addTrans("black",50))
    #}
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)
    dev.off()






 
    png(file=paste(c("all",file,"_byDeme_edge.png"),sep="", collapse=""), width = 1000, height = 500)
    plot((x),edge_mean, type="l", col="gray", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(0,1))
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
    plot((x),edge_mean, type="l", col="gray", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(0,1))
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


###### PE and DI singles #####
nrep <- 100
deme_plot <- 1000
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
file_list <- c("_dipIndex")
edge_lab <- "Dip Index"

    mod <- "/work/users/w/w/wwbooker/polyploid_expansion_load/output/revision_run/dip/diploidization-diff_dipLambda-100_remDipMuts-1_meioticfFitness-1_pe_inflection-85_pe_slope-1_bd_dr_exp_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.002_rho-1.0e-06_bs-0.005_ds--0.005_g-999999_start-5001"
    file <- file_list[1]
    file2 <- "_pairEfficiency"
    setwd(mod)
    d_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    p_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    rep_count <- 0
    for(i in 1:100){
        if(rep_count == nrep){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = NULL)
        if(length(which(duplicated(dat[,1]))) > 0){
            next
        }
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        dat2 <- read.csv(paste(c(i,"/",i,file2,".csv"),sep="", collapse=""), row.names = 1)
        dat2 <- dat2[1:(deme_plot/10),]       
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]
        if(length(na.omit(pops[,1])) < 100){
            next
        }
        for(gen in 1:length(pops[,1])){
            d_mat[gen,i] <- dat[gen,(max(which(pops[gen,] > 0)))]
            p_mat[gen,i] <- dat2[gen,(max(which(pops[gen,] > 0)))]

        }
        rep_count <- rep_count + 1
        print(i)
        x <- as.integer(row.names(dat))
        png(paste(c(i,"/",i,"_dipIndex_PE_single.png"),sep="", collapse=""), width = 500, height = 500)
        par( mar=c(3,5,3,5)) 
        plot((x),d_mat[,i], type="l",las=1, col="gray", lwd=5, xlab="Generation", ylab="Pairing Efficiency", ylim = c(0,1))
        #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
        #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
        lines(x, p_mat[,i],col="black",lwd=5)
        axis(side = 4, at = seq(0,1,0.2), las=1)
        mtext("Diploidization Index", side = 4, line = 3)   
        #legend(1, 450, legend=c("Diploid", "Auto", "Allo"),
        #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
        #legend(1, 450, legend=c("Diploid", "Auto"),
        #     col=c("Black", "Blue"), lty=1, cex=0.8)

        dev.off()

        pdf(paste(c(i,"/",i,"_dipIndex_PE_single.pdf"),sep="", collapse=""))
        par( mar=c(3,5,3,5)) 
        plot((x),d_mat[,i], type="l",las=1, col="gray", lwd=5, xlab="Generation", ylab="Pairing Efficiency", ylim = c(0,1))
        #plot((x),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=edge_lab, ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd,allo_mean-allo_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd,allo_mean+allo_sd))))
        #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
        lines(x, p_mat[,i],col="black",lwd=5)
        axis(side = 4, at = seq(0,1,0.2), las=1)
        mtext("Diploidization Index", side = 4, line = 3)   
        #    col=c("Black", "#2B6F59", "#DCA92E"), lty=1, cex=0.8)
        #legend(1, 450, legend=c("Diploid", "Auto"),
        #     col=c("Black", "Blue"), lty=1, cex=0.8)

        dev.off()
    
    }

