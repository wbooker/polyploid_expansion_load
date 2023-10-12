require(plotrix)

nrep <- 10
deme_plot <- 2000

file_list <- c("_coreSiteFrequencies")
#file_list <- c("_fixedMutations","_meanBenMutationsPerGenome","_meanBenMutationsPerInd","_unscaledFitness","_maxPhs","_heterozygosity","_fixedMutations","_meanDelMutationsPerGenome","_meanDelMutationsPerInd","_meanFitnessScaled")


for(file in file_list){
    print(file)
    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/single_pop/auto_recessive_K-10000_m-0.005_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-1000")
    core_beg <- NULL
    core_end <- NULL
    for(i in 1:nrep){
        dat <- paste(c(i,"/",i,file,".csv"),sep="", collapse="")
        fields <- count.fields(dat, sep=",")
        dat <- read.table(dat, header = FALSE, sep = ",", col.names = paste0("V",seq_len(max(fields))), fill = TRUE)
        core_beg <- c(core_beg,as.numeric(dat[1,1:fields[1]]))
        core_end <- c(core_end,as.numeric(dat[length(dat[,1]),2:fields[length(fields)]]))
    }
}

auto_beg <- core_beg
auto_end <- core_end

for(file in file_list){
    print(file)
    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/single_pop/diploid_recessive_K-10000_m-0.005_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-2.5e-08_bs-0.0_ds--0.005_g-999999_start-1000")
    core_beg <- NULL
    core_end <- NULL
    for(i in 1:nrep){
        dat <- paste(c(i,"/",i,file,".csv"),sep="", collapse="")
        fields <- count.fields(dat, sep=",")
        dat <- read.table(dat, header = FALSE, sep = ",", col.names = paste0("V",seq_len(max(fields))), fill = TRUE)
        core_beg <- c(core_beg,as.numeric(dat[1,1:fields[1]]))
        core_end <- c(core_end,as.numeric(dat[length(dat[,1]),2:fields[length(fields)]]))
    }
}

dip_beg <- core_beg
dip_end <- core_end
len <- min(length(dip_beg), length(dip_end), length(auto_beg), length(auto_end))
l <- list(sample(dip_beg, len, replace = FALSE),sample(dip_end, len, replace = FALSE),sample(auto_beg, len, replace = FALSE),sample(auto_end, len, replace = FALSE))

    
png(file=paste(c("core_sfs_autoK10000_dipK10000.png"),sep="", collapse=""), width = 1000, height = 800)
#plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
multhist(l, col = c("black", "gray", "#0E3434", "#2B6F59"), main = "20k gens, single pop, auto k = 10000, dip k = 10000") 
legend("topright", legend=c("Diploid (start)", "Diploid (end)", "Autopolyploid (start)", "Autopolyploid (end)"),
    fill = c("black", "gray", "#0E3434", "#2B6F59"), cex=1.3)
#legend(1, 450, legend=c("Diploid", "Auto"),
#     col=c("Black", "Blue"), lty=1, cex=0.8)

dev.off()



nrep <- 5
deme_plot <- 19000 

file_list <- c("_meanDelMutationsPerGenome", "_meanDelMutationsPerInd")
#file_list <- c("_fixedMutations","_meanBenMutationsPerGenome","_meanBenMutationsPerInd","_unscaledFitness","_maxPhs","_heterozygosity","_fixedMutations","_meanDelMutationsPerGenome","_meanDelMutationsPerInd","_meanFitnessScaled")


for(file in file_list){
    print(file)
    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/single_pop/auto_recessive_K-500_m-0.005_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.0_ds--0.005_g-999999_start-1000")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            #edge_mat[gen,i] <- dat[gen,max(which(pops[gen,] > 0))]
            core_mat[gen,i] <- dat[gen,1]
        }
        print(i)
    }

    auto_mean <- apply(core_mat, 1, mean, na.rm = TRUE)
    auto_sd <- apply(core_mat, 1, sd, na.rm = TRUE)
    auto_mat <- core_mat
    auto_mat[is.na(auto_mat)] <- 0
    setwd("/proj/dschridelab/wwbooker/polyploid_expansion_load/output/single_pop/diploid_recessive_K-500_m-0.005_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_rho-1.0e-06_bs-0.0_ds--0.005_g-999999_start-1000")
    edge_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    core_mat <- matrix(nrow=deme_plot/10,ncol = nrep)
    for(i in 1:nrep){
        dat <- read.csv(paste(c(i,"/",i,file,".csv"),sep="", collapse=""), row.names = 1)
        dat <- dat[1:(deme_plot/10),]
        pops <- read.csv(paste(c(i,"/",i,"_popSize",".csv"),sep="", collapse=""), row.names = 1)
        pops <- pops[1:(deme_plot/10),]

        for(gen in 1:length(pops[,1])){
            #edge_mat[gen,i] <- dat[gen,max(which(pops[gen,] > 0))]
            core_mat[gen,i] <- dat[gen,1]
        }
        print(i)
    }

    dip_mean <- apply(core_mat, 1, mean, na.rm = TRUE)
    dip_sd <- apply(core_mat, 1, sd, na.rm = TRUE)
    dip_mat <- core_mat
    dip_mat[is.na(dip_mat)] <- 0
    
    auto_mean <- auto_mean[1:max(which(!is.nan(auto_mean)))]
    auto_sd <- auto_sd[1:max(which(!is.nan(auto_mean)))]
    auto_sd[is.na(auto_sd)] <- 0    
    dip_mean <- dip_mean[1:max(which(!is.nan(dip_mean)))]
    dip_sd <- dip_sd[1:max(which(!is.nan(dip_mean)))]
    dip_sd[is.na(dip_sd)] <- 0
    x1 <- as.integer(seq(0,20000,10))[1:length(auto_mean)]
    x2 <- as.integer(seq(0,20000,10))[1:length(dip_mean)]
    
    png(file=paste(c("SinglePopAutoDipK500",file,"_byDeme.png"),sep="", collapse=""), width = 1000, height = 500)
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    plot((x1),auto_mean, type="l", col="#2B6F59", lwd=5, xlab="Generation", ylab=file, ylim = c(min(c(auto_mean-auto_sd,dip_mean-dip_sd)),max(c(auto_mean+auto_sd,dip_mean+dip_sd))))
    lines(x1, auto_mean+auto_sd, lty = 'dashed', col = '#2B6F59')
    lines(x1, auto_mean-auto_sd, lty = 'dashed', col = '#2B6F59')
    #plot((x_diploid),y_diploid, type="l", col="black", lwd=5, xlab="Generation", ylab="Deme", main="Deme by Generation")
    lines(x2,dip_mean,col="black",lwd=5)
    lines(x2, dip_mean+dip_sd, lty = 'dashed', col = 'black')
    lines(x2, dip_mean-dip_sd, lty = 'dashed', col = 'black')
    legend("topright", legend=c("Auto", "Diploid"),
        col=c("#2B6F59", "Black"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)

    dev.off()

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

    x <- as.integer(seq(0,deme_plot-10,10))
    png(file=paste(c("SinglePopAutoDipK500_singles",file,"_byDeme.png"),sep="", collapse=""), width = 1000, height = 500)
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    plot(x,auto_mat[,1], type="l", col=addTrans("#2B6F59",50), xlab="Generation", ylab=file, ylim = c(min(c(auto_mat,dip_mat)),max(c(auto_mat,dip_mat))))
    for(gen in 2:length(auto_mat[1,])){
            lines(x,auto_mat[,gen], col = addTrans("#2B6F59",50))
    }
    for(gen in 1:length(dip_mat[1,])){
            lines(x,dip_mat[,gen], col = addTrans("black",50))
    }
    legend("topright", legend=c("Auto", "Diploid"),
        col=c("#2B6F59", "Black"), lty=1, cex=0.8)
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)
    dev.off()


}

 
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


nrep <- 10
deme_plot <- 1000
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
file_list <- c("_pairEfficiency")
edge_lab <- "Pair Efficiency"
mod_dir <- list.dirs(path = "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/pair_fit_dip_FIXED", full.names = TRUE, recursive = FALSE)
#file_list <- c("_fixedMutations","_meanBenMutationsPerGenome","_meanBenMutationsPerInd","_unscaledFitness","_maxPhs","_heterozygosity","_fixedMutations","_meanDelMutationsPerGenome","_meanDelMutationsPerInd","_meanFitnessScaled")
for(mod in mod_dir){
    #mod <- "/proj/dschridelab/wwbooker/polyploid_expansion_load/output/paper_run/pair_fit_dip_FIXED/test_recomb/no_peFit/diploidization-diff_dipLambda-100_remDipMuts-1_meioticfFitness-0_recessive_K-100_m-0.05_r-0.693147_u_del-2.5e-08_u_ben-2.5e-09_u_dip-0.001_rho-1.0e-06_bs-0.0_ds--0.01_g-999999_start-2501"
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


}






    png(file=paste(c("all",file,"core_byDeme.png"),sep="", collapse=""), width = 1000, height = 500)
    #plot((x),auto_mean, type="l", col="blue", lwd=5, xlab="Generation", ylab="Edge Value", ylim = c(min(c(auto_mean-auto_sd,diploid_mean-diploid_sd)),max(c(auto_mean+auto_sd,diploid_mean+diploid_sd))))
    plot(x,core_mat[,1], type="l", col=addTrans("black",50), xlab="Generation", ylab="Pairing Efficiency", ylim = c(0,1))
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
    plot(x,core_mat[,1], type="l", col=addTrans("black",50), xlab="Generation", ylab="Pairing Efficiency", ylim = c(0,1))
    for(gen in 2:length(core_mat[1,])){
            lines(x,core_mat[,gen], col = addTrans("black",50))
    }
    #for(gen in 1:length(core_mat[1,])){
    #        lines(x,core_mat[,gen], col = addTrans("black",50))
    #}
    #legend(1, 450, legend=c("Diploid", "Auto"),
    #     col=c("Black", "Blue"), lty=1, cex=0.8)
    dev.off()



}


 
    png(file=paste(c("all",file,"_byDeme_edge.png"),sep="", collapse=""), width = 1000, height = 500)
    plot((x),edge_mean, type="l", col="gray", lwd=5, xlab="Generation", ylab="Pairing Efficiency", ylim = c(0,1))
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
    plot((x),edge_mean, type="l", col="gray", lwd=5, xlab="Generation", ylab="Pairing Efficiency", ylim = c(0,1))
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