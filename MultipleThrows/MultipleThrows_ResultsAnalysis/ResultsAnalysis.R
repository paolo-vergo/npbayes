######################### RESULTS ANALYSIS ############################################

setwd("YOUR PATH/MultipleThrows/MultipleThrows_resulAnalysis")

#import MCMC
MU_MCMC<-read.csv(file="Mu_Estimates.csv", header=F)
SIGMA_MCMC<-read.csv(file="Sigma_Estimates.csv", header=F)
SIGMA_MCMC<-as.numeric(SIGMA_MCMC)

## Import saved structures
load("Activity_matrix.RData")
load("Times.RData")
load("True_Mean_matrix.RData")
load("Complete_Data.RData")

# compute estimates and store them
MU<-colMeans(MU_MCMC)
SIGMA<-mean(SIGMA_MCMC)

MU_complete<- matrix(data = NA, nrow = dim(Activity_matrix)[1], ncol = dim(Activity_matrix)[2])
k<-1
for( i in 1:dim(Activity_matrix)[2]){
  for ( j in 1:dim(Activity_matrix)[1]){
    if(Activity_matrix[j,i]==TRUE){
      MU_complete[j,i] <- MU[k]
      k<-k+1
    }
  }
}


# Approximated T intervals
alpha<-0.05
conf.int <- list() 
for (i in 1:dim(MU_MCMC)[2]){
  pred<-MU[i]
  stddev<-sd(MU_MCMC[,i])
  vec<- (-MU_MCMC[,i]+pred)/stddev
  r.quant<- quantile(vec, 1-alpha/2)
  l.quant<- quantile(vec, alpha/2)
  inf<- pred-r.quant*sqrt(SIGMA)
  sup<- pred-l.quant*sqrt(SIGMA)
  conf.int[[i]]<-c(inf,sup)
}
CONF<-list()
k<-1
for( i in 1:dim(Activity_matrix)[2]){
  CONF[[i]]<-list()
  for ( j in 1:dim(Activity_matrix)[1]){
    if(Activity_matrix[j,i]==TRUE){
      CONF[[i]][[j]] <- conf.int[[k]]
      k<-k+1
    }
    else{
      CONF[[i]][[j]]<-NA
    }
  }
}

## percentile intervals for sigma
r.quant<-quantile(SIGMA_MCMC,1-alpha/2)
l.quant<-quantile(SIGMA_MCMC,alpha/2)
sigma.conf<-c(l.quant,r.quant)


## Visualization: first 10 athletes

for (i in 1:10){
ts  <- c(Times[[1]][[i]],Times[[2]][[i]],Times[[3]][[i]],Times[[4]][[i]],Times[[5]][[i]])   # easily generalizable to athlete i
throws <- c(Complete_Data[[1]][[i]],Complete_Data[[2]][[i]],Complete_Data[[3]][[i]],Complete_Data[[4]][[i]],Complete_Data[[5]][[i]])

plot(ts,throws, xlim=c(0,1), ylim=c(-0.1,12), xlab='time',ylab='trials', main = paste("Athlete ",i))
abline(v=c(0.0,0.2,0.4,0.6,0.8,1.0),col='blue')
segments(0.0,True_Mean_matrix[i,1],0.2,True_Mean_matrix[i,1],col='green')
segments(0.2,True_Mean_matrix[i,2],0.4,True_Mean_matrix[i,2],col='green')
segments(0.4,True_Mean_matrix[i,3],0.6,True_Mean_matrix[i,3],col='green')
segments(0.6,True_Mean_matrix[i,4],0.8,True_Mean_matrix[i,4],col='green')
segments(0.8,True_Mean_matrix[i,5],1.0,True_Mean_matrix[i,5],col='green')

segments(0.0,MU_complete[i,1],0.2,MU_complete[i,1],col='red')
segments(0.2,MU_complete[i,2],0.4,MU_complete[i,2],col='red')
segments(0.4,MU_complete[i,3],0.6,MU_complete[i,3],col='red')
segments(0.6,MU_complete[i,4],0.8,MU_complete[i,4],col='red')
segments(0.8,MU_complete[i,5],1.0,MU_complete[i,5],col='red')

segments(0.0,CONF[[1]][[i]],0.2,CONF[[1]][[i]],col='blue',lty=2)
segments(0.2,CONF[[2]][[i]],0.4,CONF[[2]][[i]],col='blue',lty=2)
segments(0.4,CONF[[3]][[i]],0.6,CONF[[3]][[i]],col='blue',lty=2)
segments(0.6,CONF[[4]][[i]],0.8,CONF[[4]][[i]],col='blue',lty=2)
segments(0.8,CONF[[5]][[i]],1.0,CONF[[5]][[i]],col='blue',lty=2)

}

## Clustering between seasons
#devtools::install_github("sarawade/mcclust.ext")

clust_MCMC<-read.csv(file="Clust_Estimates.csv", header=F)
clust_MCMC<-clust_MCMC+1
library(mcclust.ext)

dim(clust_MCMC)
n <- dim(clust_MCMC)[1]

psm<-comp.psm(as.matrix(clust_MCMC))

image(psm)

heatmap(psm)

post_clust <- minVI(psm)

plot(post_clust$cl)

cl<-post_clust$cl

## post_clust_1<-minbinder(psm) best to use VI

#  Visualization

Clust_complete<- matrix(data = NA, nrow = dim(Activity_matrix)[1], ncol = dim(Activity_matrix)[2])
k<-1
for( i in 1:dim(Activity_matrix)[2]){
  for ( j in 1:dim(Activity_matrix)[1]){
    if(Activity_matrix[j,i]==TRUE){
      Clust_complete[j,i] <- cl[k]
      k<-k+1
    }
  }
}

## Plot of athletes 5-6-7-8-9
ts<-vector()
throws<-vector()
actual_cl<-vector()
pchs<-vector()
for ( i in 5:9){
  temp1 <- c(Times[[1]][[i]],Times[[2]][[i]],Times[[3]][[i]],Times[[4]][[i]],Times[[5]][[i]])
  ts    <- c(ts,temp1)
  temp2 <- c(Complete_Data[[1]][[i]],Complete_Data[[2]][[i]],Complete_Data[[3]][[i]],Complete_Data[[4]][[i]],Complete_Data[[5]][[i]]) 
  throws<- c(throws,temp2)
  for( j in 1:5){
    tmp<-rep(Clust_complete[i,j],length(Times[[j]][[i]]))
    actual_cl<-c(actual_cl,tmp)
    pchs<-c(pchs,rep(i,length(Times[[j]][[i]])))
  }
}

plot(ts,throws, xlim=c(0,1), ylim=c(-0.1,12.1),pch=pchs,col=actual_cl,xlab='time',ylab='trials', main = 'Athletes 5-9 with between-season clustering')
abline(v=c(0.0,0.2,0.4,0.6,0.8,1.0),col='darkgrey')


## Comparison of Athlete 1 vs Athlete 2

ts<-vector()
throws<-vector()
actual_cl<-vector()
pchs<-vector()
for ( i in 1:2){
  temp1 <- c(Times[[1]][[i]],Times[[2]][[i]],Times[[3]][[i]],Times[[4]][[i]],Times[[5]][[i]])
  ts    <- c(ts,temp1)
  temp2 <- c(Complete_Data[[1]][[i]],Complete_Data[[2]][[i]],Complete_Data[[3]][[i]],Complete_Data[[4]][[i]],Complete_Data[[5]][[i]]) 
  throws<- c(throws,temp2)
  for( j in 1:5){
    tmp<-rep(Clust_complete[i,j],length(Times[[j]][[i]]))
    actual_cl<-c(actual_cl,tmp)
    pchs<-c(pchs,rep(i,length(Times[[j]][[i]])))
  }
}

plot(ts,throws, xlim=c(0,1), ylim=c(-0.1,12.1),pch=pchs,col=actual_cl,xlab='time',ylab='trials', main = 'Athletes 1 vs Athlete 2')
abline(v=c(0.0,0.2,0.4,0.6,0.8,1.0),col='darkgrey')

segments(0.0,True_Mean_matrix[1,1]+0.05,0.2,True_Mean_matrix[1,1]+0.05,col='darkgreen',lty='dashed')
segments(0.2,True_Mean_matrix[1,2]+0.05,0.4,True_Mean_matrix[1,2]+0.05,col='darkgreen',lty='dashed')
segments(0.4,True_Mean_matrix[1,3]+0.05,0.6,True_Mean_matrix[1,3]+0.05,col='darkgreen',lty='dashed')
segments(0.6,True_Mean_matrix[1,4]+0.05,0.8,True_Mean_matrix[1,4]+0.05,col='darkgreen',lty='dashed')
segments(0.8,True_Mean_matrix[1,5]+0.05,1.0,True_Mean_matrix[1,5]+0.05,col='darkgreen',lty='dashed')

segments(0.0,True_Mean_matrix[2,1],0.2,True_Mean_matrix[2,1],col='green',lty='dashed')
segments(0.2,True_Mean_matrix[2,2],0.4,True_Mean_matrix[2,2],col='green',lty='dashed')
segments(0.4,True_Mean_matrix[2,3],0.6,True_Mean_matrix[2,3],col='green',lty='dashed')
segments(0.6,True_Mean_matrix[2,4],0.8,True_Mean_matrix[2,4],col='green',lty='dashed')
segments(0.8,True_Mean_matrix[2,5],1.0,True_Mean_matrix[2,5],col='green',lty='dashed')





