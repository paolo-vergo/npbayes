######################### RESULTS ANALYSIS ############################################

setwd("C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/SAVED")

#import MCMC
MU_MCMC<-read.csv(file="Mu_Estimates.csv", header=F)
SIGMA_MCMC<-read.csv(file="Sigma_Estimates.csv", header=F)
SIGMA_MCMC<-as.numeric(SIGMA_MCMC)
clust_MCMC<-read.csv(file="Clust_Estimates.csv", header=F)
clust_MCMC<-clust_MCMC+1

## Import saved structures
load("Activity_matrix.RData")
load("Times.RData")
load("True_Mean_matrix.RData")
load("Complete_Data.RData")



##Structures for traceplots
NAlist<-vector("list",1086)
for(i in 1:1086)
  NAlist[[i]]<-NA
MUlist <- list( "Season_1" = NAlist,"Season_2" = NAlist,"Season_3" = NAlist,"Season_4" = NAlist,"Season_5" = NAlist,
                "Season_6" = NAlist,"Season_7" = NAlist,"Season_8" = NAlist,"Season_9" = NAlist,"Season_10" = NAlist)
k<-1
for( i in 1:dim(Activity_matrix)[2]){
  for ( j in 1:dim(Activity_matrix)[1]){
    if(Activity_matrix[j,i]==TRUE){
      MUlist[[i]][[j]] <- MU_MCMC[,k] 
      k<-k+1
    }
  }
}

#Traceplots

#Sigma

plot(SIGMA_MCMC, type = 'l', col='darkgrey', main = 'Sigma Traceplot - Burnin: 200', xlab = '', ylab = 'Sigma')

##MU of chosen Athlete
j<- 20                                                     #athlete
par(mfrow=c(3,2))
for(i in 1:5){
  
  if(Activity_matrix[j,i]==TRUE){
    plot(MUlist[[i]][[j]],type = 'l',col='darkgrey', main = paste('Mean Traceplot -','Athlete ',j,'Season ',i), xlab='index',ylab='Mu')
  }
  else{
    plot(1:90,rep(17,90),type = 'n', main = paste('Mean Traceplot -','Athlete ',j,'Season ',i), xlab='index',ylab='Mu')
    text(45,17,labels = 'Athlete not partecipating in this season')
  }
}


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

window95pc<-TRUE

if(window95pc){
  
  # Gaussian Intervals: should contain 95% of data (wrt to single throws)
  alpha<-0.05
  conf.int <- list() 
  for (i in 1:dim(MU_MCMC)[2]){
    pred<-MU[i]
    inf<- pred-2*sqrt(SIGMA)
    sup<- pred+2*sqrt(SIGMA)
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
}  else {
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
}

## percentile intervals for sigma
r.quant<-quantile(SIGMA_MCMC,1-alpha/2)
l.quant<-quantile(SIGMA_MCMC,alpha/2)
sigma.conf<-c(l.quant,r.quant)

plot(SIGMA_MCMC, type = 'l', col='darkgrey', main = 'Sigma: Estimates & bands',sub = "alpha = 0.5   gamma = 0.7", xlab = '', ylab = 'Sigma')
abline(h=c(l.quant,SIGMA,r.quant),col=c("dodgerblue2","green","dodgerblue2"))

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

library(mcclust.ext)

dim(clust_MCMC)
n <- dim(clust_MCMC)[1]

psm<-comp.psm(as.matrix(clust_MCMC))

image(psm)

heatmap(psm)

post_clust <- minVI(psm)

plot(post_clust$cl)

cl<-post_clust$cl

ClMeans<-rep(NA,length(levels(as.factor(cl))))
for(i in 1:length(ClMeans)){
  s<-0
  l<-0
  for (j in 1:length(MU)){
    if(cl[j]==i){
      s<-s+MU[j]
      l<-l+1
    }
  }
  ClMeans[i]<-s/l
}

ClMeans<-round(ClMeans,digits=2)

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

## Plot of clustered points
ts<-vector()
throws<-vector()
actual_cl<-vector()
for ( i in 1:100){
  temp1 <- c(Times[[1]][[i]],Times[[2]][[i]],Times[[3]][[i]],Times[[4]][[i]],Times[[5]][[i]])
  ts    <- c(ts,temp1)
  temp2 <- c(Complete_Data[[1]][[i]],Complete_Data[[2]][[i]],Complete_Data[[3]][[i]],Complete_Data[[4]][[i]],Complete_Data[[5]][[i]]) 
  throws<- c(throws,temp2)
  for( j in 1:5){
    tmp<-rep(Clust_complete[i,j],length(Times[[j]][[i]]))
    actual_cl<-c(actual_cl,tmp)
  }
}

plot(ts,throws,col=actual_cl,cex=0.4,xlim=c(0,1),ylim=c(0,13),xlab='time',ylab='trials', main = 'All points - Clustering')
abline(v=c(0,0.2,0.4,0.6,0.8,1.0),col='darkgrey',lty='dashed')
legend('bottomright', fill = levels(as.factor(actual_cl)),cex=0.73,legend = c(paste("Cl. 1 - ", ClMeans[1]),paste("Cl. 2 - ", ClMeans[2]),paste("Cl. 3 - ", ClMeans[3]),paste("Cl. 4 - ", ClMeans[4]),paste("Cl. 5 - ", ClMeans[5])))
text(0.1,13,labels = "Season 1", col = 'darkgrey',cex = 0.63)
text(0.3,13,labels = "Season 2", col = 'darkgrey',cex = 0.63)
text(0.5,13,labels = "Season 3", col = 'darkgrey',cex = 0.63)
text(0.7,13,labels = "Season 4", col = 'darkgrey',cex = 0.63)
text(0.9,13,labels = "Season 5", col = 'darkgrey',cex = 0.63) 


## Intra seasonal histograms
library(ggplot2)
j=1

season<-NULL
for(i in 1:100)
  if(!is.na(Complete_Data[[j]][[i]]))
    season<-c(season, Complete_Data[[j]][[i]])

cl_season<-rep(NA,length(season))
k<-1
l<-0
for (i in 1:100){
  if(!is.na(Complete_Data[[j]][[i]])){
  l<-length(Complete_Data[[j]][[i]])
  cl_season[k:(k+l-1)]<-rep(Clust_complete[i,j],l)
  k<-k+l
  }
}
d1<-data.frame(throws = season, cl=factor(cl_season))
p1<-ggplot(d1, aes(x = throws, fill = cl)) +                       
  geom_histogram(position = "identity", alpha = 0.2, bins = 30) +
  labs(title=paste("Clustered Data Distributions in Season ",j), x="Throw", y="") +
  theme(plot.title = element_text(size=12))

j=2

season<-NULL
for(i in 1:100)
  if(!is.na(Complete_Data[[j]][[i]]))
    season<-c(season, Complete_Data[[j]][[i]])

cl_season<-rep(NA,length(season))
k<-1
l<-0
for (i in 1:100){
  if(!is.na(Complete_Data[[j]][[i]])){
    l<-length(Complete_Data[[j]][[i]])
    cl_season[k:(k+l-1)]<-rep(Clust_complete[i,j],l)
    k<-k+l
  }
}
d2<-data.frame(throws = season, cl=factor(cl_season))
p2<-ggplot(d2, aes(x = throws, fill = cl)) +                       
  geom_histogram(position = "identity", alpha = 0.2, bins = 30) +
  labs(title=paste("Clustered Data Distributions in Season ",j), x="Throw", y="") +
  theme(plot.title = element_text(size=12))

j=3

season<-NULL
for(i in 1:100)
  if(!is.na(Complete_Data[[j]][[i]]))
    season<-c(season, Complete_Data[[j]][[i]])

cl_season<-rep(NA,length(season))
k<-1
l<-0
for (i in 1:100){
  if(!is.na(Complete_Data[[j]][[i]])){
    l<-length(Complete_Data[[j]][[i]])
    cl_season[k:(k+l-1)]<-rep(Clust_complete[i,j],l)
    k<-k+l
  }
}
d3<-data.frame(throws = season, cl=factor(cl_season))
p3<-ggplot(d3, aes(x = throws, fill = cl)) +                       
  geom_histogram(position = "identity", alpha = 0.2, bins = 30) +
  labs(title=paste("Clustered Data Distributions in Season ",j), x="Throw", y="") +
  theme(plot.title = element_text(size=12))

j=4

season<-NULL
for(i in 1:100)
  if(!is.na(Complete_Data[[j]][[i]]))
    season<-c(season, Complete_Data[[j]][[i]])

cl_season<-rep(NA,length(season))
k<-1
l<-0
for (i in 1:100){
  if(!is.na(Complete_Data[[j]][[i]])){
    l<-length(Complete_Data[[j]][[i]])
    cl_season[k:(k+l-1)]<-rep(Clust_complete[i,j],l)
    k<-k+l
  }
}
d4<-data.frame(throws = season, cl=factor(cl_season))
p4<-ggplot(d4, aes(x = throws, fill = cl)) +                       
  geom_histogram(position = "identity", alpha = 0.2, bins = 30) +
  labs(title=paste("Clustered Data Distributions in Season ",j), x="Throw", y="") +
  theme(plot.title = element_text(size=12))

j=5

season<-NULL
for(i in 1:100)
  if(!is.na(Complete_Data[[j]][[i]]))
     season<-c(season, Complete_Data[[j]][[i]])

cl_season<-rep(NA,length(season))
k<-1
l<-0
for (i in 1:100){
  if(!is.na(Complete_Data[[j]][[i]])){
    l<-length(Complete_Data[[j]][[i]])
    cl_season[k:(k+l-1)]<-rep(Clust_complete[i,j],l)
    k<-k+l
  }
}
d5<-data.frame(throws = season, cl=factor(cl_season))
p5<-ggplot(d5, aes(x = throws, fill = cl)) +                       
  geom_histogram(position = "identity", alpha = 0.2, bins = 30) +
  labs(title=paste("Clustered Data Distributions in Season ",j), x="Throw", y="") +
  theme(plot.title = element_text(size=12))


library(ggpubr)
ggarrange(p1,p2,p3,p4,p5, labels = NULL,
          nrow = 3,ncol = 2)



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





