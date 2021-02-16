######################### RESULTS ANALYSIS ############################################

setwd("C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/Dataset/Second Run")

#import MCMC
MU_MCMC<-read.csv(file="Mu_Estimates.csv", header=F)
SIGMA_MCMC<-read.csv(file="Sigma_Estimates.csv", header=F)
SIGMA_MCMC<-as.numeric(SIGMA_MCMC)


## Import saved structures
setwd("C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/Dataset")
load("Activity_matrix.RData")
load("Times.RData")
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
plot(SIGMA_MCMC, type = 'l', col='darkgrey', main = 'Sigma Traceplot', sub = "alpha = 0.5    gamma = 0.7", xlab = '', ylab = 'Sigma')
abline(v=10,col='red')
text(5,3.0, labels=" Burnin: 10", col = 'red')
SIGMA_MCMC<-SIGMA_MCMC[-seq(1,200,by=1)]
plot(SIGMA_MCMC, type = 'l', col='darkgrey', main = 'Sigma Traceplot - Burnin: 10',sub = "alpha = 0.5   gamma = 0.7", xlab = '', ylab = 'Sigma')

##MU of chosen Athlete
j<- 500                                                     #athlete
x11()
par(mfrow=c(5,2))
for(i in 1:10){
  
  if(Activity_matrix[j,i]==T){
     plot(MUlist[[i]][[j]],type = 'l',col='darkgrey', main = paste('Mean Traceplot -','Athlete ',j,'Season ',i), xlab='index',ylab='Mu')
     abline(v=3,col='red')
  }
  else{
     plot(1:90,rep(17,90),type = 'n', main = paste('Mean Traceplot -','Athlete ',j,'Season ',i), xlab='index',ylab='Mu')
     text(45,17,labels = 'Athlete not partecipating in this season')
  }
}

MU_MCMC<-MU_MCMC[-seq(1,200,by=1),]
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


#Choose intervals to plot: 2 times sd or approximated t intervals

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

## Times Normalization 
for (i in 1:1086){
  for(k in 2:10){
    if(Activity_matrix[i,k]){
      if(min(Times[[k]][[i]])<365*(k-1)){
      addterm<- 366*(k-1)-min(Times[[k]][[i]])
         for (j in 1:10){
             Times[[j]][[i]]<-Times[[j]][[i]]+addterm
         }
      }
    }
  }
}
  
  
## Visualization

for (i in c(1,100,500)){
  ts  <- c(Times[[1]][[i]],Times[[2]][[i]],Times[[3]][[i]],Times[[4]][[i]],Times[[5]][[i]],
           Times[[6]][[i]],Times[[7]][[i]],Times[[8]][[i]],Times[[9]][[i]],Times[[10]][[i]])
  
  throws <- c(Complete_Data[[1]][[i]],Complete_Data[[2]][[i]],Complete_Data[[3]][[i]],Complete_Data[[4]][[i]],Complete_Data[[5]][[i]],
              Complete_Data[[6]][[i]],Complete_Data[[7]][[i]],Complete_Data[[8]][[i]],Complete_Data[[9]][[i]],Complete_Data[[10]][[i]])
  
  plot(ts,throws, xlim=c(-0.1, 3655), ylim=c(10,25), xlab='time',ylab='trials', main = paste("Athlete ",i))
  abline(v=c(0,365,365*2,365*3,365*4,365*5,365*6,365*7,365*8,365*9,365*10),col='darkgrey',lty='dashed')
  legend('bottomright', legend = c("Throws","Guessed Mean","95pc band"),cex=0.7, col = c('black','red','blue'), pch = c('o','-','--'))
  
  text(365/2,25,labels = "Season 1", col = 'darkgrey',cex = 0.7)
  text(365*3/2,25,labels = "Season 2", col = 'darkgrey',cex = 0.7)
  text(365*5/2,25,labels = "Season 3", col = 'darkgrey',cex = 0.7)
  text(365*7/2,25,labels = "Season 4", col = 'darkgrey',cex = 0.7)
  text(365*9/2,25,labels = "Season 5", col = 'darkgrey',cex = 0.7) 
  text(365*11/2,25,labels = "Season 6", col = 'darkgrey',cex = 0.7)
  text(365*13/2,25,labels = "Season 7", col = 'darkgrey',cex = 0.7)
  text(365*15/2,25,labels = "Season 8", col = 'darkgrey',cex = 0.7)
  text(365*17/2,25,labels = "Season 9", col = 'darkgrey',cex = 0.7)
  text(365*19/2,25,labels = "Season 10", col = 'darkgrey',cex = 0.7)
  
  segments(0,MU_complete[i,1],365,MU_complete[i,1],col='red')
  segments(365,MU_complete[i,2],365*2,MU_complete[i,2],col='red')
  segments(365*2,MU_complete[i,3],365*3,MU_complete[i,3],col='red')
  segments(365*3,MU_complete[i,4],365*4,MU_complete[i,4],col='red')
  segments(365*4,MU_complete[i,5],365*5,MU_complete[i,5],col='red')
  segments(365*5,MU_complete[i,6],365*6,MU_complete[i,6],col='red')
  segments(365*6,MU_complete[i,7],365*7,MU_complete[i,7],col='red')
  segments(365*7,MU_complete[i,8],365*8,MU_complete[i,8],col='red')
  segments(365*8,MU_complete[i,9],365*9,MU_complete[i,9],col='red')
  segments(365*9,MU_complete[i,10],365*10,MU_complete[i,10],col='red')
  
  segments(0,CONF[[1]][[i]],365,CONF[[1]][[i]],col='blue',lty=2)
  segments(365,CONF[[2]][[i]],365*2,CONF[[2]][[i]],col='blue',lty=2)
  segments(365*2,CONF[[3]][[i]],365*3,CONF[[3]][[i]],col='blue',lty=2)
  segments(365*3,CONF[[4]][[i]],365*4,CONF[[4]][[i]],col='blue',lty=2)
  segments(365*4,CONF[[5]][[i]],365*5,CONF[[5]][[i]],col='blue',lty=2)
  segments(365*5,CONF[[6]][[i]],365*6,CONF[[6]][[i]],col='blue',lty=2)
  segments(365*6,CONF[[7]][[i]],365*7,CONF[[7]][[i]],col='blue',lty=2)
  segments(365*7,CONF[[8]][[i]],365*8,CONF[[8]][[i]],col='blue',lty=2)
  segments(365*8,CONF[[9]][[i]],365*9,CONF[[9]][[i]],col='blue',lty=2)
  segments(365*9,CONF[[10]][[i]],365*10,CONF[[10]][[i]],col='blue',lty=2)
  
}

## Clustering between seasons
#devtools::install_github("sarawade/mcclust.ext")

setwd("C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/Dataset/Third Run")

clust_MCMC<-read.csv(file="Clust_Estimates.csv", header=F)
clust_MCMC<-clust_MCMC+1
library(mcclust.ext)

dim(clust_MCMC)
n <- dim(clust_MCMC)[1]

psm<-comp.psm(as.matrix(clust_MCMC))

#image(psm)

#heatmap(psm)

post_clust <- minVI(psm)

load("ClusterPal.Rdata")
palette(PAL)
plot(post_clust$cl, col=post_clust$cl,xlab='',ylab='', pch=4, cex=0.8, main='Cluster levels')
legend('bottomright', fill = palette(),cex=0.75,legend = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6","Cluster 7","Cluster 8","Cluster 9","Cluster 10","Cluster 11"))

cl<-post_clust$cl
x11()
## Cluster means
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
for ( i in 1:1086){
  temp1 <- c(Times[[1]][[i]],Times[[2]][[i]],Times[[3]][[i]],Times[[4]][[i]],Times[[5]][[i]],Times[[6]][[i]],Times[[7]][[i]],Times[[8]][[i]],Times[[9]][[i]],Times[[10]][[i]])
  ts    <- c(ts,temp1)
  temp2 <- c(Complete_Data[[1]][[i]],Complete_Data[[2]][[i]],Complete_Data[[3]][[i]],Complete_Data[[4]][[i]],Complete_Data[[5]][[i]],Complete_Data[[6]][[i]],Complete_Data[[7]][[i]],Complete_Data[[8]][[i]],Complete_Data[[9]][[i]],Complete_Data[[10]][[i]]) 
  throws<- c(throws,temp2)
  for( j in 1:10){
    tmp<-rep(Clust_complete[i,j],length(Times[[j]][[i]]))
    actual_cl<-c(actual_cl,tmp)
  }
}

plot(ts,throws,col=actual_cl,cex=0.4,xlim=c(0,4350),ylim=c(10,25),xlab='time',ylab='trials', main = 'All points - Clustering')
abline(v=c(0,365,365*2,365*3,365*4,365*5,365*6,365*7,365*8,365*9,365*10),col='darkgrey',lty='dashed')
legend('bottomright', fill = levels(as.factor(actual_cl)),cex=0.73,legend = c(paste("Cl. 1 - ", ClMeans[1]),paste("Cl. 2 - ", ClMeans[2]),paste("Cl. 3 - ", ClMeans[3]),paste("Cl. 4 - ", ClMeans[4]),paste("Cl. 5 - ", ClMeans[5]),paste("Cl. 6 - ", ClMeans[6]),paste("Cl. 7 - ", ClMeans[7]),paste("Cl. 8 - ", ClMeans[8]),paste("Cl. 9 - ", ClMeans[9]),paste("Cl. 10 - ", ClMeans[10]),paste("Cl. 11 - ", ClMeans[11])))
text(365/2,25,labels = "Season 1", col = 'darkgrey',cex = 0.63)
text(365*3/2,25,labels = "Season 2", col = 'darkgrey',cex = 0.63)
text(365*5/2,25,labels = "Season 3", col = 'darkgrey',cex = 0.63)
text(365*7/2,25,labels = "Season 4", col = 'darkgrey',cex = 0.63)
text(365*9/2,25,labels = "Season 5", col = 'darkgrey',cex = 0.63) 
text(365*11/2,25,labels = "Season 6", col = 'darkgrey',cex = 0.63)
text(365*13/2,25,labels = "Season 7", col = 'darkgrey',cex = 0.63)
text(365*15/2,25,labels = "Season 8", col = 'darkgrey',cex = 0.63)
text(365*17/2,25,labels = "Season 9", col = 'darkgrey',cex = 0.63)
text(365*19/2,25,labels = "Season 10", col = 'darkgrey',cex = 0.63)

## Histograms of different clusters in each season
library(ggplot2)
library(ggpubr)

setwd("C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/Dataset")
data = read.table("dataset_comp_with_n_seas.txt", sep = " ", header = T)
data$n_seas=data$n_seas+1
data<-data[which(data$n_seas<=10),]

setwd("C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/Dataset/second Run")

## select season and plot
j=1

  season<-data[which(data$n_seas==j),]
  cl_season<-rep(NA,length(season$result))
  for (i in 1:length(cl_season)){
    cl_season[i]<-Clust_complete[season$id[i],j]
  }
  d<-data.frame(throws = season$result, cl=factor(cl_season))
  p1<-ggplot(d, aes(x = throws, fill = cl)) +                       
    geom_histogram(position = "identity", alpha = 0.2, bins = 30) +
    labs(title=paste("Clustered Data Distributions in Season ",j), x="Throw", y="") +
    theme(plot.title = element_text(size=12))
  
  p2<-ggplot(d, aes(x = throws, fill = cl)) +                       
    geom_density(position = "identity", alpha = 0.2) +
    labs(title=paste("Cluster Estimated Distributions in Season ",j), x="Throw", y="")+
    theme(plot.title = element_text(size=12))
  ggarrange(p1,p2,ncol=2,nrow = 1)


par(mfrow=c(3,1))
for (i in c(1,100,150)){
  ts  <- c(Times[[1]][[i]],Times[[2]][[i]],Times[[3]][[i]],Times[[4]][[i]],Times[[5]][[i]],
           Times[[6]][[i]],Times[[7]][[i]],Times[[8]][[i]],Times[[9]][[i]],Times[[10]][[i]])
  
  throws <- c(Complete_Data[[1]][[i]],Complete_Data[[2]][[i]],Complete_Data[[3]][[i]],Complete_Data[[4]][[i]],Complete_Data[[5]][[i]],
              Complete_Data[[6]][[i]],Complete_Data[[7]][[i]],Complete_Data[[8]][[i]],Complete_Data[[9]][[i]],Complete_Data[[10]][[i]])
  
  plot(ts,throws, xlim=c(-0.1, 3655), ylim=c(10,25), xlab='time',ylab='trials', main = paste("Athlete ",i))
  abline(v=c(0,365,365*2,365*3,365*4,365*5,365*6,365*7,365*8,365*9,365*10),col='darkgrey',lty='dashed')
  
  text(365/2,25,labels = "Season 1", col = 'darkgrey',cex = 0.7)
  text(365*3/2,25,labels = "Season 2", col = 'darkgrey',cex = 0.7)
  text(365*5/2,25,labels = "Season 3", col = 'darkgrey',cex = 0.7)
  text(365*7/2,25,labels = "Season 4", col = 'darkgrey',cex = 0.7)
  text(365*9/2,25,labels = "Season 5", col = 'darkgrey',cex = 0.7) 
  text(365*11/2,25,labels = "Season 6", col = 'darkgrey',cex = 0.7)
  text(365*13/2,25,labels = "Season 7", col = 'darkgrey',cex = 0.7)
  text(365*15/2,25,labels = "Season 8", col = 'darkgrey',cex = 0.7)
  text(365*17/2,25,labels = "Season 9", col = 'darkgrey',cex = 0.7)
  text(365*19/2,25,labels = "Season 10", col = 'darkgrey',cex = 0.7)
  
  segments(0,MU_complete[i,1],365,MU_complete[i,1],col=Clust_complete[i,1])
  segments(365,MU_complete[i,2],365*2,MU_complete[i,2],col=Clust_complete[i,2])
  segments(365*2,MU_complete[i,3],365*3,MU_complete[i,3],col=Clust_complete[i,3])
  segments(365*3,MU_complete[i,4],365*4,MU_complete[i,4],col=Clust_complete[i,4])
  segments(365*4,MU_complete[i,5],365*5,MU_complete[i,5],col=Clust_complete[i,5])
  segments(365*5,MU_complete[i,6],365*6,MU_complete[i,6],col=Clust_complete[i,6])
  segments(365*6,MU_complete[i,7],365*7,MU_complete[i,7],col=Clust_complete[i,7])
  segments(365*7,MU_complete[i,8],365*8,MU_complete[i,8],col=Clust_complete[i,8])
  segments(365*8,MU_complete[i,9],365*9,MU_complete[i,9],col=Clust_complete[i,9])
  segments(365*9,MU_complete[i,10],365*10,MU_complete[i,10],col=Clust_complete[i,10])
}


## Step visualization for 3 specific athletes: 98,500, 800
pl<-800
ts  <- c(Times[[1]][[pl]],Times[[2]][[pl]],Times[[3]][[pl]],Times[[4]][[pl]],Times[[5]][[pl]],
         Times[[6]][[pl]],Times[[7]][[pl]],Times[[8]][[pl]],Times[[9]][[pl]],Times[[10]][[pl]])

throws <- c(Complete_Data[[1]][[pl]],Complete_Data[[2]][[pl]],Complete_Data[[3]][[pl]],Complete_Data[[4]][[pl]],Complete_Data[[5]][[pl]],
            Complete_Data[[6]][[pl]],Complete_Data[[7]][[pl]],Complete_Data[[8]][[pl]],Complete_Data[[9]][[pl]],Complete_Data[[10]][[pl]])


cl2<-c(rep(Clust_complete[pl,1],length(Times[[1]][[pl]])),rep(Clust_complete[pl,2],length(Times[[2]][[pl]])),
      rep(Clust_complete[pl,3],length(Times[[3]][[pl]])),rep(Clust_complete[pl,4],length(Times[[4]][[pl]])),
      rep(Clust_complete[pl,5],length(Times[[5]][[pl]])),rep(Clust_complete[pl,6],length(Times[[6]][[pl]])),
      rep(Clust_complete[pl,7],length(Times[[7]][[pl]])),rep(Clust_complete[pl,8],length(Times[[8]][[pl]])),
      rep(Clust_complete[pl,9],length(Times[[9]][[pl]])),rep(Clust_complete[pl,10],length(Times[[10]][[pl]])))

pp2<-data.frame(x=ts,y=throws)
d2=data.frame(x=c(0,365,365,365*2,365*2,365*3,365*3,365*4,365*4,365*5,365*5,365*6,365*6,365*7,365*7,365*8,365*8,365*9,365*9,365*10), 
             y=c(MU_complete[pl,1],MU_complete[pl,1],MU_complete[pl,2],MU_complete[pl,2],MU_complete[pl,3],MU_complete[pl,3],MU_complete[pl,4],MU_complete[pl,4],MU_complete[pl,5],MU_complete[pl,5],MU_complete[pl,6],MU_complete[pl,6],
                 MU_complete[pl,7],MU_complete[pl,7],MU_complete[pl,8],MU_complete[pl,8],MU_complete[pl,9],MU_complete[pl,9],MU_complete[pl,10],MU_complete[pl,10]))
P<-ggplot() +
  geom_step(data=d2, mapping=aes(x=x, y=y),size=1.0) +
  geom_vline(xintercept = c(0,365,365*2,365*3,365*4,365*5,365*6,365*7,365*8,365*9,365*10), linetype="dashed", 
               color = "darkgrey", size=0.5) +
  geom_point(data=pp2,aes(x=x, y=y, colour = factor(cl2))) +
  scale_color_manual("Clust",values=PAL[sort(unique(cl2[!is.na(cl2)]))])+
  ggtitle(paste("Athlete ",pl)) +
  xlab("Time [days]") +  
  ylab("Length of throw [m]")
nam<-paste("P",pl,sep = "")
assign(nam, P)


library(dplyr)

# create legend Chart
palette(PAL)
nn<-c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6","Cluster 7","Cluster 8","Cluster 9","Cluster 10","Cluster 11")
x<-data.frame(Cluster=nn, Approx_Mean=ClMeans,Color=rep('',11))
Leg<-ggtexttable(x,rows = NULL)
Leg <-table_cell_bg(Leg,row = 2,column = 3,fill = palette()[1])
Leg <-table_cell_bg(Leg,row = 3,column = 3,fill = palette()[2])
Leg <-table_cell_bg(Leg,row = 4,column = 3,fill = palette()[3])
Leg <-table_cell_bg(Leg,row = 5,column = 3,fill = palette()[4])
Leg <-table_cell_bg(Leg,row = 6,column = 3,fill = palette()[5])
Leg <-table_cell_bg(Leg,row = 7,column = 3,fill = palette()[6])
Leg <-table_cell_bg(Leg,row = 8,column = 3,fill = palette()[7])
Leg <-table_cell_bg(Leg,row = 9,column = 3,fill = palette()[8])
Leg <-table_cell_bg(Leg,row = 10,column = 3,fill = palette()[9])
Leg <-table_cell_bg(Leg,row = 11,column = 3,fill = palette()[10])
Leg <-table_cell_bg(Leg,row = 12,column = 3,fill = palette()[11])

ggarrange(P98,P500,P800,Leg, nrow = 2, ncol = 2,
          labels = c("","","",""))


##
# traceplot of unique dish values
n_clust<-rep(NA,1000)
for (i in 1:1000){
   n_clust[i]<-length(levels(as.factor(clust_MCMC[i,])))
}
plot(n_clust, type='l',main='Number of unique values across iterations',xlab='',sub = "Alpha = 0.5 - Gamma=0.7")

# traceplot of unique dish values between different seasons

Dims <- read.csv(file="C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/Dataset/Dims.csv",header = FALSE)
Dims<-as.numeric(Dims)
x11()
par(mfrow=c(5,2))
for( j in 1:10){
 n_clust<-rep(NA,200)
 for (i in 1:200){
  n_clust[i]<-length(levels(as.factor(clust_MCMC[i,1:Dims[j]])))
 }
 plot(n_clust, type='l',main=paste('Number of unique values across iterations - Season',j),xlab='',sub = "Alpha = 0.5 - Gamma=0.7")
 
}

#LPML
library(mvtnorm)
Data <- read.csv(file="C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/Dataset/Data.csv",header = FALSE)
Data<-as.numeric(Data)
cpo<-NULL
z<-0
count<-0
for (j in 1:10){                                  #season
  for (i in 1:1086 ){                            #athlete
    if (!is.na(Complete_Data[[j]][[i]])){
        dd<-NULL
        z<-z+1
        dato<-Complete_Data[[j]][[i]]
        for(k in 1:990){            #iteration
        mu<-rep(MU_MCMC[k,z],length(dato))    
        sig<-diag(length(dato))*SIGMA_MCMC[k]
        dd<-c(dd,1/dmvnorm(dato,mu,sig))
        count<-count+length(which(dd>1e1))
        dd<-dd[which(dd<1e1)]
        }
        cpo<-c(cpo,1/mean(dd))
    }
  }
}
LPML<-sum(log(cpo[!is.na(cpo)]))-count*1e-3

