######################### DENSITY ESTIMATIOn ############################

### LOAD DATA
setwd("C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/NN_DataGeneration")

par(mfrow=c(1,3))
Data<-read.csv("Data2.csv",header=F)
Data<-as.numeric(Data)
hist(Data[1:150], breaks=20, col="grey",prob=T, main="First season",xlim=c(-6,6))
hist(Data[151:250], breaks=20, col="grey",prob=T, main="Second season",xlim=c(-6,6))
hist(Data[251:325], breaks=20, col="grey",prob=T, main="Third season",xlim=c(-6,6))

setwd("C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/NN_ResultsAnalysis")

## Density Estimation 
res<-read.csv("Predictive.csv",header=F)

grid<-read.csv("Grid.csv",header=F)
grid<-as.numeric(grid)

layout(1)
plot(grid, res[1,],type='l')
plot(grid, res[2,],type='l')
plot(grid, res[3,],type='l')


hist(Data[1:150], breaks=30, col="grey",prob=T, main="First season",xlim=c(-6,6),ylim=c(0,0.9))
lines(grid, res[1,],col='red')
hist(Data[151:250], breaks=30, col="grey",prob=T, main="Second season",xlim=c(-6,6),ylim=c(0,0.9))
lines(grid, res[2,],col='red')
hist(Data[251:325], breaks=30, col="grey",prob=T, main="Third season",xlim=c(-6,6),ylim=c(0,0.9))
lines(grid, res[3,],col='red')
