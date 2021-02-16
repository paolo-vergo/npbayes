library(distr)

setwd("C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/NIG_DataGeneration")

## simulated data:case 1
# X[1,1]...X[1,150]~0.2N(-3,0.1)+0.8N(0,0.5)
# X[2,1]...X[2,100]~0.1N(0,0.5)+0.9N(1,1.5))

mix1 <- UnivarMixingDistribution(Norm(mean=-3, sd=sqrt(0.1)), 
                                 Norm(mean=0, sd=sqrt(0.5)),
                                 mixCoeff=c( 0.2,0.8))
mix2 <- UnivarMixingDistribution(Norm(mean=0, sd=sqrt(0.5)), 
                                 Norm(mean=1, sd=sqrt(1.5)),
                                 mixCoeff=c(0.1, 0.9))

rmix1<-r(mix1)
rmix2<-r(mix2)

X1<-rmix1(100)
X2<-rmix2(100)

par(mfrow=c(2,2))
hist(X1, breaks=20, col="grey",prob=T, main="First season",xlim=c(-5,5))
hist(X2, breaks=20, col="grey",prob=T, main="Second season",xlim=c(-5,5))
plot(mix1,to.draw.arg='d',mfColRow=FALSE,xlim=c(-5,5))
plot(mix2,to.draw.arg='d',mfColRow=FALSE,xlim=c(-5,5))

dim1<-as.numeric(length(X1))
dim2<-as.numeric(length(X2))

Data<-c(X1,X2)
Dims<-c(dim1,dim2)

write.table(rbind(Data), file = "Data.csv",row.names= F, col.names= F,sep = ',')
write.table(rbind(Dims), file = "Dims.csv",row.names= F, col.names= F,sep = ',')
