library(distr)

setwd("C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/NN_DataGeneration")

## simulated data:case 1
# X[1,1]...X[1,150]~0.33N(0,0.3)+0.33N(-2,0.3)+0.33N(2,0.3)
# X[2,1]...X[2,100]~0.5N(0,0.3)+0.5N(-2,0.3))
# X[3,1]...X[3,75]~0.5N(0,0.3)+0.5N(2,0.3))

mix1 <- UnivarMixingDistribution(Norm(mean=0, sd=sqrt(0.3)), 
                               Norm(mean=-2, sd=sqrt(0.3)),
                               Norm(mean=2, sd=sqrt(0.3)),
                                  mixCoeff=c( 1/3,1/3, 1/3))
mix2 <- UnivarMixingDistribution(Norm(mean=0, sd=sqrt(0.3)), 
                               Norm(mean=-2, sd=sqrt(0.3)),
                               mixCoeff=c(1/2, 1/2))
mix3 <- UnivarMixingDistribution(Norm(mean=0, sd=sqrt(0.3)), 
                               Norm(mean=2, sd=sqrt(0.3)),
                               mixCoeff=c(1/2, 1/2))
rmix1<-r(mix1)
rmix2<-r(mix2)
rmix3<-r(mix3)

X1<-rmix1(150)
X2<-rmix2(100)
X3<-rmix3(75)

#######################    data visualization   ################################

par(mfrow=c(2,3))
hist(X1, breaks=20, col="grey",prob=T, main="First season",xlim=c(-6,6))
hist(X2, breaks=20, col="grey",prob=T, main="Second season",xlim=c(-6,6))
hist(X3, breaks=20, col="grey",prob=T, main="Third season",xlim=c(-6,6))
plot(mix1,to.draw.arg='d',mfColRow=FALSE,xlim=c(-6,6))
plot(mix2,to.draw.arg='d',mfColRow=FALSE,xlim=c(-6,6))
plot(mix3,to.draw.arg='d',mfColRow=FALSE,xlim=c(-6,6))

##############################   save data      #################################

dim1<-as.numeric(length(X1))
dim2<-as.numeric(length(X2))
dim3<-as.numeric(length(X3))

Data<-c(X1,X2,X3)
Dims<-c(dim1,dim2,dim3)

write.table(rbind(Data), file = "Data2.csv",row.names= F, col.names= F,sep = ',')
write.table(rbind(Dims), file = "Dims2.csv",row.names= F, col.names= F,sep = ',')


