######################### DENSITY ESTIMATIOn ############################

### setwd and load libraries
setwd("Your_path/NN/NN_DataGeneration")

library(distr)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(wesanderson)

col_pal=c("lightgreen", "darksalmon","steelblue3")


par(mfrow=c(1,3))

# Import Data

Data<-read.csv("Data2.csv",header=F)
Data<-as.numeric(Data)
data<- data.frame(
  type = c( rep("1ST REST", 150), rep("2ND REST", 100), rep("3RD REST",75) ),
  value = c( Data[1:150], Data[151:250],Data[251:325])
)


# Import Result: Predictive

setwd("Your_path/NN/NN_ResultsAnalysis")
res<-read.csv("Predictive.csv",header=F)
grid<-read.csv("Grid.csv",header=F)
grid<-as.numeric(grid)
fit<-data.frame("GRID"=grid,"FIT1"=unlist(res[1,]),"FIT2"=unlist(res[2,]),"FIT3"=unlist(res[3,]))



p <- data %>%
  ggplot( ) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',binwidth = 0.3,aes(x=value, fill=type,y = ..density..)) +
  scale_fill_manual(values=col_pal) +geom_line(data = fit, aes(x=GRID,y=FIT1), color= col_pal[1], size=1.1)+labs(fill="")+geom_line(data = fit, aes(x = GRID, y = FIT2), col = col_pal[2],size=1.1)+geom_line(data = fit, aes(x = GRID, y = FIT3), col =col_pal[3],size=1.1)+ ggtitle("alpha = 1, gamma=1") +
  theme(plot.title = element_text(hjust = 0.5))

p


