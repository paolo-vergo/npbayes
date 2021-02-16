############################### RESULT ANALYSIS ####################################
# library
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(wesanderson)
hrbrthemes::import_roboto_condensed()
col_pal=c("steelblue3", "darksalmon","lightgreen")

# Set wd
setwd("C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/NIG_ResultsAnalysis")

# Import Data
Data<-read.csv("C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/NIG_DataGeneration/Data.csv",header=F)
Data<-as.numeric(Data)
data<- data.frame(
  type = c( rep("1ST REST", 100), rep("2ND REST", 100)),
  value = c( Data[1:100], Data[101:200])
)

# Import Result: Predictive
res<-read.csv("Predictive.csv",header=F)

grid<-read.csv("Grid.csv",header=F)
grid<-as.numeric(grid)

predictive<-data.frame(GRID=grid, X1=t(res[1,]),X2=t(res[2,]))

p <-data %>%
  ggplot( ) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=20,aes(x=value, fill=type,y = ..density..)) +
  scale_fill_manual(values=col_pal) +
  geom_line(data = predictive, aes(x=GRID,y=X1), color= col_pal[1], size=1.1)+
  geom_line(data = predictive, aes(x = GRID, y = X2), col = col_pal[2],size=1.1) +
  ggtitle("alpha=10, gamma=5") +
  labs(fill="") +
  theme(plot.title = element_text(hjust = 0.5))

p

# Import Result: Predictive New Season

res<-read.csv("PredictiveNewRestaurant.csv",header=F)

data_pred<- data.frame(
  type = c( rep("DATA",length(Data) )),
  value = c(Data)
)

predictiveNew<-data.frame(GRID=grid,X1=t(res))


q <-data_pred %>%
    ggplot( ) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',binwidth = 0.4,aes(x=value, fill=type,y = ..density..)) +
    scale_fill_manual(values='darkorchid1') + 
    geom_line(data = predictiveNew, aes(x=GRID,y=X1), color= 'darkgoldenrod1', size=1.1) +
    labs(fill="") + 
    ggtitle("Predictive density for a new season") +
    theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")

q
