################### Real data preprocessing ###########################################

setwd("C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/Dataset")

data = read.table("dataset_comp_with_n_seas.txt", sep = " ", header = T)
data$n_seas=data$n_seas+1
dim(data)                                #55390 9
colnames(data)      #"id" "tij" "result" "gender" "age"  "dop" "enveironment" "season"       "n_seas" 
head(data)

attach(data)

levels(as.factor(id))         #1086 athletes
levels(as.factor(n_seas))     #27 seasons
table(as.factor(n_seas))

detach(data)

# Focus on first 10 seasons
data<-data[which(data$n_seas<=10),]

table(as.factor(data$n_seas))

levels(as.factor(data$id)) 

head(data)

max(data$tij)

## Data
Data<-NULL
for(i in 1:10){
  season<-data[which(data$n_seas==i),]
  Data<-c(Data, season$result)
}
## Dims
Dims<-rep(NA,10)
for(i in 1:10){
  season<-data[which(data$n_seas==i),]
  Dims[i]<-length(levels(as.factor(season$id)))
}

## per athlete dims
per_athlete_dims<-rep(NA, sum(Dims))
for(i in 1:10){
  season<-data[which(data$n_seas==i),]
  for(j in 1:Dims[i] ){
   per_athlete_dims[sum(Dims[1:i])-Dims[i]+j]<-table(as.factor(season$id))[[j]]
  }
}

## Generating data structure for the sampler
write.table(rbind(Data), file = "Data.csv",row.names= F, col.names= F,sep = ',')
write.table(rbind(Dims), file = "Dims.csv",row.names= F, col.names= F,sep = ',')
write.table(rbind(per_athlete_dims), file = "per_athlete_dims.csv",row.names= F, col.names= F,sep = ',')

## Activity matrix
Activity_matrix<-matrix(FALSE,nrow = max(Dims), ncol = length(Dims))
for(i in 1:10){
  season<-data[which(data$n_seas==i),]
  idx<-unique(season$id)
  Activity_matrix[idx,i]<-TRUE
}

## Times
NAlist<-vector("list",1086)
for(i in 1:1086)
  NAlist[[i]]<-NA
Times <- list( "Season_1" = NAlist,"Season_2" = NAlist,"Season_3" = NAlist,"Season_4" = NAlist,"Season_5" = NAlist,
               "Season_6" = NAlist,"Season_7" = NAlist,"Season_8" = NAlist,"Season_9" = NAlist,"Season_10" = NAlist)
for(i in 1:10){
  season<-data[which(data$n_seas==i),]
  idx<-unique(season$id)
  for (j in idx){
    athlete<-season[which(season$id==j),]
    Times[[i]][[j]]<-athlete$tij
  }
}

## Complete Data
NAlist<-vector("list",1086)
for(i in 1:1086)
  NAlist[[i]]<-NA
Complete_Data <- list( "Season_1" = NAlist,"Season_2" = NAlist,"Season_3" = NAlist,"Season_4" = NAlist,"Season_5" = NAlist,
               "Season_6" = NAlist,"Season_7" = NAlist,"Season_8" = NAlist,"Season_9" = NAlist,"Season_10" = NAlist)
for(i in 1:10){
  season<-data[which(data$n_seas==i),]
  idx<-unique(season$id)
  for (j in idx){
    athlete<-season[which(season$id==j),]
    Complete_Data[[i]][[j]]<-athlete$result
  }
}

# Visualization: Athlete 2
ts  <- c(Times[[1]][[2]],Times[[2]][[2]],Times[[3]][[2]],Times[[4]][[2]],Times[[5]][[2]],
         Times[[6]][[2]],Times[[7]][[2]],Times[[8]][[2]],Times[[9]][[2]],Times[[10]][[2]])

throws <- c(Complete_Data[[1]][[2]],Complete_Data[[2]][[2]],Complete_Data[[3]][[2]],Complete_Data[[4]][[2]],Complete_Data[[5]][[2]],
            Complete_Data[[6]][[2]],Complete_Data[[7]][[2]],Complete_Data[[8]][[2]],Complete_Data[[9]][[2]],Complete_Data[[10]][[2]]) 
plot(ts,throws, xlim=c(-0.1, 3655), ylim=c(0,20), xlab='time',ylab='trials', main = 'Athlete 2')
abline(v=c(0,365,365*2,365*3,365*4,365*5,365*6,365*7,365*8,365*9,365*10),col='blue')

#Save structures for Results Analysis
save(Activity_matrix, file="Activity_matrix.RData")
save(Times, file="Times.RData")

save(Complete_Data, file="Complete_Data.RData")

# Dataset quantities relevant for parameters setting

mean(data$result)
hist(data$result)
sd(data$result)
var(data$result)
