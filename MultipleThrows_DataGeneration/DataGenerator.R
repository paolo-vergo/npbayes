############################# DATA GENERATION ################################################

setwd("C:/Users/aughi/Desktop/Hierarchical Bayesian Nonparametric models to smooth functional data/MultipleThrowsDataGeneration")


# Hypothesis:
# 5 seasons; total time is represented by 1, each season equispaced: length of 0.2
# 100 athletes;
# Number of trials per season are specified by a poisson of parameter 1.5
  # Season 1 and 2 are forced to contain at least one trial per athlete: Pois(1.5)+1
# times of trials pick values in the conrespondent season interval according to a uniform

# Output required by the Sampler:
#   1. Row vector containing all data, ordered by season.
#   2. Row vector containing the dimensions of each season (number of athletes participants)
#   3. Row vector specifying the number of trials per athlete in each season

# Note: 2 and 3 are necessary to parse data in the vector<vector<double>> structure of the sampler
# Note: List which keep track of which athletes is actually present (trials > 0) in
#       each season is necessary for post-proccesing of the result of the sampler.


############################# SET SEED ########################################


set.seed(210197)


# Underlying distributions in each season:


#    Season 1: Mixture of five gaussians:
#                     C=(2,4,6,8,10) 
#                     V=(v1,v2,v3) with vi iid U(0.5 +/- 0.05)
#                     W=(0.2,0.2,0.2,0.2,0.2)

N_trials_1 <- rpois(100,2.5)+1
X1         <- numeric(sum(N_trials_1))
Dim1       <- sum(N_trials_1>0)
Active_1   <- N_trials_1>0
ReducedN_1 <- numeric(Dim1)
mu_1       <- rep(NA,100) 
sigma_1    <- rep(NA,100)
times_1    <- list()
comp_data_1<- list()

j=1
for (i in 1:100){
  if (N_trials_1[i]==0){
    times_1[[i]]<- NA
    comp_data_1[[i]] <- NA
  }
  else { 
    
    mu_1[i] = sample( x = c(2,4,6,8,10),size = 1, prob = c(0.2,0.2,0.2,0.2,0.2) )
    sigma_1[i] = runif( 1, 0.45, 0.55)
    
    ReducedN_1[j]<-N_trials_1[i]
    j=j+1
    
    times_1[[i]] <- runif(N_trials_1[i], min=0.1, max= 0.19)
    
    if ( i==1){
      X1[i : N_trials_1[i]] = rnorm(N_trials_1[i], mean=mu_1[i], sd=sqrt(sigma_1[i]))
      comp_data_1[[i]] <- X1[i : N_trials_1[i]]
    }
    else {
      X1[sum(N_trials_1[1:i-1])+1 : N_trials_1[i]] = rnorm(N_trials_1[i], mean=mu_1[i], sd=sqrt(sigma_1[i]))
      comp_data_1[[i]] <- X1[sum(N_trials_1[1:i-1])+1 : N_trials_1[i]]
    }
  }
}


#    Season 2: Mixture of four gaussians:
#                     C=(2,4,6,10) 
#                     V=(v1,v2,v3,v4) with vi iid U(0.5 +/- 0.05)
#                     W=(0.2,0.2,0.3,0.3)

# specify n_trials for the season
N_trials_2 <- rpois(100,2.5)+1
X2         <- numeric(sum(N_trials_2))
Dim2       <- sum(N_trials_2>0)
Active_2   <- N_trials_2>0
ReducedN_2 <- numeric(Dim2)
mu_2       <- rep(NA,100) 
sigma_2    <- rep(NA,100)
times_2    <- list()
comp_data_2<- list()

j=1

for (i in 1:100){
  if (N_trials_2[i]==0){
    times_2[[i]] <- NA
    comp_data_2[[i]] <-NA
  }
  else { 
    
    mu_2[i] = sample( x = c(2,4,6,10),size = 1, prob = c(0.2,0.2,0.3,0.3) )
    sigma_2[i] = runif( 1, 0.45, 0.55)
    
    ReducedN_2[j]<-N_trials_2[i]
    j=j+1
    
    times_2[[i]] <- runif(N_trials_2[i],min = 0.21, max= 0.39)
    
    if ( i==1){
      X2[i : N_trials_2[i]] = rnorm(N_trials_2[i], mean=mu_2[i], sd=sqrt(sigma_2[i]))
      comp_data_2[[i]] <- X2[i : N_trials_2[i]]
    }
    else {
      X2[sum(N_trials_2[1:i-1])+1 : N_trials_2[i]] = rnorm(N_trials_2[i], mean=mu_2[i], sd=sqrt(sigma_2[i]))
      comp_data_2[[i]] <- X2[sum(N_trials_2[1:i-1])+1 : N_trials_2[i]]
    }
  }
}


#    Season 3: Mixture of three gaussians:
#                     C=(4,6,10) 
#                     V=(v1,v2,v3) with vi iid U(0.5 +/- 0.05)
#                     W=(0.3,0.4,0.3)

# specify n_trials for the season: from now on the athlete can make zero trials!
N_trials_3 <- rpois(100,2.5)
X3         <- numeric(sum(N_trials_3))
Dim3       <- sum(N_trials_3>0)
Active_3   <- N_trials_3>0
ReducedN_3 <- numeric(Dim3)
mu_3       <- rep(NA,100) 
sigma_3    <- rep(NA,100)
times_3    <- list()
comp_data_3<- list()

j=1
for (i in 1:100){
  if (N_trials_3[i]==0){
    times_3[[i]] <- NA
    comp_data_3[[i]] <-NA
  }
  else { 
    
    mu_3[i] = sample( x = c(4,6,10),size = 1, prob = c(0.3,0.4,0.3) )
    sigma_3[i] = runif( 1, 0.45, 0.55)
    
    ReducedN_3[j]<-N_trials_3[i]
    j=j+1
    
    times_3[[i]] <- runif(N_trials_3[i], min=0.41, max= 0.59)
    
    if ( i==1){
       X3[i : N_trials_3[i]] = rnorm(N_trials_3[i], mean=mu_3[i], sd=sqrt(sigma_3[i]))
       comp_data_3[[i]] <- X3[i : N_trials_3[i]]
    }
    else {
       X3[sum(N_trials_3[1:i-1])+1 : N_trials_3[i]] = rnorm(N_trials_3[i], mean=mu_3[i], sd=sqrt(sigma_3[i]))
       comp_data_3[[i]] <- X3[sum(N_trials_3[1:i-1])+1 : N_trials_3[i]]
    }
  }
}


#    Season 4: Mixture of three gaussians:
#                     C=(4,6,8) 
#                     V=(v1,v2,v3) with vi iid U(0.5 +/- 0.05)
#                     W=(0.4,0.2,0.4)

# specify n_trials for the season: from now on the athlete can make zero trials!
N_trials_4 <- rpois(100,2.5)
X4         <- numeric(sum(N_trials_4))
Dim4       <- sum(N_trials_4>0)
Active_4   <- N_trials_4>0
ReducedN_4 <- numeric(Dim4)
mu_4       <- rep(NA,100) 
sigma_4    <- rep(NA,100)
times_4    <- list()
comp_data_4<- list()

j=1
for (i in 1:100){
  if (N_trials_4[i]==0){
    times_4[[i]] <- NA
    comp_data_4[[i]] <-NA
  }
  else { 
    
    mu_4[i] = sample( x = c(4,6,8),size = 1, prob = c(0.4,0.2,0.4) )
    sigma_4[i] = runif( 1, 0.45, 0.55)
    
    ReducedN_4[j]<-N_trials_4[i]
    j=j+1
    
    times_4[[i]] <- runif(N_trials_4[i], min=0.61, max= 0.79)
    
    if ( i==1){
      X4[i : N_trials_4[i]] = rnorm(N_trials_4[i], mean=mu_4[i], sd=sqrt(sigma_4[i]))
      comp_data_4[[i]] <- X4[i : N_trials_4[i]]
    }
    else {
      X4[sum(N_trials_4[1:i-1])+1 : N_trials_4[i]] = rnorm(N_trials_4[i], mean=mu_4[i], sd=sqrt(sigma_4[i]))
      comp_data_4[[i]] <- X4[sum(N_trials_4[1:i-1])+1 : N_trials_4[i]]
    }
  }
}


#    Season 5: Mixture of five gaussians:
#                     C=(2,4,6,8,10) 
#                     V=(v1,v2,v3,v4,v5) with vi iid U(0.5 +/- 0.05)
#                     W=(0.1,0.2,0.4,0.2,0.1)

# specify n_trials for the season: from now on the athlete can make zero trials!
N_trials_5 <- rpois(100,2.5)
X5         <- numeric(sum(N_trials_5))
Dim5       <- sum(N_trials_5>0)
Active_5   <- N_trials_5>0
ReducedN_5 <- numeric(Dim5)
mu_5       <- rep(NA,100) 
sigma_5    <- rep(NA,100)
times_5    <-list()
comp_data_5<-list()

j=1
for (i in 1:100){
  if (N_trials_5[i]==0){
    times_5[[i]] <- NA
    comp_data_5[[i]] <- NA
  }
  else { 
    
    mu_5[i] = sample( x = c(2,4,6,8,10),size = 1, prob = c(0.1,0.2,0.4,0.2,0.1) )
    sigma_5[i] = runif( 1, 0.45, 0.55)
    
    ReducedN_5[j]<-N_trials_5[i]
    j=j+1
    
    times_5[[i]] <- runif(N_trials_5[i], min=0.81, max= 0.99)
    
    if ( i==1){
      X5[i : N_trials_5[i]] = rnorm(N_trials_5[i], mean=mu_5[i], sd=sqrt(sigma_5[i]))
      comp_data_5[[i]] <- X5[i : N_trials_5[i]]
    }
    else {
      X5[sum(N_trials_5[1:i-1])+1 : N_trials_5[i]] = rnorm(N_trials_5[i], mean=mu_5[i], sd=sqrt(sigma_5[i]))
      comp_data_5[[i]] <- X5[sum(N_trials_5[1:i-1])+1 : N_trials_5[i]]
    }
  }
}


### Build Data structure for the sampler and save them

Data             <- c(X1,X2,X3,X4,X5)
Dims             <- c(Dim1,Dim2,Dim3,Dim4,Dim5)
per_athlete_dims <- c(ReducedN_1,ReducedN_2,ReducedN_3,ReducedN_4,ReducedN_5)

write.table(rbind(Data), file = "Data.csv",row.names= F, col.names= F,sep = ',')
write.table(rbind(Dims), file = "Dims.csv",row.names= F, col.names= F,sep = ',')
write.table(rbind(per_athlete_dims), file = "per_athlete_dims.csv",row.names= F, col.names= F,sep = ',')


## Build Data structure to handle presence 

Activity_matrix <- data.frame(t(rbind(Active_1,Active_2,Active_3,Active_4,Active_5)))
names(Activity_matrix)<- c("Season_1","Season_2","Season_3","Season_4","Season_5")

## Build data structure to handle time of trials in each season per each athlete

Times <- list( "Season_1" = times_1,"Season_2" = times_2,"Season_3" = times_3,"Season_4" = times_4,"Season_5" = times_5)

#       example on how to access data structure:
        Times_per_athlete_season1 = Times["Season_1"]
        Times_athlete1_season1    = Times[[1]][[1]]
        Times_athlete1            = c(Times[[1]][[1]],Times[[2]][[1]],Times[[3]][[1]],Times[[4]][[1]],Times[[5]][[1]])

        
## Data structure for true means (and variances..)
        
True_Mean_matrix <- data.frame(t(rbind(mu_1,mu_2,mu_3,mu_4,mu_5)))
names(True_Mean_matrix)<- c("Season_1","Season_2","Season_3","Season_4","Season_5")

True_Var_matrix <- data.frame(t(rbind(sigma_1,sigma_2,sigma_3,sigma_4,sigma_5)))
names(True_Var_matrix)<- c("Season_1","Season_2","Season_3","Season_4","Season_5")

## Data structure for handling N_trials (with zeros, easier to access here in R)

N_trials_matrix <- data.frame(t(rbind(N_trials_1,N_trials_2,N_trials_3,N_trials_4,N_trials_5)))
names(True_Mean_matrix)<- c("Season_1","Season_2","Season_3","Season_4","Season_5")

## Data structure with complete Data

Complete_Data <- list( "Season_1" = comp_data_1,"Season_2" = comp_data_2,"Season_3" = comp_data_3,"Season_4" = comp_data_4,"Season_5" = comp_data_5)

## Plot: athlete 3

# data extraction:

Activity_matrix[3,]     # present in all seasons!
ts  <- c(Times[[1]][[3]],Times[[2]][[3]],Times[[3]][[3]],Times[[4]][[3]],Times[[5]][[3]])   # easily generalizable to athlete i
throws <- c(Complete_Data[[1]][[3]],Complete_Data[[2]][[3]],Complete_Data[[3]][[3]],Complete_Data[[4]][[3]],Complete_Data[[5]][[3]])

plot(ts,throws, xlim=c(0,1), ylim=c(-0.1,10.1), xlab='time',ylab='trials', main = 'Athlete 3')
abline(v=c(0.0,0.2,0.4,0.6,0.8,1.0),col='blue')
segments(0.0,True_Mean_matrix[3,1],0.2,True_Mean_matrix[3,1],col='red')
segments(0.2,True_Mean_matrix[3,2],0.4,True_Mean_matrix[3,2],col='red')
segments(0.4,True_Mean_matrix[3,3],0.6,True_Mean_matrix[3,3],col='red')
segments(0.6,True_Mean_matrix[3,4],0.8,True_Mean_matrix[3,4],col='red')
segments(0.8,True_Mean_matrix[3,5],1.0,True_Mean_matrix[3,5],col='red')


## save objects
save(Activity_matrix, file="Activity_matrix.RData")
save(Times, file="Times.RData")
save(True_Mean_matrix, file="True_Mean_matrix.RData")
save(Complete_Data, file="Complete_Data.RData")
