library("lattice")
library(JMbayes)
library(JMbayes2)
library(dplyr)
library(readxl)
library(pec)
library(caret)
library(truncdist)
Range_T_linr <- function(TALL, DIST, X, U, variation, snrhigh){  # TALL = c(0,TS)
  u = U
  x1 = X$X1
  x2 = X$X2
  x3 = X$X3
  x4 = X$X4
  x5 = X$X5
  x6 = X$X6
  if (variation) {
    x3 = 0.5
    x4 = 0.5
    x6 = 1
  }
  
  if (snrhigh){
    Beta <- c(5,5,-5,-5,2.5,-2.5,-4)
  } else {
    Beta <- c(1,1,-1,-1,0.5,-0.84,0.5)
  }
  tlen = length(TALL)
  
  if(DIST == "Exp"){
    Lambda = 0.003
    Alpha = 0
    V = 0
    
    R0 = exp(Beta[1]*x1 +Beta[2]*x2 + Beta[3]*x3 + Beta[4]*x4 + Beta[5]*x5 + Beta[6]*x6 + Beta[7])
    R = Lambda*R0[-tlen] * (TALL[-1] - TALL[-tlen])
    
    VEC = c(0,cumsum(R),Inf)
    
    R.ID <- findInterval(-log(u), VEC)
    TT = -log(u) - VEC[R.ID] 
    TT = TT/Lambda/R0[R.ID] + TALL[R.ID]
  }else if(DIST == "WD"){
    Lambda = 0.012
    V = 0.8
    Alpha = 0
    R0 = exp(Beta[1]*x1 +Beta[2]*x2 + Beta[3]*x3 + Beta[4]*x4 + Beta[5]*x5 + Beta[6]*x6+ Beta[7])
    R = Lambda*R0[-tlen]*(TALL[-1]^V - TALL[-tlen]^V)
    
    VEC = c(0,cumsum(R),Inf)
    
    R.ID <- findInterval(-log(u), VEC)
    
    TT = -log(u) - VEC[R.ID] 
    TT = (TT/Lambda/R0[R.ID] + TALL[R.ID]^V)^(1/V)
    
  }else if(DIST == "WI"){
    Lambda = 0.001
    V = 2
    Alpha = 0
    R0 = exp(Beta[1]*x1 +Beta[2]*x2 + Beta[3]*x3 + Beta[4]*x4 + Beta[5]*x5 + Beta[6]*x6+ Beta[7])
    R = Lambda*R0[-tlen] * (TALL[-1]^V - TALL[-tlen]^V)
    
    VEC = c(0,cumsum(R),Inf)
    
    R.ID <- findInterval(-log(u), VEC)
    
    TT = -log(u) - VEC[R.ID] 
    TT = (TT/Lambda/R0[R.ID] + TALL[R.ID]^V)^(1/V)
    
  }else if(DIST == "Gtz"){
    Alpha = 0.1
    Lambda = 0.008
    V = 0
    R0 = exp(Beta[1]*x1 +Beta[2]*x2 + Beta[3]*x3 + Beta[4]*x4 + Beta[5]*x5 + Beta[6]*x6+ Beta[7])
    R = Lambda*R0[-tlen]/Alpha*(exp(Alpha*TALL[-1])-exp(Alpha*TALL[-tlen]))
    
    VEC = c(0,cumsum(R),Inf)
    
    R.ID <- findInterval(-log(u), VEC)
    TT = Alpha*(-log(u)-VEC[R.ID])/Lambda/R0[R.ID]
    TT = log(TT + exp(Alpha*TALL[R.ID]))/Alpha
    
  }
  
  result = list(Time = TT, Row = R.ID, Lambda = Lambda, Beta = Beta, Alpha = Alpha, V = V, Xi = R0)
  return(result)
}


#####################======== Large number of pseudo-subjects ===============#####################
Timevarying_PH_linear_gnrt <- function(N = 200, Distribution = "WI", censor.rate = 1,
                                       partial = TRUE, variation = FALSE, snrhigh = FALSE){
  npseu = 11
  Data <- as.data.frame(matrix(NA,npseu*N,27))
  names(Data)<-c("I","ID","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10", 
                 "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20",
                 "Start","Stop","C","Event","Xi")
  Data$ID <- rep(1:N,each=npseu)
  ## time-invariant
  Data$X1 <- rep(sample(c(0,1),N,replace=TRUE),each=npseu)
  Data$X2 <- rep(runif(N,0,1),each=npseu)
  
  Data$X7 <- rep(runif(N,0,1),each=npseu)
  Data$X8 <- rep(runif(N,1,2),each=npseu)
  Data$X9 <- rep(sample(c(1:3),N,replace=TRUE),each=npseu)
  Data$X10 <- rep(runif(N,0,1),each=npseu)
  Data$X11 <- rep(sample(c(0,1),N,replace=TRUE),each=npseu)
  Data$X12 <- rep(sample(c(0,1,2),N,replace=TRUE),each=npseu) 
  ## time-varying
  Data$X3 <- rbeta(npseu*N,2,5)
  Data$X4 <- runif(npseu*N,0,1)
  Data$X5 <- rbinom(npseu*N,5,0.46)
  Data$X5 <- sample(c(1:3),npseu*N,replace = TRUE)
  Data$X6 <- rnorm(npseu*N,0,1)
  Data$X14 <- sample(c(1:5),npseu*N,replace = TRUE)
  Data$X15 <- runif(npseu*N,0,1)
  Data$X17 <- runif(npseu*N,0,1)
  Data$X19 <- sample(c(0,1),npseu*N,replace=TRUE) 
  
  Count = 1
  tall <- rep(0, N)
  while(Count <= N){
    TS <- rep(0, npseu-1)
    if (Distribution == "Exp"){
      while(any(diff(sort(TS))<0.4)){
        TS <- sort(rtrunc(npseu-1, spec="beta", a=0.0001, b = 1, shape1 = 0.01, shape2 = 2))*900
      }
    } else if (Distribution == "WD"){
      while(any(diff(sort(TS))<0.4)){
        TS <- sort(rtrunc(npseu-1, spec="beta", a=0.0001, b = 1, shape1 = 0.001, shape2 = 2))*800
      }
    } else if (Distribution == "WI"){
      while(any(diff(sort(TS))<0.4)){
        TS <- sort(rtrunc(npseu-1, spec="beta", a=0.001, b = 1, shape1 = 0.05, shape2 = 5))*120
      }
    } else if (Distribution == "Gtz"){
      while(any(diff(sort(TS))<0.4)){
        TS <- sort(rtrunc(npseu-1, spec="beta", a=0.006, b = 5, shape1 = 0.05, shape2 = 5))*80
      }
    } else {
      stop("Wrong dististribution is given.")
    }   
    
    u = runif(1)
    k = runif(2)
    
    
#    x61 <- sample(c(0,1,2),1)
#    x62 <- sample(npseu,1)-1
#    Data[Data$ID==Count,]$X6[1:(1+x62)] <- rep(x61*(x61!=2) + 2*(x61==2), x62+1)
#    if (x61 == 0 && x62 != (npseu - 1)){
#      x63 <- sample(npseu-(1+x62),1)
#      Data[Data$ID==Count,]$X6[(1+x62+1):(1+x62+x63)] <- rep(1,x63)
#    }
    
    x13r <- sample(npseu-1,1)
    Data[Data$ID==Count,]$X13 <- c(numeric(x13r), rep(1,npseu-x13r))
    
    x16r <- sample(npseu-1,1)
    x161 <- sample(c(0,1),1)
    Data[Data$ID==Count,]$X16 <- c(rep(x161, x16r), rep(x161+1, npseu-x16r))
    
    x18r <- sort(sample(npseu-1,2))
    Data[Data$ID==Count,]$X18 <- c(rep(0, x18r[1]), rep(1,x18r[2]-x18r[1]), rep(2,npseu-x18r[2]))
    
    Data[Data$ID==Count,]$X20 <- k[1] * c(0,TS) + k[2]
    
    Data[Data$ID==Count,]$Start <- c(0,TS)
    Data[Data$ID==Count,]$Stop <- c(TS,NA)
    
    RT <- Range_T_linr(TALL = c(0, TS), 
                       DIST = Distribution, 
                       X = Data[Data$ID == Count, ], 
                       U = u, 
                       variation = variation, 
                       snrhigh = snrhigh)
    t = RT$Time
    rID = RT$Row
    Data[Data$ID==Count,]$Xi = RT$Xi
    tall[Count] = t
    if(rID==1){
      Data[Data$ID==Count,][1,]$Stop = t
      Data[Data$ID==Count,][1,]$Event = 1
    }else{
      Data[Data$ID==Count,][1:(rID-1),]$Event=0
      Data[Data$ID==Count,][rID,]$Event = 1
      Data[Data$ID==Count,][rID,]$Stop = t
    }
    #print(Count)
    Count = Count + 1
    
  }## end of the while loop
  DATA <- Data[!is.na(Data$Event),]
  rm(Data)
  gc()
  ###================== Add Censoring =========================================
  DATA$C <- 0
  
  if(length(unique(DATA$ID)) != N){
    stop("ID length NOT equal to N")
  }
  RET <- NULL
  
  if (variation){ 
    if (snrhigh) {
      if(Distribution == "WI"){
        if (censor.rate == 0){
          Censor.time <- rep(Inf, N)
        } else if (censor.rate == 1){
          Censor.time = rexp(N, rate = 1/63) ## NEW
        } else if (censor.rate == 2) { 
          Censor.time = rexp(N,rate = 1/7.5) ## NEW
        } else {
          stop("Wrong censoring type")
        } 
      } else {
        stop("Wrong distribution")
      }
    } else {
      if(Distribution == "WI"){
        if (censor.rate == 0){
          Censor.time <- rep(Inf, N)
        } else if (censor.rate == 1){
          Censor.time = rexp(N, rate = 1/78)
        } else if (censor.rate == 2) {
          Censor.time = rexp(N,rate = 1/20) ## NEW
        } else {
          stop("Wrong censoring type")
        } 
      } else {
        stop("Wrong distribution")
      }
    }
    
  } else { # Basic scenario
    if (snrhigh) {
      if(Distribution == "WI"){
        if (censor.rate == 0) {
          Censor.time <- rep(Inf,N)
        } else if (censor.rate == 1){
          Censor.time = rexp(N,rate = 1/78) 
        } else if (censor.rate == 2){
          Censor.time = rexp(N,rate = 1/6) # NEW
        } else{
          stop("Wrong censoring type")
        }
      } else {
        stop("Wrong distribution")
      }
    } else { # SNR original
      if (censor.rate == 0){
        Censor.time <- rep(Inf,N)
      }else if(censor.rate == 1){
        if(Distribution == "Exp"){
          Censor.time = rexp(N,rate = 1/805)
        }else if(Distribution == "WD"){
          Censor.time = rexp(N,rate = 1/440)
        }else if(Distribution == "WI"){
          Censor.time = rexp(N,rate = 1/82)
        }else if(Distribution == "Gtz"){
          Censor.time = rexp(N,rate = 1/75)
        }else {
          stop("Wrong distribution")
        }
      }else if(censor.rate == 2){ 
        if(Distribution == "WI"){
          Censor.time = rexp(N,rate = 1/22)
        } else {
          stop("Wrong distribution")
        }
      }else{
        stop("Wrong censoring type")
      }
    }
  }
  
  
  for( j in 1:length(unique(DATA$ID)) ){
    Vec <- c(0,DATA[DATA$ID==j,]$Stop,Inf)
    ID <- findInterval(Censor.time[j], Vec)
    
    if( ID <= nrow(DATA[DATA$ID==j,]) ){
      DATA[DATA$ID==j,][ID,]$C = 1
      DATA[DATA$ID==j,][ID,]$Event = 0
      DATA[DATA$ID==j,][ID,]$Stop = Censor.time[j]
      if( ID != nrow(DATA[DATA$ID==j,]) ){
        DATA[DATA$ID==j,][(ID+1):nrow(DATA[DATA$ID==j,]),]$Event = NA
      }
    }
  }
  
  Data <- DATA[!is.na(DATA$Event),]
  rm(DATA)
  gc()
  if(length(unique(Data$ID))!=N){
    stop("ID length NOT equal to N")
  }
  Data$I <- 1:nrow(Data)
  RET$fullData <- Data
  RET$fullData$Start <- round(RET$fullData$Start, 3)
  RET$fullData$Stop <- round(RET$fullData$Stop, 3)
  RET$fullData$Stop[RET$fullData$Start == RET$fullData$Stop] = RET$fullData$Stop[RET$fullData$Start == RET$fullData$Stop] + 0.001
  RET$Info = list(Coeff=list(Lambda = RT$Lambda, Alpha = RT$Alpha, Beta = RT$Beta, V = RT$V),
                  Dist=Distribution,
                  Set = "PH")
  
  DATA = data.frame(matrix(0,nrow=N,ncol=ncol(Data)))
  names(DATA) = names(Data)
  for (ii in 1:N){
    DATA[ii,] = Data[Data$ID==ii,][1,]
    ni = nrow(Data[Data$ID==ii,])
    if (ni >1){
      DATA[ii,]$Stop = Data[Data$ID==ii,]$Stop[ni]
      DATA[ii,]$Event = Data[Data$ID==ii,]$Event[ni]
    }
  }
  DATA$I = 1:nrow(DATA)
  RET$baselineData = DATA
  rm(DATA)
  
  if (partial){
    IDrev = rev(unique(Data$ID))
    # partialInfo = data.frame(matrix(0,nrow = N,ncol = 7))
    # names(partialInfo) = c("n","n_unobv","percT_unobv","wXiabsDiff","n_mto","percT_unobv_mto","wXiabsDiff_mto")
    for (Count in IDrev) {
      pN = sum(Data$ID == Count)
      I_count = Data[Data$ID==Count,]$I
      ## at least 50% of the times are not observable
      n_unobv = floor(pN*0.5)
      # partialInfo[Count,1]=pN # number of pseudo-subject in the full dataset
      # partialInfo[Count,2]=n_unobv # number of unobserved in the partial dataset
      if (pN > 1){
        unobv = sort(sample(1:(pN-1),n_unobv))+1
        tLast = Data[Data$ID==Count,]$Stop[pN]
        EventLast = Data[Data$ID==Count,]$Event[pN]
        
        T_unobv = Data[I_count[unobv],]$Stop
        dt_unobv = Data[I_count[unobv],]$Stop - Data[I_count[unobv],]$Start
        Xi_unobv = Data[I_count[unobv],]$Xi
        # partialInfo[Count,3]= sum(dt_unobv)/tLast
        
        Data = Data[-I_count[unobv],]
        if (pN > 2){
          Data[Data$ID==Count,]$Stop[1:(pN-n_unobv-1)] = Data[Data$ID==Count,]$Start[2:(pN-n_unobv)]
        } 
        Data[Data$ID==Count,]$Stop[pN-n_unobv] = tLast
        Data[Data$ID==Count,]$Event[pN-n_unobv] = EventLast
        
        # Xi_new_rID = findInterval(T_unobv, c(0,Data[Data$ID==Count,]$Stop),
        #                           left.open = TRUE, rightmost.closed = TRUE)
        # partialInfo[Count,4]=sum(abs(Data[Data$ID==Count,]$Xi[Xi_new_rID]-Xi_unobv)/Xi_unobv *dt_unobv)/sum(dt_unobv)
      }    
    }
    Data$I = 1:nrow(Data)
    RET$partialData = Data
    RET$partialData$Start <- round(RET$partialData$Start, 3)
    RET$partialData$Stop <- round(RET$partialData$Stop, 3)
    # RET$partialInfo = round(colMeans(partialInfo),digits = 3)
    # RET$partialInfo[5] = sum(partialInfo$n_unobv==0)
    # RET$partialInfo[6] = round(mean(partialInfo$percT_unobv[partialInfo$percT_unobv!=0]),digits = 3)
    # RET$partialInfo[7] = round(mean(partialInfo$wXiabsDiff[partialInfo$wXiabsDiff!=0]),digits = 3)
  }
  
  rm(Data)
  gc()
  
  return(RET)
  # return(tall)
}
set.seed(123)
df1 <- Timevarying_PH_linear_gnrt(N=200, Distribution = "WI", censor.rate = 1,
                                  partial = FALSE,variation = FALSE,snrhigh = FALSE)
sum(df1$baselineData$Event)/200
simdata <- df1$fullData[,c(2:8,23:26)]
ttp <- aggregate(Stop ~ ID, data = simdata, FUN = max)
simdata <- merge(simdata,ttp,by = "ID")
names(simdata)[names(simdata) == "Stop.x"] <- "time"
names(simdata)[names(simdata) == "Stop.y"] <- "ttp"

########## Hot encoding categorical variables ##############
simdata$X5 <- factor(simdata$X5)
dmy <- dummyVars( ~ X5,data = simdata)
simdata <- data.frame(simdata %>% select(-X5), predict(dmy, newdata = simdata))



######## Data split ############################################################
# Aggregate to get unique IDs with their event status
unique_individuals <- simdata %>%
  group_by(ID) %>%
  summarize(Event = max(Event), .groups = 'drop')  # Take max to ensure at least one event is considered

# Separate IDs into those with events and those without events
event_ids <- unique_individuals %>% filter(Event == 1)
non_event_ids <- unique_individuals %>% filter(Event == 0)

# Split event IDs into training and testing (75% training, 25% testing)
set.seed(132)  # Ensure reproducibility
train_event_ids_initial <- createDataPartition(event_ids$ID, p = 0.75, list = FALSE)  # 75% for initial training

# IDs for initial training and testing
train_event_ids_initial <- event_ids$ID[train_event_ids_initial]
test_event_ids <- event_ids$ID[!(event_ids$ID %in% train_event_ids_initial)]

# Combine initial training IDs with all non-event IDs
train_ids_initial <- c(train_event_ids_initial, non_event_ids$ID)
test_ids <- c(test_event_ids, setdiff(unique_individuals$ID, train_ids_initial))

# Split initial training data into final training (75%) and validation (25%)
train_ids_final <- createDataPartition(train_ids_initial, p = 0.75, list = FALSE)

# Final training and validation IDs
train_ids <- train_ids_initial[train_ids_final]
val_ids <- train_ids_initial[!(train_ids_initial %in% train_ids)]

# Create training, validation, and test datasets
training_data <- simdata %>% filter(ID %in% train_ids)
validation_data <- simdata %>% filter(ID %in% val_ids)
test_data <- simdata %>% filter(ID %in% test_ids)

# Check proportions of event in train, validation, and test data
sum(training_data$Event == 1) / sum(simdata$Event == 1)
sum(validation_data$Event == 1) / sum(simdata$Event == 1)
sum(test_data$Event == 1) / sum(simdata$Event == 1)

training_data$Event <- as.numeric(training_data$Event)
training_data.id <- training_data %>%
  group_by(ID) %>%
  filter(row_number() == n()) %>%
  select(ID, X1, X2, ttp, Event)

test_data.id <- test_data %>%
  group_by(ID) %>%
  filter(row_number() == n()) %>%
  select(ID, X1, X2, ttp, Event)

validation_data.id <- validation_data %>%
  group_by(ID) %>%
  filter(row_number() == n()) %>%
  select(ID, X1, X2, ttp, Event)

######### linear mixed model ###################################################

#training_data$time <- scale(training_data$time) # have to check if really 
# scaling the time variable is really a good idea.

fm1 <- lme(X4 ~ ns(time, 3) * X1 + ns(time,3)*X2, data = training_data,
           random = ~ 1 | ID, control = lmeControl(opt = 'optim'))

fm2 <- lme(X6 ~ ns(time, 3) * X1 + ns(time,3)*X2, data = training_data,
           random = ~ 1 | ID, control = lmeControl(opt = 'optim'))


fm3 <- lme(X3 ~ ns(time, 3) * X1 + ns(time,3)*X2, data = training_data,
           random = ~ 1 | ID, control = lmeControl(opt = 'optim'))

# fm4 <- mixed_model(X5 ~ time * X1 + time*X2, data = training_data,
#                    random = ~ 1 | ID, family = binomial())

fm4.1 <- mixed_model(X5.1 ~ time * X1 + time*X2, data = training_data,
                     random = ~ 1 | ID, family = binomial())

fm4.2 <- mixed_model(X5.2 ~ time * X1 + time*X2, data = training_data,
                    random = ~ 1 | ID, family = binomial())

fm4.3 <- mixed_model(X5.3 ~ time * X1 + time*X2, data = training_data,
                    random = ~ 1 | ID, family = binomial())

########## survival model ######################################################
training_data$Event <- as.numeric(training_data$Event)
training_data.id <- training_data %>%
  group_by(ID) %>%
  filter(row_number() == n()) %>%
  select(ID, X1, X2, ttp, Event)

test_data.id <- test_data %>%
  group_by(ID) %>%
  filter(row_number() == n()) %>%
  select(ID, X1, X2, ttp, Event)

CoxFit <- coxph(Surv(ttp, Event) ~ X1 + X2, data = training_data.id)

jointFit <- jm(CoxFit, list(fm1, fm2, fm3, fm4.1, fm4.2, fm4.3), time_var = "time")
summary(jointFit)

################################################################################
# Prediction
mean(test_data$time)
max_time <- 50
#max_time <- max(test_data$time)
summary(test_data$time)
plot(test_data$time,type ="l")
t0 <- mean(test_data.id$ttp) # follow up time up until this point is used 
# We want to do predictions for the longitudinal and survival outcomes for 
#the 40 individuals in the test data

ND <- test_data
ND <- ND[ND$time < t0,]
ND$Event <- 0  # We assume that all the individuals with time to e
# event less than the mean time haven't observed the event of interest
ND$ttp <- t0

# predictions for the longitudinal outcomes are 
predlong1 <- predict(jointFit, newdata = ND, return_newdata =TRUE) # this is the
# prediction of the longitudinal outcomes over the period we cut off
plot(predlong1)
# we now do the predictions over the enter time interval but can choose your own
# period
predlong2 <- predict(jointFit,newdata = ND, 
                     times = seq(t0, max_time,length.out = 165),
                     return_newdata = TRUE)
# plot for the 3 individuals 
for (i in test_data.id$ID[1:3]) {
  plot(predlong2, outcomes = 1, subject = i)
}

# We continue with the prediction for the event outcome
predSurv <- predict(jointFit, newdata = ND, process = "event",
                    times = seq(t0, max_time,length.out = 165),
                    return_newdata = TRUE)
plot(predlong2,predSurv) # As the follow up time increase, the cumulative risk 
# increases fast meaning the more likely a particular individual will experience
# the event. Note this plot is for the first subject.

cols <- c('#F25C78', '#D973B5', '#F28322')
plot(predlong2,predSurv,outcomes = 1:3, subject = 26,
     fun_event  = function(x) 1 - x,
     ylab_event = "Survival Probabilities",
     ylab_long = c("X4", "X6","X3"),
     bg = '#132743', col_points = cols, col_line_long = cols,
     col_line_event = '#F7F7FF', col_axis = "white", 
     fill_CI_long = c("#F25C7880", "#D973B580", "#F2832280"),
     fill_CI_event = "#F7F7FF80",
     pos_ylab_long = c(1.9, 1.9, 0.08))

################################################################################
# Prediction accuracy
validation_data$Event <- as.numeric(validation_data$Event != 0)
Dt <- max_time - t0
roc <- tvROC(jointFit, newdata = validation_data, Tstart = t0, Dt = Dt)
plot(roc)
tvAUC(roc)
confint(roc)

###############################################################################
roc <- tvROC(jointFit, newdata = test_data, Tstart = 2, Dt = 30)
roc
r <- tvAUC(roc)
r
IBC <- tvBrier(jointFit,newdata = validation_data,Tstart = t0, Dt = Dt, 
        integrated = TRUE, type_weights = "IPCW")
BS <- tvBrier(jointFit,newdata = validation_data,Tstart = t0, Dt = Dt, 
              integrated = FALSE, type_weights = "IPCW")