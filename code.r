library(bigmemory)
library(biganalytics)
library(survival)
library(dplyr)
library(plyr)
library(XLConnect)
library(chron)
library(survMisc)
library(glmnet)
library(MASS)
library(pec)
library(quantreg)
library(autopls)
library(jpeg)
library(rms)

##IMPORT DATA. BEFORE THIS, SET WD TO THE LOCATION WHERE FILES ARE STORED
data<-read.csv("CelticData29AUG16.csv",header=T,sep=",",na.strings=c(""," ","NA"))
train <- read.csv("train.csv",header=T,sep=',',na.strings=c(""," ","NA"))

#LET R READ DATES
data$Outcome.Date<-as.character.Date(data$Outcome.Date)
data$Assess.Date<-as.character.Date(data$Assess.Date)
data$V39_1<-as.character.Date(data$V39_1)
data$Outcome.Date<-as.Date(data$Outcome.Date,format = '%m/%d/%Y')
data$Assess.Date<-as.Date(data$Assess.Date,format = '%m/%d/%Y')
data$V39_1<-as.Date(data$V39_1,format = '%m/%d/%Y')

train$Outcome.Date<-as.character.Date(train$Outcome.Date)
train$Assess.Date<-as.character.Date(train$Assess.Date)
train$V39_1<-as.character.Date(train$V39_1)
train$Outcome.Date<-as.Date(train$Outcome.Date,format = '%m/%d/%Y')
train$Assess.Date<-as.Date(train$Assess.Date,format = '%m/%d/%Y')
train$V39_1<-as.Date(train$V39_1,format = '%m/%d/%Y')

#OTHER DATA MANIPULATIONS
train <- train[order(train[,1],train[,3],train[,2]),]
train<-unique.matrix(train)
#ANOTHER METHOD TO DO THIS
# v=NULL
# for (i in 1:(nrow(train)-1))
# {
#   if ((train[i,2]==train[i+1,2]))
#   {
#     if((train[i,3]==train[i+1,3]))
#     {
#       v<-c(v,i)
#     }
#   }
# }
# train <- train(-v,)
train$start=0
train$stop=0
train$outcome=0
train$Outcome.week <- abs(difftime(train$Outcome.Date,Sys.Date(),units="weeks"))
train$risk <- NULL
#A BACKUP FOR TRAIN
train1 = train
train = train1
unique.id=c(unique(train$X))

#DEVELOP START, STOP AND EVENT FOR THE Surv() FUNCTION
for (i in unique.id) 
{
  unique.out <- unique(train[which(train$X==i,arr.ind = T),3])
  for (z in 1:length(unique.out))
  {
    temp <-(which(train$X==i & train$Outcome.Date==unique.out[z],arr.ind = T))
    if (length(temp)==1)
    {
      if(train[temp,2]==train[temp,3])
      {
        train[temp,208]=0.1428571
        train[temp,209]=1
      }
      else
      {
        train[temp,208]=abs(difftime(train[temp,3],train[temp,2],units="weeks"))
        train[temp,209]=1
      }
    }
    else
    {
      temp1 <- temp[-length(temp)]
      for(j in temp1) 
      {
        if(train[j,2]==train[j,3])
        {
          train[j,208]=0.1428571
          train[j,209]=1
        }
        else
        {
          train[j, 207] = difftime(train[j,2], train[temp1[1],2], units = "weeks")
          train[j, 208] = difftime(train[j+1,2], train[temp1[1],2], units = "weeks")
          if(j==temp[1]&&train[j,3]<train[j,2])
          {
            train[j,209]<-1
          }
          else
          {
            if(train[j,3]>=train[j,2]&&train[j,3]<train[j+1,2])
            {
              train[j,209]<-1
            }
            else
            {
              if (j==temp[length(temp)]&&train[j+1,3]>=train[j+1,2])
              {
                train[j+1,209]<-1
              }
            }
          }
        }
        # if( train[j, 3] >= train[j, 2] )
        # {
        #   train[j, 209] = 1
        # }
        #temp_diff <- temp_diff + difftime(train[j, 208], train[j, 207], units = "weeks")
      }
      last_ind = temp[length(temp)]
      train[last_ind, 207] = difftime(train[last_ind, 2], train[temp1[1], 2], units = "weeks")
      train[last_ind, 208] = train[last_ind, 207]+(train[last_ind,207]/length(temp1))
      
    }
  }
}

#REMOVE COVARIATES THAT HAVE MORE THAN 20% MISSING DATA. THIS CAUSES BIG PROGLEMS WITH THE COX PROPORTIONAL HAZARDS 
#REGRESSION
f <- function(x) {
  sum(is.na(x)) < length(x) * 0.2
}
train1<-train[, vapply(train, f, logical(1)), drop = F]

#A BACKUP FOR TRAIN1
trainF=train1

###USER DEFINED FUNCTION FOR RISK CALCULAITON
risk = function(model, newdata, time) {
  as.numeric(1-summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = time)$surv)
}

#DEVELOP A COXPH REGRESSION MODEL AND CALCULATE A 13 FOLD MSE
coxmod <- vector("list",5)
c <- (split(c(5:119),sample(c(5:119),size=5,replace=FALSE)))
err <- NULL
t <- split(c(1:nrow(train1)), sample(1:nrow(train1), size=13, replace=FALSE))
for (i in 1:13)
{
  train1 <- trainF[-t[[i]],]
  val <- trainF[t[[i]],]
  r <- matrix(rep(0),nrow(val),6)
  for(i in 1:5)
  {
    cols <- colnames(train1[unlist(c[i])])
    my.formula <- as.formula(paste( "Surv(start,stop,outcome)", '~', paste( cols, collapse=' + ' ) ))
    coxmod[[i]] <- coxph(my.formula,data=na.omit(train1))
    step <- step(coxmod[[i]],direction="backward",trace=F)
    cols <- names(step$assign)
    my.formula <- as.formula(paste( "Surv(start,stop,outcome)", '~', paste( cols, collapse=' + ' ) ))
    coxmod[[i]] <- coxph(my.formula,data=na.omit(train1))
    for (j in 1:nrow(val))
    {
      r[j,i] <- risk(coxmod[[i]],val[j,],val[j,123])
    }
  }
  r[,6] <- apply(r[,1:5],1,function (x) mean(x,na.rm=T))
    for (k in 1:nrow(val))
    {
      ifelse(r[k,6]>=0.5,val[k,124]<-"HIGH RISK",ifelse(r[k,6]>=0.25,val[k,124]<-"MEDIUM RISK",val[k,124]<-"LOW RISK"))    
    }
  err <- c(err, length(which(val[,124]=="LOW RISK",arr.ind = T)))
}
mse_coxmod <- mean(err)/nrow(val)


##USER DEFINED FUNCTION FOR SURVIVAL CURVE PLOT (ACQUIRED FROM OTHER SOURCES)
ggsurv <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                   cens.col = 'red', lty.est = 1, lty.ci = 2,
                   cens.shape = 3, back.white = F, xlab = 'Time',
                   ylab = 'Survival', main = ''){
  
  library(ggplot2)
  strata <- ifelse(is.null(s$strata) ==T, 1, length(s$strata))
  stopifnot(length(surv.col) == 1 | length(surv.col) == strata)
  stopifnot(length(lty.est) == 1 | length(lty.est) == strata)
  
  ggsurv.s <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = ''){
    
    dat <- data.frame(time = c(0, s$time),
                      surv = c(1, s$surv),
                      up = c(1, s$upper),
                      low = c(1, s$lower),
                      cens = c(0, s$n.censor))
    dat.cens <- subset(dat, cens != 0)
    
    col <- ifelse(surv.col == 'gg.def', 'black', surv.col)
    
    pl <- ggplot(dat, aes(x = time, y = surv)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(col = col, lty = lty.est)
    
    pl <- if(CI == T | CI == 'def') {
      pl + geom_step(aes(y = up), color = col, lty = lty.ci) +
        geom_step(aes(y = low), color = col, lty = lty.ci)
    } else (pl)
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  
  ggsurv.m <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = '') {
    n <- s$strata
    
    groups <- factor(unlist(strsplit(names
                                     (s$strata), '='))[seq(2, 2*strata, by = 2)])
    gr.name <-  unlist(strsplit(names(s$strata), '='))[1]
    gr.df <- vector('list', strata)
    ind <- vector('list', strata)
    n.ind <- c(0,n); n.ind <- cumsum(n.ind)
    for(i in 1:strata) ind[[i]] <- (n.ind[i]+1):n.ind[i+1]
    
    for(i in 1:strata){
      gr.df[[i]] <- data.frame(
        time = c(0, s$time[ ind[[i]] ]),
        surv = c(1, s$surv[ ind[[i]] ]),
        up = c(1, s$upper[ ind[[i]] ]),
        low = c(1, s$lower[ ind[[i]] ]),
        cens = c(0, s$n.censor[ ind[[i]] ]),
        group = rep(groups[i], n[i] + 1))
    }
    
    dat <- do.call(rbind, gr.df)
    dat.cens <- subset(dat, cens != 0)
    
    pl <- ggplot(dat, aes(x = time, y = surv, group = group)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(aes(col = group, lty = group))
    
    col <- if(length(surv.col == 1)){
      scale_colour_manual(name = gr.name, values = rep(surv.col, strata))
    } else{
      scale_colour_manual(name = gr.name, values = surv.col)
    }
    
    pl <- if(surv.col[1] != 'gg.def'){
      pl + col
    } else {pl + scale_colour_discrete(name = gr.name)}
    
    line <- if(length(lty.est) == 1){
      scale_linetype_manual(name = gr.name, values = rep(lty.est, strata))
    } else {scale_linetype_manual(name = gr.name, values = lty.est)}
    
    pl <- pl + line
    
    pl <- if(CI == T) {
      if(length(surv.col) > 1 && length(lty.est) > 1){
        stop('Either surv.col or lty.est should be of length 1 in order
             to plot 95% CI with multiple strata')
      }else if((length(surv.col) > 1 | surv.col == 'gg.def')[1]){
        pl + geom_step(aes(y = up, color = group), lty = lty.ci) +
          geom_step(aes(y = low, color = group), lty = lty.ci)
      } else{pl +  geom_step(aes(y = up, lty = group), col = surv.col) +
          geom_step(aes(y = low,lty = group), col = surv.col)}
    } else {pl}
    
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  pl <- if(strata == 1) {ggsurv.s(s, CI , plot.cens, surv.col ,
                                  cens.col, lty.est, lty.ci,
                                  cens.shape, back.white, xlab,
                                  ylab, main)
  } else {ggsurv.m(s, CI, plot.cens, surv.col ,
                   cens.col, lty.est, lty.ci,
                   cens.shape, back.white, xlab,
                   ylab, main)}
  pl
}

##END OF SURVIVAL CURVE PLOTTING FUNCTION

##MAKE PREDICTIONS!!

##LOADING TEST DATA
test <- read.csv("test.csv",header=T,sep=',',na.strings=c(""," ","NA"))

#LET R READ DATES

test$Assess.Date<-as.character.Date(test$Assess.Date)
test$V39_1<-as.character.Date(test$V39_1)

test$Assess.Date<-as.Date(test$Assess.Date,format = '%m/%d/%Y')
test$V39_1<-as.Date(test$V39_1,format = '%m/%d/%Y')
test <- test[,-3]
test$risk = NULL
test$outcome.prob = NULL

##TEST PERIOD IN WEEKS
test.period = 8
set.seed(0827)
r <- matrix(rep(0),nrow(test),5)
train1<- trainF
coxmod <- vector("list",5)
cols <- vector("list",5)
c <- (split(c(5:119),sample(c(5:119),size=5,replace=FALSE)))
##DEVELOP A MODEL
for(i in 1:5)
{
  cols[[i]] <- c("age",colnames(train1[unlist(c[i])]))
  my.formula <- as.formula(paste( "Surv(start,stop,outcome)", '~', paste( cols[[i]], collapse=' + ' ) ))
  coxmod[[i]] <- coxph(my.formula,data=na.omit(train1))
  step <- step(coxmod[[i]],direction="backward",trace=F)
  cols[[i]] <- names(step$assign)
  my.formula <- as.formula(paste( "Surv(start,stop,outcome)", '~', paste( cols[[i]], collapse=' + ' ) ))
  coxmod[[i]] <- coxph(my.formula,data=na.omit(train1))
  ggsurv(survfit(coxmod[[i]],data=na.omit(train1)),main=paste(paste("Survival Curve - Model : Surv ~ "), paste(cols[[i]], collapse=' + ' ),sep='\n' ))
}

##PREDICICTIONS
for (i in 1:5)
{
  for (j in 1:nrow(test))
  {
    r[j,i] <- risk(coxmod[[i]],test[j,],test.period)
  }
}
test[,206] <- apply(r[,1:5],1,function (x) mean(x,na.rm=T))
for (k in 1:nrow(test))
{
  ifelse(test[k,206]>=0.5,test[k,207]<-"HIGH RISK",ifelse(test[k,206]>=0.25,test[k,207]<-"MEDIUM RISK",test[k,207]<-"LOW RISK"))    
}
    
##Export the findings to the WD
filename <- paste("./Risk_Prediction_for_Time=",test.period,"_Weeks_",format(Sys.time(), "%a%b%d%Y%H-%M-%S"),".csv",sep = "")
write.csv(test[,c(1:3,206,207)],file=filename,row.names = FALSE)
for (i in 1:5)
{
  mypath <- file.path("C:","Users","te282346","Desktop","celtic","CelticData29AUG16", paste("SurvPlot_",i,"_Period=",test.period,"_",format(Sys.time(), "%a%b%d%Y%H-%M-%S"),".jpeg", sep = ""))
  jpeg(file=mypath)
  ggsurv(survfit(coxmod[[i]]),main=paste(paste("Survival Curve - Model : Surv ~ "), paste(cols[[i]], collapse=' + ' ),sep='\n' ))
  dev.off()
}
