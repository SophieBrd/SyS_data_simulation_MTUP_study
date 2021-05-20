############################################################################################
#                                                                                          #
#                                                                                          #
#                       MODIFIABLE TEMPORAL UNIT PROMBLEM - MTUP                           #
#                                 Detection of unsual event                                #  
#                                                                                          #
#                                                                                          #
############################################################################################

## ############################################################################
##
## DISCLAIMER: 
## This script has been developed for research purposes only. 
## The script is provided without any warranty of any kind, either express or 
## implied. The entire risk arising out of the use or performance of the sample
## script and documentation remains with you. 
## In no event shall its author, or anyone else involved in the 
## creation, production, or delivery of the script be liable for any damages 
## whatsoever (including, without limitation, damages for loss of business 
## profits, business interruption, loss of business information, or other 
## pecuniary loss) arising out of the use of or inability to use the sample
## scripts or documentation, even if the author has been advised of the
## possibility of such damages. 
##
## ############################################################################
##
## DESCRIPTION
## Declaration of functions needed to detect unusual event the time series in main.R file
## 
##
## This code was written by Sophie Brilleaud
## 
## This is the R code used in the research article "Joint assessment of time series segmentation and statistical algorithms for the detection of unusual events in syndromic surveillance in a one health perspective", S. Brilleaud et al , 202X
## 
## ############################################################################

# Load packages
require(data.table)
require(dplyr)
require(tidyr)

require(lubridate)

require(surveillance)
require(zoo)

#============================================
# Four functions that implement the algorithm Farrington flexible
#============================================
# These functions were written by Angela Noufaily and her team. 
# Noufaily A, Morbey RA, Elliot AJ, Smith GE, Lake IR, McCarthy N. Comparison of statistical algorithms for daily syndromic surveillance aberration detection. 2019;35(17):3110-3118.



algo.farrington.N=function (disProgObj, control = list(range = NULL, b = 3, w = 3,
                                                     reweight = TRUE, verbose = FALSE, alpha = 0.01, trend = TRUE,
                                                     limit54 = c(5, 4), powertrans = "2/3", fitFun = "algo.farrington.fitGLM.N"))
{
  observed <- disProgObj$observed
  
  freq <- disProgObj$freq
  #epochStr <- switch(as.character(freq), `12` = "1 month",
  #    `52` = "1 week", `365` = "1 day")
  if (is.null(control$range)) {
    control$range <- (freq * control$b - control$w):length(observed)
  }
  if (is.null(control$b)) {
    control$b = 5
  }
  if (is.null(control$w)) {
    control$w = 3
  }
  if (is.null(control$reweight)) {
    control$reweight = TRUE
  }
  if (is.null(control$verbose)) {
    control$verbose = FALSE
  }
  if (is.null(control$alpha)) {
    control$alpha = 0.01
  }
  if (is.null(control$trend)) {
    control$trend = TRUE
  }
  if (is.null(control$plot)) {
    control$plot = FALSE
  }
  if (is.null(control$limit54)) {
    control$limit54 = c(5, 4)
  }
  if (is.null(control$powertrans)) {
    control$powertrans = "2/3"
  }
  if (is.null(control$fitFun)) {
    control$fitFun = "algo.farrington.fitGLM.N"
  }
  #else {
  #    control$fitFun <- match.arg(control$fitFun, c("algo.farrington.fitGLM"))
  #}
  #if (is.null(disProgObj[["epochAsDate", exact = TRUE]])) {
  #    epochAsDate <- FALSE
  #}
  #else {
  #    epochAsDate <- disProgObj[["epochAsDate", exact = TRUE]]
  #}
  if (!((control$limit54[1] >= 0) & (control$limit54[2] > 0))) {
    stop("The limit54 arguments are out of bounds: cases >= 0 and period > 0.")
  }
  # alarmall gives the results for all data sets whereas alarm identifies the
  # ones where there are less than 5 reports in the last 4 weeks
  alarmall <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  alarm<- matrix(data = 0, nrow = length(control$range), ncol = 1)
  upperbound <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  trend <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  pd <- matrix(data = 0, nrow = length(control$range), ncol = 2)
  
  n <- control$b * (2 * control$w + 1)
  
  # identify the leading zeroes in the data and cut them off
  
  #observed1=observed[observed!=NA]
  
  for (k in control$range) {
    #observed<-rollapply(data=observed[1:k], width=7, FUN=sum, ascending = F)
    #observed<-rollapply(data=a, width=7, FUN=sum, ascending = F)
    observed=observed[1:k]
    ob=rep(0,round(length(observed)/7))
    for(w in 0:(round(length(observed)/7)-1)){
      if((k-6*w-6-w)<=0){break}
      ob[w+1]=sum(observed[(k-6*w-w):(k-6*w-6-w)])
    }
    observed=ob[length(ob):1]
    kk=length(observed)
    for(z in 1:length(observed)){
      if(observed[z]==0){observed[z]=NA}
      else{break}
    }
    if (control$verbose) {
      cat("k=", k, "\n")
    }
    # If the remaining data is less than one year's worth,
    # do not fit a model and indicate that by reporting the value 1e+300
    if((kk-z)<(2*round(freq)+control$w+1)){
      upperbound[k - min(control$range) + 1] = 1e+300
      alarmall[k - min(control$range) + 1] = 1e+300
      alarm[k - min(control$range) + 1] = 1e+300
      trend[k - min(control$range) + 1] = 1e+300
    }
    else{
      #if (!epochAsDate) {
      # Assign level factors to the data (I have done that manually given 
      # that the data can be missing...)
      seas1 <- NULL
      seas2 <- NULL
      for (i in control$b:0) {
        if(i==3){seas1=append(seas1, seq(kk - round(freq * i) - 
                                           control$w, kk - round(freq * i) +control$w+1, by = 1))}
        else{seas1=append(seas1, seq(kk - round(freq * i) - 
                                       control$w, kk - round(freq * i) +control$w, by = 1))}
        seas2=append(seas2, seq(kk - round(freq * i) - 
                                  control$w-5, kk - round(freq * i) -control$w-1, by = 1))
      }
      #}
      seas3=seas2-5
      seas4=seas3-5
      seas5=seas4-5
      seas6=seas5-5
      seas7=seas6-5
      seas8=seas7-5
      seas9=seas8-5
      seas10=seas9-5
      seasgroup=rep(0,(length(observed)+control$w))
      seasgroup[seas1]=1
      seasgroup[seas2]=2
      seasgroup[seas3]=3
      seasgroup[seas4]=4
      seasgroup[seas5]=5
      seasgroup[seas6]=6
      seasgroup[seas7]=7
      seasgroup[seas8]=8
      seasgroup[seas9]=9
      seasgroup[seas10]=10
      if((kk-z)<(freq*control$b+control$w+1)){
        seasgroup=seasgroup[z:(kk-27)]
        seasgroup=as.factor(seasgroup)
        wtime= z:(kk-27) 
      }
      else{
        seasgroup=seasgroup[(kk-control$b*freq-control$w-1):(kk-27)]
        seasgroup=as.factor(seasgroup)
        wtime= (kk-control$b*freq-control$w-1):(kk-27) 
      }      
      response <- observed[wtime]   
      if (control$verbose) {
        print(response)
      }
      
      v=length(response[response>0])
      toosmall=(v<2) # indicator for when (stricktly) less than 2 weeks have non-zero counts
      
      # If stricktly less than 2 weeks have non-zero counts,
      # do not fit a model and indicate that by reporting the value 1e+400
      
      if(toosmall){
        upperbound[k - min(control$range) + 1] = 1e+400
        alarmall[k - min(control$range) + 1] = 1e+400
        alarm[k - min(control$range) + 1] = 1e+400
        trend[k - min(control$range) + 1] = 1e+400
      }
      else{
        oneyear=FALSE
        p=wtime[response>0]
        oneyear=((p[length(p)]-p[1])<52)
        # if all non-zero weeks are within one year, do not fit a trend.
        # This is the only case where we do not fit a trend.
        if(oneyear){
          model <- do.call(control$fitFun, args = list(response = response,
                                                       wtime = wtime, seasgroup=seasgroup,timeTrend = FALSE, reweight = control$reweight))
        }
        else{
          model <- do.call(control$fitFun, args = list(response = response,
                                                       wtime = wtime, seasgroup=seasgroup,timeTrend = control$trend, reweight = control$reweight))
        }
        if(is.null(model)){return(model)}
        
        doTrend <- control$trend
        if (model$T==1) {
          wp <- summary.glm(model)$coefficients["wtime", 4] # p-value after reweighting
          # In our new approach, we fit a trend always (unless all non-zero weeks are within one year).
          # For this reason, we consider the trend to be always significant (wp<=1).
          # If one wants to fit a trend only if significant, change this to, say, wp<=0.05.
          significant <- (wp <=1)
          mu0Hat <- predict.glm(model, data.frame(wtime = c(kk),seasgroup=factor(1)),type = "response")
          #atLeastThreeYears <- (control$b >= 3)
          # We remove the noExtrapolation condition
          #noExtrapolation <- mu0Hat <= max(response)
          if (!(significant)) {
            doTrend <- FALSE
            model <- do.call(control$fitFun, args = list(response = response,
                                                         wtime = wtime,seasgroup=seasgroup, timeTrend = FALSE, reweight = control$reweight))
          }
        }
        else {
          doTrend <- FALSE
        }
        if(model$phi<1){model$phi=1}
        pred <- predict.glm(model, data.frame(wtime = c(kk),seasgroup=factor(1)),
                            dispersion = model$phi, type = "response", se.fit = TRUE)
        # We use negative binomial quantiles to define the error structure rather 
        # than the ones based on the transformed Poisson and the Anscombe residuals.         
        if(model$phi==1){
          lu<-c(qpois(control$alpha,pred$fit),qpois(1-control$alpha,pred$fit))
          lu[2]=max(1,lu[2])
        }
        else{
          lu<-c(qnbinom(control$alpha,pred$fit/(model$phi-1),1/model$phi),qnbinom(1-control$alpha,pred$fit/(model$phi-1),1/model$phi))
          lu[2]=max(1,lu[2])
        }
        #if (control$plot) {
        #        data <- data.frame(wtime = seq(min(wtime), k, length = 1000))
        #        preds <- predict(model, data, type = "response",
        #            dispersion = model$phi)
        #        plot(c(wtime, k), c(response, observed[k]), ylim = range(c(observed[data$wtime],
        #            lu)), , xlab = "time", ylab = "No. infected",
        #            main = paste("Prediction at time t=", k, " with b=",
        #              control$b, ",w=", control$w, sep = ""), pch = c(rep(1,
        #              length(wtime)), 16))
        #        lines(data$wtime, preds, col = 1, pch = 2)
        #        lines(rep(k, 2), lu[1:2], col = 3, lty = 2)
        #    }
        
        enoughCases <- (sum(observed[(kk - control$limit54[2] + 1):kk]) >= control$limit54[1])
        
        X <- (observed[kk] - pred$fit)/(lu[2] - pred$fit)
        upperbound[k - min(control$range) + 1] <- lu[2]
        alarm[k - min(control$range) + 1] <- (X > 1)
        # if the last five weeks have less than 4 reports,
        # indicate that by returning the value 1e-300
        if(!enoughCases){alarm[k - min(control$range) + 1] = 1e-300}
        alarmall[k - min(control$range) + 1] <- (X > 1)
        trend[k - min(control$range) + 1] <- doTrend
      }
    }
    observed <- disProgObj$observed
  }
  control$name <- paste("farrington(", control$w, ",", 0, ",", control$b, ")", sep = "")
  control$data <- paste(deparse(substitute(disProgObj)))
  result <- list(alarm = alarm, alarmall = alarmall, trend = trend,
                 disProgObj = disProgObj, control = control,upperbound=upperbound)
  class(result) <- "survRes"
  return(result)
}

#=========================================
# (different from surveillance package)

algo.farrington.fitGLM.N=function (response, wtime,seasgroup, timeTrend = TRUE, reweight = TRUE)
{
  theModel <- as.formula(ifelse(timeTrend, "response~seasgroup+wtime",
                                "response~seasgroup"))
  model <- glm(theModel, family = quasipoisson(link = "log"))
  # In our new approach, we fit a trend always (unless all non-zero weeks are within one year).
  # For this reason, we consider the trend always significant (p<=1).
  # If one wants to fit a trend only if significant, change this to, say, p<=0.05.
  if (model$converged){
    if (timeTrend) {
      p<- summary.glm(model)$coefficients["wtime", 4] # p-value before reweighting
      if(p<=1){T=1}
      else{
        T=0
        model <- glm(response ~ seasgroup, family = quasipoisson(link = "log"))
        if(!model$converged) {
          cat("Warning: No convergence without insignificant timeTrend.\n")
          print(cbind(response, wtime))
          return(NULL)
        }
      }
    }
    if(!timeTrend){T=0}
  }
  
  else {
    T=0
    if (timeTrend) {
      model <- glm(response ~ seasgroup, family = quasipoisson(link = "log"))
      cat("Warning: No convergence with timeTrend -- trying without.\n")
    }
    if (!model$converged) {
      cat("Warning: No convergence without timeTrend.\n")
      print(cbind(response, wtime))
      return(NULL)
    }
  }
  
  phi <- max(summary(model)$dispersion, 1)
  if (reweight) {
    s <- anscombe.residuals(model, phi)
    omega <- algo.farrington.assign.weights.N(s)
    theModel <- as.formula(ifelse(T==1, "response~seasgroup+wtime", "response~seasgroup"))
    model <- glm(theModel, family = quasipoisson(link = "log"),weights = omega)
    if (!model$converged) {
      
      if (T==1) {
        model <- glm(response ~ seasgroup, family = quasipoisson(link = "log"),weights = omega)
        cat("Warning: No convergence with weights and timeTrend -- trying without.\n")
      }
      if (!model$converged) {
        cat("Warning: No convergence with weights but without timeTrend.\n")
        print(cbind(response, wtime))
        return(NULL)
      }
      T=0
    }
    
    phi <- max(summary(model)$dispersion, 1)
  }
  model$phi <- phi
  model$T=T
  return(model)
}

#=========================================

## (different from surveillance package) 

algo.farrington.assign.weights.N=function (s)
{   # we use the cut off point s>2.58 as a compromise between s>2 and s>3 
  gamma <- length(s)/(sum((s^(-2))^(s > 2.58)))
  omega <- numeric(length(s))
  omega[s > 2.58] <- gamma * (s[s > 2.58]^(-2))
  omega[s <= 2.58] <- gamma
  return(omega)
}


#=========================================


anscombe.residuals=function (m, phi) 
{
  y <- m$y
  mu <- fitted.values(m)
  a <- 3/2 * (y^(2/3) * mu^(-1/6) - mu^(1/2))
  a <- a/sqrt(phi * (1 - hatvalues(m)))
  return(a)
}


#==========================================






#### Definition methods for event detection ----------------

setGeneric(name = "detection_farrington",
           def = function(time_series){standardGeneric("detection_farrington")}
           #time_series: object with simulated time series (filling lists 
                        #alarm, alarmRollSum, alarmSegWeek if needed)
         
)


setMethod (f = "detection_farrington",
           signature = "time_series",
           definition = function (time_series){
            
                 #DAILY COUNTS
                   for (scenario in c(1:16)){
                     
                     scen <- paste("scenario", scenario, sep="" )
                     
                     
                     time_series@alarmall[["farrington"]][[scen]]=data.frame(array(rep(0,time_series@nsim*49*time_series@years),dim=c(49*time_series@years,time_series@nsim)))
                     for(sim in c(1:time_series@nsim)){
                        
                       a.disProg=create.disProg(week=weeks,observed=time_series@total[[scenario]][,sim],freq=52.18)
                        cntrl <- list(range=(nrow(time_series@total[[scenario]])-49*time_series@years+1):nrow(time_series@total[[scenario]]),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM.N")
                        a.farrington<- algo.farrington.N(a.disProg, control = cntrl)
                       
                        time_series@alarmall[["farrington"]][[scen]][,sim]=a.farrington$alarmall[,1]
                        
                        # alarmall gives the results for all data sets whereas alarm identifies the
                        # ones where there are less than 5 reports in the last 4 weeks
                        
                        
                        
                       
                      }# end for sim
                    }#end for scenario
                        
                   #7-MOVING-DAY COUNTS
                   if (length(time_series@totalRollSum)!=0){
                     for (scenario in c(1:16)){
                       
                       scen <- paste("scenario", scenario, sep="" )
                       
                       
                       time_series@alarmRollSumall[["farrington"]][[scen]]=data.frame(array(rep(0,time_series@nsim*49*time_series@years),dim=c(49*time_series@years,time_series@nsim)))
                       
                       for(sim in c(1:time_series@nsim)){
                         RollSum=time_series@totalRollSum[[scenario]][-(1:6),sim] #6 first values are NA
                         
                         a.disProg=create.disProg(week=weeks,observed=RollSum,freq=52.18)
                         cntrl <- list(range=(length(RollSum)-49*time_series@years+1):length(RollSum),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM.N")
                         a.farrington<- algo.farrington.N(a.disProg, control = cntrl)
                         time_series@alarmRollSumall[["farrington"]][[scen]][,sim]=a.farrington$alarmall[,1]
                         
                       }# end for sim
                     }#end for scenario
                   }
                   #WEEKLY COUNTS
                     if (length(time_series@totalSegweek)!=0){
                       
                       
                       
                       for (scenario in c(1:16)){
                         
                         scen <- paste("scenario", scenario, sep="" )
                         
                         time_series@alarmSegWeekall[["farrington"]][[scen]]=data.frame(array(rep(0,time_series@nsim*49),dim=c(49,time_series@nsim)))
                         
                         for(sim in c(1:time_series@nsim)){
                           
                           
                           
                           
                           #time_series@alarmSegWeekall[["farrington"]][[scen]]=data.frame(array(rep(0,time_series@nsim*49*time_series@years),dim=c(49*time_series@years,time_series@nsim)))
                           
                           a.disProg=create.disProg(week=49,observed=time_series@totalSegweek[[scen]][,sim],freq=52)
                           cntrl <- list(range=316:364,w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
                           a.farrington <- algo.farrington(a.disProg, control = cntrl)
                           
                           time_series@alarmSegWeekall[["farrington"]][[scen]][,sim]=a.farrington$alarm[,1]   

                         }# end for sim
                       }#end for scenario
                     
                    }#end if totalSegweek
                   
                
                 
                 

             return(time_series)
           
            
        } #fin definition 
 ) #fin setMethod
  




#############################################################################
# EARS C1

setGeneric(name = "detection_ears.c1",
           def = function(time_series){standardGeneric("detection_ears.c1")}
          
           
)


setMethod (f = "detection_ears.c1",
           signature = "time_series",
           definition = function (time_series){

 #DAILY COUNTS
             
             in2016 <- 2206:2548
             
             control <- list(range=in2016, method="C1", alpha=0.01)
             
             for (scenario in c(1:16)){
               
               
               scen <- paste("scenario", scenario, sep="" )
               
               time_series@alarmall[["ears.c1"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*343),dim=c(343,time_series@nsim)))
               
               # sts conversion
               totSts  <- sts(time_series@total[[scen]][1:100], start=c(2010,1), frequency=364,
                              epoch=as.numeric(as.Date(time_series@total[[scen]]$dates)), epochAsDate=TRUE)
               
               for(sim in c(1:time_series@nsim)){
                 
                 
                 
                 
                 det  <- earsC(totSts[,sim], control=control)
                 
                 
                 
                 # saving alarms 
                 time_series@alarmall[["ears.c1"]][[scen]][,sim] <- as.numeric(as.vector(unlist(det@alarm)))
                 
                 rm(det)
                 
               }#end for sim
               
               rm(totSts)
               
               
               time_series@alarmall[["ears.c1"]][[scen]][is.na(time_series@alarmall[["ears.c1"]][[scen]])] <- 0
               
               
             }#end for scenario
             
             
             
  # 7-MOVING DAY COUNTS
             if (length(time_series@totalRollSum)!=0){
             in2016 <- 2206:2548
             
             control <- list(range=in2016, method="C1", alpha=0.01)
             
             for (scenario in c(1:16)){
               
               
               scen <- paste("scenario", scenario, sep="" )
               
               time_series@alarmRollSumall[["ears.c1"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*343),dim=c(343,time_series@nsim)))
               
               # stsconversion
               
               totSts  <- sts(time_series@totalRollSum[[scen]], start=c(2010,1), frequency=364,
                              epoch=as.numeric(as.Date(time_series@total[[scen]]$dates)), epochAsDate=TRUE)
               
               for(sim in c(1:time_series@nsim)){
                 
                 
 
                 
                 det  <- earsC(totSts[,sim], control=control)
                 
                 
                 
                 # saving alarms 
                 time_series@alarmRollSumall[["ears.c1"]][[scen]][,sim] <- as.numeric(as.vector(unlist(det@alarm)))
                 
                 rm(det)
                 
               }#end for sim
               
               rm(totSts)  
               
               time_series@alarmRollSumall[["ears.c1"]][[scen]][is.na(time_series@alarmRollSumall[["ears.c1"]][[scen]])] <- 0
             }#end for scenario
             }#end if 7-moving-day
  
             #WEEKLY COUNTS
             
             if (length(time_series@totalSegweek)!=0){ 
             
             in2016 <- 316:364
             
             control <- list(range=in2016, method="C1", alpha=0.01)
             
             for (scenario in c(1:16)){
               
               
               scen <- paste("scenario", scenario, sep="" )
               
               time_series@alarmSegWeekall[["ears.c1"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49),dim=c(49,time_series@nsim)))
               
               # sts conversion
               totSts  <- sts(time_series@totalSegweek[[scen]][1:100], start=c(2010,1), frequency=49,
                              epoch=as.numeric(as.Date(time_series@totalSegweek[[scen]]$dates)), epochAsDate=TRUE)
               
               
               for(sim in c(1:time_series@nsim)){
                 
                 
 
                 
                 det  <- earsC(totSts[,sim], control=control)
                 
                 
                 
                 #saving alarms : 
                 time_series@alarmSegWeekall[["ears.c1"]][[scen]][,sim] <- as.numeric(as.vector(unlist(det@alarm)))
                 
                 rm(det)
                 
               }#end for sim
               
               
               rm(totSts)
               
               
               
               time_series@alarmSegWeekall[["ears.c1"]][[scen]][is.na(time_series@alarmSegWeekall[["ears.c1"]][[scen]])] <- 0
             }#end for scenario
             
             }#end if segmentation exist
             
             
             
             return(time_series)     
             
})# fin setmethod



############################################################

setGeneric(name = "detection_ears.c2",
           def = function(time_series){standardGeneric("detection_ears.c2")}

)


setMethod (f = "detection_ears.c2",
           signature = "time_series",
           definition = function (time_series){
            #DAILY COUNTS
             
             in2016 <- 2206:2548
             
             control <- list(range=in2016, method="C2", alpha=0.01)
             
             for (scenario in c(1:16)){
               
               
               scen <- paste("scenario", scenario, sep="" )
               
               time_series@alarmall[["ears.c2"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*343),dim=c(343,time_series@nsim)))
               
               # sts conversion
               totSts  <- sts(time_series@total[[scen]][1:100], start=c(2010,1), frequency=364,
                              epoch=as.numeric(as.Date(time_series@total[[scen]]$dates)), epochAsDate=TRUE)
               
               for(sim in c(1:time_series@nsim)){
                 
                 
 
                 
                 det  <- earsC(totSts[,sim], control=control)
                 
                 
                 
                 # saving alarms 
                 time_series@alarmall[["ears.c2"]][[scen]][,sim] <- as.numeric(as.vector(unlist(det@alarm)))
                 
                 rm(det)
                 
               }#end for sim
               
               rm(totSts)
               
               
               time_series@alarmall[["ears.c2"]][[scen]][is.na(time_series@alarmall[["ears.c2"]][[scen]])] <- 0
               
               
             }#end for scenario
             
             
             
             # 7-MOVING DAY COUNTS
             if (length(time_series@totalRollSum)!=0){
               in2016 <- 2206:2548
               
               control <- list(range=in2016, method="C2", alpha=0.01)
               
               for (scenario in c(1:16)){
                 
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 time_series@alarmRollSumall[["ears.c2"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*343),dim=c(343,time_series@nsim)))
                 
                 # sts conversion
                 
                 totSts  <- sts(time_series@totalRollSum[[scen]], start=c(2010,1), frequency=364,
                                epoch=as.numeric(as.Date(time_series@total[[scen]]$dates)), epochAsDate=TRUE)
                 
                 for(sim in c(1:time_series@nsim)){
                   
                   
   
                   
                   det  <- earsC(totSts[,sim], control=control)
                   
                   
                   
                   # saving alarms 
                   time_series@alarmRollSumall[["ears.c2"]][[scen]][,sim] <- as.numeric(as.vector(unlist(det@alarm)))
                   
                   rm(det)
                   
                 }#end for sim
                 
                 rm(totSts)  
                 
                 
                 
                 
                 
                 time_series@alarmRollSumall[["ears.c2"]][[scen]][is.na(time_series@alarmRollSumall[["ears.c2"]][[scen]])] <- 0
               }#end for scenario
             }#end if 7-moving-day
             
             #WEEKLY COUNTS
             
             if (length(time_series@totalSegweek)!=0){ 
               
               in2016 <- 316:364
               
               control <- list(range=in2016, method="C2", alpha=0.01)
               
               for (scenario in c(1:16)){
                 
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 time_series@alarmSegWeekall[["ears.c2"]][[scen]]  <- data.frame(array(rep(0,time_series@nsim*49),dim=c(49,time_series@nsim)))
                 
                 # sts conversion
                 totSts  <- sts(time_series@totalSegweek[[scen]][1:100], start=c(2010,1), frequency=49,
                                epoch=as.numeric(as.Date(time_series@totalSegweek[[scen]]$dates)), epochAsDate=TRUE)
                 
                 
                 for(sim in c(1:time_series@nsim)){
                   
                   
   
                   
                   det  <- earsC(totSts[,sim], control=control)
                   
                   
                   
                   # saving alarms 
                   time_series@alarmSegWeekall[["ears.c2"]][[scen]][,sim] <- as.numeric(as.vector(unlist(det@alarm)))
                   
                   rm(det)
                   
                 }#end for sim
                 
                 
                 rm(totSts)
                 
                 time_series@alarmSegWeekall[["ears.c2"]][[scen]][is.na(time_series@alarmSegWeekall[["ears.c2"]][[scen]])] <- 0
                 
                 
               }#end for scenario
               
             }#end if segementation existe
             
             
             
             return(time_series)     
             
           })# end setmethod


############################################################

setGeneric(name = "detection_ears.c3",
           def = function(time_series){standardGeneric("detection_ears.c3")}
           
           
)


setMethod (f = "detection_ears.c3",
           signature = "time_series",
           definition = function (time_series){
            
             #DAILY COUNTS
             
             in2016 <- 2206:2548
             
             control <- list(range=in2016, method="C3", alpha=0.01)
             
             for (scenario in c(1:16)){
               
               
               scen <- paste("scenario", scenario, sep="" )
               
               time_series@alarmall[["ears.c3"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*343),dim=c(343,time_series@nsim)))
               
               # sts conversion
               totSts  <- sts(time_series@total[[scen]][1:100], start=c(2010,1), frequency=364,
                              epoch=as.numeric(as.Date(time_series@total[[scen]]$dates)), epochAsDate=TRUE)
               
               for(sim in c(1:time_series@nsim)){
                 
                 
 
                 
                 det  <- earsC(totSts[,sim], control=control)
                 
                 
                 
                 # saving alarms 
                 time_series@alarmall[["ears.c3"]][[scen]][,sim] <- as.numeric(as.vector(unlist(det@alarm)))
                 
                 rm(det)
                 
               }#end for sim
               
               rm(totSts)
               
               
               time_series@alarmall[["ears.c3"]][[scen]][is.na(time_series@alarmall[["ears.c3"]][[scen]])] <- 0
               
               
               
             }#end for scenario
             
             
             
             # 7-MOVING DAY COUNTS
             if (length(time_series@totalRollSum)!=0){
               in2016 <- 2206:2548
               
               control <- list(range=in2016, method="C3", alpha=0.01)
               
               for (scenario in c(1:16)){
                 
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 time_series@alarmRollSumall[["ears.c3"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*343),dim=c(343,time_series@nsim)))
                 
                 # sts conversion
                 
                 totSts  <- sts(time_series@totalRollSum[[scen]], start=c(2010,1), frequency=364,
                                epoch=as.numeric(as.Date(time_series@total[[scen]]$dates)), epochAsDate=TRUE)
                 
                 for(sim in c(1:time_series@nsim)){
                   
                   
   
                   
                   det  <- earsC(totSts[,sim], control=control)
                   
                   
                   
                   # saving alarms 
                   time_series@alarmRollSumall[["ears.c3"]][[scen]][,sim] <- as.numeric(as.vector(unlist(det@alarm)))
                   
                   rm(det)
                   
                 }#end for sim
                 
                 rm(totSts) 
                 
                 
                 time_series@alarmRollSumall[["ears.c3"]][[scen]][is.na(time_series@alarmRollSumall[["ears.c3"]][[scen]])] <- 0
                 
                 
                 
               }#end for scenario
             }#end if 7-moving-day
             #WEEKLY COUNTS
             
             if (length(time_series@totalSegweek)!=0){ 
               
               in2016 <- 316:364
               
               control <- list(range=in2016, method="C3", alpha=0.01)
               
               for (scenario in c(1:16)){
                 
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 time_series@alarmSegWeekall[["ears.c3"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49),dim=c(49,time_series@nsim)))
                 
                 # sts conversion
                 totSts  <- sts(time_series@totalSegweek[[scen]][1:100], start=c(2010,1), frequency=49,
                                epoch=as.numeric(as.Date(time_series@totalSegweek[[scen]]$dates)), epochAsDate=TRUE)
                 
                 
                 for(sim in c(1:time_series@nsim)){
                   
                   
   
                   
                   det  <- earsC(totSts[,sim], control=control)
                   
                   
                   
                   # saving alarms 
                   time_series@alarmSegWeekall[["ears.c3"]][[scen]][,sim] <- as.numeric(as.vector(unlist(det@alarm)))
                   
                   rm(det)
                   
                 }#end for sim
                 
                 
                 rm(totSts)
  
                 time_series@alarmSegWeekall[["ears.c3"]][[scen]][is.na(time_series@alarmSegWeekall[["ears.c3"]][[scen]])] <- 0
                 
               }#end for scenario
               
             }#end if segementation existe
             
             
             
             return(time_series)     
             
           })# fin setmethod
################################################################################################################################

setGeneric(name = "detection_ears.NB",
           def = function(time_series){standardGeneric("detection_ears.NB")}
         
)


setMethod (f = "detection_ears.NB",
           signature = "time_series",
           definition = function (time_series){
            #DAILY COUNTS
             
             in2016 <- 2206:2548
             
             control <- list(range=in2016, c.ARL = 5, 
                             mu0 = NULL, alpha = 0, Mtilde = 1, M = -1, change = "intercept", 
                             theta = NULL)
             
             for (scenario in c(1:16)){
               
               
               scen <- paste("scenario", scenario, sep="" )
               
               time_series@alarmall[["ears.NB"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*343),dim=c(343,time_series@nsim)))
               
               # sts conversion
               totSts  <- sts(time_series@total[[scen]][1:100], start=c(2010,1), frequency=364,
                              epoch=as.numeric(as.Date(time_series@total[[scen]]$dates)), epochAsDate=TRUE)
               
               for(sim in c(1:time_series@nsim)){
                 
                 
 
                 
                 det  <- algo.glrnb(sts2disProg(totSts[,sim]), control=control)
                 
                 
                 # saving alarms 
                 time_series@alarmall[["ears.NB"]][[scen]][,sim] <-  as.numeric(as.vector(unlist(det$alarm)))
                 
                 rm(det)
                 
               }#end for sim
               
               rm(totSts)
               
               time_series@alarmall[["ears.NB"]][[scen]][is.na(time_series@alarmall[["ears.NB"]][[scen]])] <- 0
               
               
               
               
             }#end for scenario
             
             
             
             # 7-MOVING DAY COUNTS
             if (length(time_series@totalRollSum)!=0){
               in2016 <- 2206:2548
               
               control <- list(range=in2016, c.ARL = 5, 
                               mu0 = NULL, alpha = 0, Mtilde = 1, M = -1, change = "intercept", 
                               theta = NULL)
               
               for (scenario in c(1:16)){
                 
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 time_series@alarmRollSumall[["ears.NB"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*343),dim=c(343,time_series@nsim)))
                 
                 # sts conversion
                 
                 totSts  <- sts(time_series@totalRollSum[[scen]], start=c(2010,1), frequency=364,
                                epoch=as.numeric(as.Date(time_series@total[[scen]]$dates)), epochAsDate=TRUE)
                 
                 for(sim in c(1:time_series@nsim)){
                   
                   
   
                   
                   det  <- algo.glrnb(sts2disProg(totSts[,sim]), control=control)
                   
                   
                   
                   
                   # saving alarms 
                   time_series@alarmRollSumall[["ears.NB"]][[scen]][,sim] <-  as.numeric(as.vector(unlist(det$alarm)))
                   
                   rm(det)
                   
                 }#end for sim
                 
                 rm(totSts) 
                 
                 time_series@alarmRollSumall[["ears.NB"]][[scen]][is.na(time_series@alarmRollSumall[["ears.NB"]][[scen]])] <- 0
                 
                 
                 
               }#end for scenario
             }#end if 7-moving-day
             #WEEKLY COUNTS
             
             if (length(time_series@totalSegweek)!=0){ 
               
               in2016 <- 316:364
               
               control <- list(range=in2016, c.ARL = 5, 
                               mu0 = NULL, alpha = 0, Mtilde = 1, M = -1, change = "intercept", 
                               theta = NULL)
               
               for (scenario in c(1:16)){
                 
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 time_series@alarmSegWeekall[["ears.NB"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49),dim=c(49,time_series@nsim)))
                 
                 # sts conversion
                 totSts  <- sts(time_series@totalSegweek[[scen]][1:100], start=c(2010,1), frequency=49,
                                epoch=as.numeric(as.Date(time_series@totalSegweek[[scen]]$dates)), epochAsDate=TRUE)
                 
                 
                 for(sim in c(1:time_series@nsim)){
                   
                   
   
                   
                   det  <- algo.glrnb(sts2disProg(totSts[,sim]), control=control)
                   
                   
                   
                   # saving alarms 
                   time_series@alarmSegWeekall[["ears.NB"]][[scen]][,sim] <-  as.numeric(as.vector(unlist(det$alarm)))
                   
                   rm(det)
                   
                 }#end for sim
                 
                 
                 rm(totSts)
 
                 time_series@alarmSegWeekall[["ears.NB"]][[scen]][is.na(time_series@alarmSegWeekall[["ears.NB"]][[scen]])] <- 0
                 
                 
                                 
               }#end for scenario
               
             }#end if segementation existe
             
             
             
             return(time_series)     
             
           })# fin setmethod
##########################################################################################################

setGeneric(name = "detection_rki1",
           def = function(time_series){standardGeneric("detection_rki1")}
        
)


setMethod (f = "detection_rki1",
           signature = "time_series",
           definition = function (time_series){
            #DAILY COUNTS
             
             
             #DAILY COUNTS
             for (scenario in c(1:16)){
               
               scen <- paste("scenario", scenario, sep="" )
               
               
               time_series@alarmall[["rki1"]][[scen]]<-data.frame(array(rep(0,time_series@nsim*49*time_series@years),dim=c(49*time_series@years,time_series@nsim)))
               
               for(sim in c(1:time_series@nsim)){
                 
                 a.disProg <- create.disProg(week=weeks,observed=time_series@total[[scenario]][,sim],freq=364)
                 cntrl <- list(range=(nrow(time_series@total[[scenario]])-49*time_series@years+1):nrow(time_series@total[[scenario]]),w = (6*7), b = 0, alpha = 0.01, actY=TRUE)
                 a.rki<- algo.rki1(a.disProg, control = cntrl)
                 
                 time_series@alarmall[["rki1"]][[scen]][,sim] <- a.rki$alarm[,1]
                 
                 # alarmall gives the results for all data sets whereas alarm identifies the
                 # ones where there are less than 5 reports in the last 4 weeks
                 
                 
                 
                 
               }# end for sim
             }#end for scenario
             
             #7-MOVING-DAY COUNTS
             if (length(time_series@totalRollSum)!=0){
               for (scenario in c(1:16)){
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 
                 time_series@alarmRollSumall[["rki1"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49*time_series@years),dim=c(49*time_series@years,time_series@nsim)))
                 
                 for(sim in c(1:time_series@nsim)){
                   RollSum <- time_series@totalRollSum[[scenario]][-(1:6),sim] #les 6 premières valeurs sont NA
                   
                   a.disProg <- create.disProg(week=weeks,observed=RollSum,freq=364)
                   cntrl <- list(range=(length(RollSum)-49*time_series@years+1):length(RollSum),w = 6, b = 0, alpha = 0.01, actY=TRUE)
                   a.rki1<- algo.rki1(a.disProg, control = cntrl)
                   time_series@alarmRollSumall[["rki1"]][[scen]][,sim] <- a.rki1$alarm[,1]
                   
                 }# end for sim
               }#end for scenario
             }
             #WEEKLY COUNTS
             if (length(time_series@totalSegweek)!=0){
               
               
               
               for (scenario in c(1:16)){
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 time_series@alarmSegWeekall[["rki1"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49),dim=c(49,time_series@nsim)))
                 
                 for(sim in c(1:time_series@nsim)){
                   
                   
                   a.disProg <- create.disProg(week=49,observed=time_series@totalSegweek[[scen]][,sim],freq=52)
                   cntrl <- list(range=316:364,w = 6, b = 0, alpha = 0.01, actY=TRUE)
                   a.rki1 <- algo.rki1(a.disProg, control = cntrl)
                   
                   time_series@alarmSegWeekall[["rki1"]][[scen]][,sim] <- a.rki1$alarm[,1]   
                   
                 }# end for sim
               }#end for scenario
               
             }#end if totalSegweek
             
             
             
             
             
             return(time_series)
             
             
           } #fin definition 
) #fin setMethod
#################################################################################################################

setGeneric(name = "detection_rki2",
           def = function(time_series){standardGeneric("detection_rki2")}
         
)


setMethod (f = "detection_rki2",
           signature = "time_series",
           definition = function (time_series){
            #DAILY COUNTS
             
             
             #DAILY COUNTS
             for (scenario in c(1:16)){
               
               scen <- paste("scenario", scenario, sep="" )
               
               
               time_series@alarmall[["rki2"]][[scen]]<-data.frame(array(rep(0,time_series@nsim*49*time_series@years),dim=c(49*time_series@years,time_series@nsim)))
               
               for(sim in c(1:time_series@nsim)){
                 
                 a.disProg <- create.disProg(week=weeks,observed=time_series@total[[scenario]][,sim],freq=364)
                 cntrl <- list(range=(nrow(time_series@total[[scenario]])-49*time_series@years+1):nrow(time_series@total[[scenario]]),w = (6*7), b = 1, alpha = 0.01, actY=TRUE)
                 a.rki<- algo.rki2(a.disProg, control = cntrl)
                 
                 time_series@alarmall[["rki2"]][[scen]][,sim] <- a.rki$alarm[,1]
                 
                 # alarmall gives the results for all data sets whereas alarm identifies the
                 # ones where there are less than 5 reports in the last 4 weeks
                 
                 
                 
                 
               }# end for sim
             }#end for scenario
             
             #7-MOVING-DAY COUNTS
             if (length(time_series@totalRollSum)!=0){
               for (scenario in c(1:16)){
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 
                 time_series@alarmRollSumall[["rki2"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49*time_series@years),dim=c(49*time_series@years,time_series@nsim)))
                 
                 for(sim in c(1:time_series@nsim)){
                   RollSum <- time_series@totalRollSum[[scenario]][-(1:6),sim] #les 6 premières valeurs sont NA
                   
                   a.disProg <- create.disProg(week=weeks,observed=RollSum,freq=364)
                   cntrl <- list(range=(length(RollSum)-49*time_series@years+1):length(RollSum),w = 6, b = 1, alpha = 0.01, actY=TRUE)
                   a.rki2<- algo.rki2(a.disProg, control = cntrl)
                   time_series@alarmRollSumall[["rki2"]][[scen]][,sim] <- a.rki2$alarm[,1]
                   
                 }# end for sim
               }#end for scenario
             }
             #WEEKLY COUNTS
             if (length(time_series@totalSegweek)!=0){
               
               
               
               for (scenario in c(1:16)){
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 time_series@alarmSegWeekall[["rki2"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49),dim=c(49,time_series@nsim)))
                 
                 for(sim in c(1:time_series@nsim)){
                   
                   
                   a.disProg <- create.disProg(week=49,observed=time_series@totalSegweek[[scen]][,sim],freq=52)
                   cntrl <- list(range=316:364, w = 6, b = 1, alpha = 0.01, trend=TRUE, actY=TRUE)
                   a.rki2 <- algo.rki2(a.disProg, control = cntrl)
                   
                   time_series@alarmSegWeekall[["rki2"]][[scen]][,sim] <- a.rki2$alarm[,1]   
                   
                 }# end for sim
               }#end for scenario
               
             }#end if totalSegweek
             
             
             
             
             
             return(time_series)
             
             
           } #fin definition 
) #fin setMethod

###########################################################################################################

setGeneric(name = "detection_rki3",
           def = function(time_series){standardGeneric("detection_rki3")}
         
           
)


setMethod (f = "detection_rki3",
           signature = "time_series",
           definition = function (time_series){
            #DAILY COUNTS
             
             
             #DAILY COUNTS

             for (scenario in c(1:16)){
               
               scen <- paste("scenario", scenario, sep="" )
               
               
               time_series@alarmall[["rki3"]][[scen]]<-data.frame(array(rep(0,time_series@nsim*49*time_series@years),dim=c(49*time_series@years,time_series@nsim)))
               
               for(sim in c(1:time_series@nsim)){
                
                  a.disProg <- create.disProg(week=weeks,observed=time_series@total[[scenario]][,sim],freq=364)
                 cntrl <- list(range=(nrow(time_series@total[[scenario]])-49*time_series@years+1):nrow(time_series@total[[scenario]]),w = (4*7), b = 2, alpha = 0.01, actY=FALSE)
                 a.rki<- algo.rki3(a.disProg, control = cntrl)
                 
                 time_series@alarmall[["rki3"]][[scen]][,sim] <- a.rki$alarm[,1]
                 
                 # alarmall gives the results for all data sets whereas alarm identifies the
                 # ones where there are less than 5 reports in the last 4 weeks
                 
                 
                 
                 
               }# end for sim
             }#end for scenario
             
             #7-MOVING-DAY COUNTS
             if (length(time_series@totalRollSum)!=0){
               for (scenario in c(1:16)){
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 
                 time_series@alarmRollSumall[["rki3"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49*time_series@years),dim=c(49*time_series@years,time_series@nsim)))
                 
                 for(sim in c(1:time_series@nsim)){
                   RollSum <- time_series@totalRollSum[[scenario]][-(1:6),sim] #les 6 premières valeurs sont NA
                   
                   a.disProg <- create.disProg(week=weeks,observed=RollSum,freq=364)
                   cntrl <- list(range=(length(RollSum)-49*time_series@years+1):length(RollSum),w = 4, b = 2, alpha = 0.01, actY=FALSE)
                   a.rki3<- algo.rki3(a.disProg, control = cntrl)
                   time_series@alarmRollSumall[["rki3"]][[scen]][,sim] <- a.rki3$alarm[,1]
                   
                 }# end for sim
               }#end for scenario
             }
             #WEEKLY COUNTS
             if (length(time_series@totalSegweek)!=0){
               
               
               
               for (scenario in c(1:16)){
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 time_series@alarmSegWeekall[["rki3"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49),dim=c(49,time_series@nsim)))
                 
                 for(sim in c(1:time_series@nsim)){
                   
                   
                   a.disProg <- create.disProg(week=49,observed=time_series@totalSegweek[[scen]][,sim],freq=52)
                   cntrl <- list(range=316:364,w = 4, b = 2 , alpha = 0.01, actY=FALSE)
                   a.rki3 <- algo.rki3(a.disProg, control = cntrl)
                   
                   time_series@alarmSegWeekall[["rki3"]][[scen]][,sim] <- a.rki3$alarm[,1]   
                   
                 }# end for sim
               }#end for scenario
               
             }#end if totalSegweek
             
             
             
             
             
             return(time_series)
             
             
           } #fin definition 
) #fin setMethod

################################################################################################################################


setGeneric(name = "detection_bayes1",
           def = function(time_series){standardGeneric("detection_bayes1")}
          
)

 
setMethod (f = "detection_bayes1",
           signature = "time_series",
           definition = function (time_series){
            
             
             
             #DAILY COUNTS
             for (scenario in c(1:16)){
               
               scen <- paste("scenario", scenario, sep="" )
               
               
               time_series@alarmall[["bayes1"]][[scen]]<-data.frame(array(rep(0,time_series@nsim*49*7),dim=c(49*7,time_series@nsim)))
               
               for(sim in c(1:time_series@nsim)){
                 
                 a.disProg <- create.disProg(week=weeks,observed=time_series@total[[scenario]][,sim],freq=364)
                 cntrl <- list(range=(nrow(time_series@total[[scenario]]) - (49*7)+1):nrow(time_series@total[[scenario]]),w = 42, b = 0, alpha = 0.05, actY=TRUE)
                 #w=6*7=42 parce que 6 valeurs journalières pas assez (idem que pour rki)
                 
                 a.bayes <- algo.bayes(a.disProg, control = cntrl)
                 
                 time_series@alarmall[["bayes1"]][[scen]][,sim] <- a.bayes$alarm[,1]
                 
                 # alarmall gives the results for all data sets whereas alarm identifies the
                 # ones where there are less than 5 reports in the last 4 weeks
                 
                 
                 
                 
               }# end for sim
             }#end for scenario
             
             #7-MOVING-DAY COUNTS
             if (length(time_series@totalRollSum)!=0){
               for (scenario in c(1:16)){
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 
                 time_series@alarmRollSumall[["bayes1"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49*time_series@years),dim=c(49*time_series@years,time_series@nsim)))
                 
                 for(sim in c(1:time_series@nsim)){
                   RollSum <- time_series@totalRollSum[[scenario]][-(1:6),sim] #les 6 premières valeurs sont NA
                   
                   a.disProg <- create.disProg(week=weeks,observed=RollSum,freq=364)
                   cntrl <- list(range=(length(RollSum)-49*time_series@years+1):length(RollSum),w = 6, b = 0, alpha = 0.05, actY=TRUE)
                   
                   a.bayes <- algo.bayes(a.disProg, control = cntrl)
                   
                   time_series@alarmRollSumall[["bayes1"]][[scen]][,sim] <- a.bayes$alarm[,1]
                   
                 }# end for sim
               }#end for scenario
             }
             #WEEKLY COUNTS
             if (length(time_series@totalSegweek)!=0){
               
               
               
               for (scenario in c(1:16)){
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 time_series@alarmSegWeekall[["bayes1"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49),dim=c(49,time_series@nsim)))
                 
                 for(sim in c(1:time_series@nsim)){
                   
                   
                   a.disProg <- create.disProg(week=49,observed=time_series@totalSegweek[[scen]][,sim],freq=52)
                   
                   cntrl <- list(range=316:364,w = 6, b = 0, alpha = 0.05, actY=TRUE)
                   
                   a.bayes <- algo.bayes(a.disProg, control = cntrl)
                   
                   time_series@alarmSegWeekall[["bayes1"]][[scen]][,sim] <- a.bayes$alarm[,1]   
                   
                 }# end for sim
               }#end for scenario
               
             }#end if totalSegweek
             
             
             
             
             
             return(time_series)
             
             
           } #fin definition 
) #fin setMethod

#######################################################################################################################################

setGeneric(name = "detection_bayes2",
           def = function(time_series){standardGeneric("detection_bayes2")}
     
)


setMethod (f = "detection_bayes2",
           signature = "time_series",
           definition = function (time_series){
            #DAILY COUNTS
             
             
             #DAILY COUNTS
             for (scenario in c(1:16)){
               
               scen <- paste("scenario", scenario, sep="" )
               
               
               time_series@alarmall[["bayes2"]][[scen]]<-data.frame(array(rep(0,time_series@nsim*49*time_series@years),dim=c(49*time_series@years,time_series@nsim)))
               
               for(sim in c(1:time_series@nsim)){
                 
                 a.disProg <- create.disProg(week=weeks,observed=time_series@total[[scenario]][,sim],freq=364)
                 cntrl <- list(range=(nrow(time_series@total[[scenario]])-49*time_series@years+1):nrow(time_series@total[[scenario]]),w = (6*7), b = 1, alpha = 0.05, actY=TRUE)
                 
                 a.bayes <- algo.bayes(a.disProg, control = cntrl)
                 
                 time_series@alarmall[["bayes2"]][[scen]][,sim] <- a.bayes$alarm[,1]
                 
                 # alarmall gives the results for all data sets whereas alarm identifies the
                 # ones where there are less than 5 reports in the last 4 weeks
                 
                 
                 
                 
               }# end for sim
             }#end for scenario
             
             #7-MOVING-DAY COUNTS
             if (length(time_series@totalRollSum)!=0){
               for (scenario in c(1:16)){
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 
                 time_series@alarmRollSumall[["bayes2"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49*time_series@years),dim=c(49*time_series@years,time_series@nsim)))
                 
                 for(sim in c(1:time_series@nsim)){
                   RollSum <- time_series@totalRollSum[[scenario]][-(1:6),sim] #les 6 premières valeurs sont NA
                   
                   a.disProg <- create.disProg(week=weeks,observed=RollSum,freq=364)
                   cntrl <- list(range=(length(RollSum)-49*time_series@years+1):length(RollSum),w = 6, b = 1, alpha = 0.05, actY=TRUE)
                   
                   a.bayes <- algo.bayes(a.disProg, control = cntrl)
                   
                   time_series@alarmRollSumall[["bayes2"]][[scen]][,sim] <- a.bayes$alarm[,1]
                   
                 }# end for sim
               }#end for scenario
             }
             #WEEKLY COUNTS
             if (length(time_series@totalSegweek)!=0){
               
               
               
               for (scenario in c(1:16)){
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 time_series@alarmSegWeekall[["bayes2"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49),dim=c(49,time_series@nsim)))
                 
                 for(sim in c(1:time_series@nsim)){
                   
                   
                   a.disProg <- create.disProg(week=49,observed=time_series@totalSegweek[[scen]][,sim],freq=52)
                   cntrl <- list(range=316:364,w = 6, b = 1, alpha = 0.05, actY=TRUE)
                   a.bayes2 <- algo.bayes(a.disProg, control = cntrl)
                   
                   time_series@alarmSegWeekall[["bayes2"]][[scen]][,sim] <- a.bayes2$alarm[,1]   
                   
                 }# end for sim
               }#end for scenario
               
             }#end if totalSegweek
             
             
             
             
             
             return(time_series)
             
             
           } #fin definition 
) #fin setMethod

##########################################################################################################################

setGeneric(name = "detection_bayes3",
           def = function(time_series){standardGeneric("detection_bayes3")}

)


setMethod (f = "detection_bayes3",
           signature = "time_series",
           definition = function (time_series){
            #DAILY COUNTS
             
             
             #DAILY COUNTS
             for (scenario in c(1:16)){
               
               scen <- paste("scenario", scenario, sep="" )
               
               
               time_series@alarmall[["bayes3"]][[scen]]<-data.frame(array(rep(0,time_series@nsim*49*time_series@years),dim=c(49*time_series@years,time_series@nsim)))
               
               for(sim in c(1:time_series@nsim)){
                 
                 a.disProg <- create.disProg(week=weeks,observed=time_series@total[[scenario]][,sim],freq=364)
                 cntrl <- list(range=(nrow(time_series@total[[scenario]])-49*time_series@years+1):nrow(time_series@total[[scenario]]),w = (4*7), b = 2, alpha = 0.05, actY = FALSE)
                 
                 a.bayes <- algo.bayes(a.disProg, control = cntrl)
                 
                 time_series@alarmall[["bayes3"]][[scen]][,sim] <- a.bayes$alarm[,1]
                 
                 # alarmall gives the results for all data sets whereas alarm identifies the
                 # ones where there are less than 5 reports in the last 4 weeks
                 
                 
                 
                 
               }# end for sim
             }#end for scenario
             
             #7-MOVING-DAY COUNTS
             if (length(time_series@totalRollSum)!=0){
               for (scenario in c(1:16)){
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 
                 time_series@alarmRollSumall[["bayes3"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49*time_series@years),dim=c(49*time_series@years,time_series@nsim)))
                 
                 for(sim in c(1:time_series@nsim)){
                   RollSum <- time_series@totalRollSum[[scenario]][-(1:6),sim] #les 6 premières valeurs sont NA
                   
                   a.disProg <- create.disProg(week=weeks,observed=RollSum,freq=364)
                   cntrl <- list(range=(length(RollSum)-49*time_series@years+1):length(RollSum),w = 4, b = 2, alpha = 0.05, actY = FALSE)
                   
                   a.bayes <- algo.bayes(a.disProg, control = cntrl)
                   
                   time_series@alarmRollSumall[["bayes3"]][[scen]][,sim] <- a.bayes$alarm[,1]
                   
                 }# end for sim
               }#end for scenario
             }
             #WEEKLY COUNTS
             if (length(time_series@totalSegweek)!=0){
               
               
               
               for (scenario in c(1:16)){
                 
                 scen <- paste("scenario", scenario, sep="" )
                 
                 time_series@alarmSegWeekall[["bayes3"]][[scen]] <- data.frame(array(rep(0,time_series@nsim*49),dim=c(49,time_series@nsim)))
                 
                 for(sim in c(1:time_series@nsim)){
                   
                   
                   a.disProg <- create.disProg(week=49,observed=time_series@totalSegweek[[scen]][,sim],freq=52)
                   cntrl <- list(range=316:364,w = 4, b = 2, alpha = 0.05, actY = FALSE)
                   a.bayes3 <- algo.bayes(a.disProg, control = cntrl)
                   
                   time_series@alarmSegWeekall[["bayes3"]][[scen]][,sim] <- a.bayes3$alarm[,1]   
                   
                 }# end for sim
               }#end for scenario
               
             }#end if totalSegweek
             
             
             
             
             
             return(time_series)
             
             
           } #end definition 
) #end setMethod


# -----------------
# End of file
# -----------------