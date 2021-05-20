############################################################################################
#                                                                                          #
#                                                                                          #
#                       MODIFIABLE TEMPORAL UNIT PROMBLEM - MTUP                           #
#                               Time series simulation                                     #  
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
## Declaration of functions needed to simulates time series baselines and outbreaks in main.R file
## 
##
## This code, modified by Sophie Brilleaud, was originally written
## by Angela Noufaily and her team. 
## 
## This is the R code used in the research article "Joint assessment of time series segmentation and statistical algorithms for the detection of unusual events in syndromic surveillance in a one health perspective", S. Brilleaud et al , 202X
## 
## ############################################################################

# Load packages
require(data.table)
require(dplyr)
require(tidyr)
require(surveillance)
require(lubridate)
require(zoo)

# "not in" function 

"%ni%" = Negate("%in%")


#==============
# 5 day systems
#==============

h1=function(N,k,k2,theta,beta,gama1,gama2,gama3,gama4,shift,shift2){
  t=1:N
  if(k==0 & k2==0){h1=theta+beta*t}
  else{
    if(k==0)
    {
      l=1:k2
      h1=rep(0,N)
      for(i in 1:N){
        h1[i]=theta+beta*(t[i]+shift)+sum(gama3*cos((2*pi*l*(t[i]+shift))/5)+gama4*sin((2*pi*l*(t[i]+shift))/5))
      }
    }
    else{
      j=1:k
      l=1:k2
      h1=rep(0,N)
      for(i in 1:N){
        h1[i]=theta+beta*(t[i]+shift)+sum(gama1*cos((2*pi*j*(t[i]+shift))/(52*5))+gama2*sin((2*pi*j*(t[i]+shift2))/(52*5)))+sum(gama3*cos((2*pi*l*(t[i]+shift))/5)+gama4*sin((2*pi*l*(t[i]+shift))/5))
      }
    }
  }
  h1
}

negbinNoise1=function(N,k,k2,theta,beta,gama1,gama2,gama3,gama4,phi,shift,shift2){
  mu <- exp(h1(N,k,k2,theta,beta,gama1,gama2,gama3,gama4,shift,shift2))
  if(phi==1){yi <- rpois(N,mu)}
  else{
    prob <- 1/phi 
    size <- mu/(phi-1) 
    yi <- rnbinom(N,size=size,prob=prob)
  }
  yi
}


outbreak5=function(currentday,weeklength,wtime,yi,interval,k,k2,theta,beta,gama1,gama2,gama3,gama4,shift,shift2,phi,numoutbk,peakoutbk,meanlog,sdlog, correction){
  # theta, beta, gama1 and gama2 are the parameters of the equation for mu in Section 3.1
  
  N=length(yi)
  t=1:N
  
  mu <- exp(h1(N,k,k2,theta,beta,gama1,gama2,gama3,gama4,shift,shift2))/correction
  s=sqrt(mu*phi)
  
  #wtime = (currentday-49*5+1):currentday         # current outbreaks
  
  # GENERATING OUTBREAKS
  
  # STARTING TIMES OF OUTBREAKS
  
  startoutbk <- sample(wtime, numoutbk, replace = FALSE)
  
  # OUTBREAK SIZE OF CASES
  
  sizeoutbk=rep(0,numoutbk)
  for(i in 1:numoutbk){
    set.seed(i)
    soutbk=1
    sou=1
    while(soutbk<2){
      set.seed(sou)
      soutbk=rpois(1,s[startoutbk[i]]*peakoutbk)
      sou=sou+1
    }
    sizeoutbk[i]=soutbk
  }
  
  # DISTRIBUTE THESE CASES OVER TIME USING LOGNORMAL
  
  outbreak=rep(0,2*N)
  for( j in 1:numoutbk){
    set.seed(j)
    outbk <-rlnorm(sizeoutbk[j], meanlog = meanlog, sdlog = sdlog) 
    #outbk <-rnorm(sizeoutbk[j], mean = meanlog2, sd = sdlog) 
    #h<- hist(outbk,breaks=seq(0,ceiling(max(outbk)),1),plot=FALSE)
    h<- hist(outbk,breaks=seq(0,ceiling(max(outbk)),interval),plot=FALSE)
    cases <- h$counts
    weight=rep(0,length(cases))
    duration<-startoutbk:(startoutbk+length(cases)-1)
    dayofweek<-duration%%5 # 0 is friday; 1 is monday; 2 is tuesday etc.
    for(i in 1:length(cases)){
      if(dayofweek[i]==0){weight[i]=1.1}
      if(dayofweek[i]==1){weight[i]=1.5}
      if(dayofweek[i]==2){weight[i]=1.1}
      if(dayofweek[i]==3){weight[i]=1}
      if(dayofweek[i]==4){weight[i]=1}
    }
    cases2 <- cases*weight
    for (l in 1:(length(cases2))){
      outbreak[startoutbk[j]+(l-1)]= cases2[l]+outbreak[startoutbk[j]+(l-1)]
    }# l loop 
  }# j loop
  
  #for(v in 1:(currentday-49*5)){if(outbreak[v]>0){outbreak[v]=0}}
  for(v in currentday:(currentday+100)){if(outbreak[v]>0){outbreak[v]=0}}
  outbreak=outbreak[1:N]
  
  # ADD NOISE AND OUTBREAKS 
  
  yitot=yi+outbreak
  result=list(yitot=yitot,outbreak=outbreak,startoutbk=startoutbk,sizeoutbk=sizeoutbk,sd=s,mean=mu)
  #return(result)
}

#==============
# 7 days systems 
#==============

h2=function(N,k,k2,theta,beta,gama1,gama2,gama3,gama4,shift){
  t=1:N
  if(k==0 & k2==0){h2=theta+beta*t}
  else{
    if(k==0)
    {
      l=1:k2
      h2=rep(0,N)
      for(i in 1:N){
        h2[i]=theta+beta*(t[i]+shift)+sum(gama3*cos((2*pi*l*(t[i]+shift))/7)+gama4*sin((2*pi*l*(t[i]+shift))/7))
      }
    }
    else{
      j=1:k
      l=1:k2
      h2=rep(0,N)
      for(i in 1:N){
        h2[i]=theta+beta*(t[i]+shift)+sum(gama1*cos((2*pi*j*(t[i]+shift))/(52*7))+gama2*sin((2*pi*j*(t[i]+shift))/(52*7)))+sum(gama3*cos((2*pi*l*(t[i]+shift))/7)+gama4*sin((2*pi*l*(t[i]+shift))/7))
      }
    }
  }
  h2
}


negbinNoise2=function(N,k,k2,theta,beta,gama1,gama2,gama3,gama4,phi,shift){
  mu <- exp(h2(N,k,k2,theta,beta,gama1,gama2,gama3,gama4,shift))
  if(phi==1){yi <- rpois(N,mu)}
  else{
    prob <- 1/phi 
    size <- mu/(phi-1) 
    yi <- rnbinom(N,size=size,prob=prob)
  }
  yi
}


outbreak7=function(currentday,weeklength,wtime,yi,interval,k,k2,theta,beta,gama1,gama2,gama3,gama4,shift,phi,numoutbk,peakoutbk,meanlog,sdlog, correction){
  # theta, beta, gama1 and gama2 are the parameters of the equation for mu in Section 3.1
  
  N=length(yi)
  t=1:N
  
  mu <- exp(h2(N,k,k2,theta,beta,gama1,gama2,gama3,gama4,shift))/correction
  s=sqrt(mu*phi)
  
  
  #wtime = (currentday-49*7+1):currentday         # current outbreaks
  # wtime = 350*1:7         # current outbreaks
  
  # GENERATING OUTBREAKS
  
  # STARTING TIMES OF OUTBREAKS
  
  startoutbk <- sample(wtime, numoutbk, replace = FALSE)
  
  # OUTBREAK SIZE OF CASES
  
  sizeoutbk=rep(0,numoutbk)
  for(i in 1:numoutbk){
    set.seed(i)
    soutbk=1
    sou=1
    while(soutbk<2){
      set.seed(sou)
      soutbk=rpois(1,s[startoutbk[i]]*peakoutbk)
      sou=sou+1
    }
    sizeoutbk[i]=soutbk
  }
  
  # DISTRIBUTE THESE CASES OVER TIME USING LOGNORMAL
  
  outbreak=rep(0,2*N)
  for( j in 1:numoutbk){
    set.seed(j)
    while(sum(outbreak[1:N])==0)   {
      
    outbk <-rlnorm(sizeoutbk[j], meanlog = meanlog, sdlog = sdlog) 
    #outbk <-rnorm(sizeoutbk[j], mean = meanlog2, sd = sdlog) 
    #h<- hist(outbk,breaks=seq(0,ceiling(max(outbk)),1),plot=FALSE)
    h<- hist(outbk,breaks=seq(0,ceiling(max(outbk)),interval),plot=FALSE)
    cases <- h$counts
    weight=rep(0,length(cases))
    duration<-startoutbk:(startoutbk+length(cases)-1)
    dayofweek<-duration%%7 # 0 is sunday; 1 is monday; 2 is tuesday etc.
    for(i in 1:length(cases)){
      if(dayofweek[i]==0){weight[i]=2}
      if(dayofweek[i]==1){weight[i]=1}
      if(dayofweek[i]==2){weight[i]=1}
      if(dayofweek[i]==3){weight[i]=1}
      if(dayofweek[i]==4){weight[i]=1}
      if(dayofweek[i]==5){weight[i]=1}
      if(dayofweek[i]==6){weight[i]=2}
    }
    cases2 <- cases*weight
    for (l in 1:(length(cases2))){
      outbreak[startoutbk[j]+(l-1)]= cases2[l]+outbreak[startoutbk[j]+(l-1)]
    }# l loop 
  }# j loop
}#while

#for(v in (currentday-49*7):currentday){if(outbreak[v]>0){outbreak[v]=0}}
#for(v in currentday:(currentday+100)){if(outbreak[v]>0){outbreak[v]=0}}


  outbreak=outbreak[1:N]
  
  # ADD NOISE AND OUTBREAKS 
  
  yitot=yi+outbreak
  result=list(yitot=yitot,outbreak=outbreak,startoutbk=startoutbk,sizeoutbk=sizeoutbk,sd=s,mean=mu)
  #return(result)
}

#==========================
# Specify the bank holidays
#==========================



bankholidays=read.csv("Bankholidays.csv")

#fix(bankholidays)
bankhols7=bankholidays$bankhol
bankhols7=as.numeric(bankhols7)
length(bankhols7)
#fix(bankhols7)

bankhols5=bankhols7[-seq(6,length(bankhols7),7)] 
bankhols5=bankhols5[-seq(6,length(bankhols5),6)]
bankhols5=as.numeric(bankhols5)
length(bankhols5)
#fix(bankhols5)

#=======================
# Class definition


setClass(Class = "time_series",
         representation = representation(
           
           type = "numeric", #outbreak's size (very small=2, small=3, medium=5, big=10  ) 
           nsim = "numeric", #number of simulation by scenario
           years = "numeric", #time series' size in years
           
           # level 1 : scenario, level 2 : simulations
           baseline = "list", 
           outbreaks = "list", 
           seasoutbreaks = "list", 
           
           total = "list", #no aggregation=daily counts
           totalRollSum = "list",  # 7-day-moving day counts
           totalSegweek = "list", # weekly counts
           outbreaksWeek="list",
           
           #level 1 : detection algorithm, level 2 : scenarii, level 3 : simulations
            # daily counts
           alarmall= "list",
           
           # 7-day-moving day counts
           alarmRollSumall = "list",
           
           # weekly counts
           alarmSegWeekall = "list",
           
           #level 1 : time segmentation type, level 2 : scenarii
           indic_perf = "list"
         ),
          prototype = prototype(
            
            type = 0,
            nsim = 0,
            years = 0,
            
            
            baseline = list(),
            outbreaks = list(),
            seasoutbreaks = list(),
            
            total = list(),
            totalRollSum=list(),
            totalSegweek=list(),
            outbreaksWeek = list(),

            alarmall= list(),
            alarmRollSumall = list(),
            alarmSegWeekall = list(),
            indic_perf = list()
          )
         
         
         )



#Definition method to complete object of class "time_series" :

# method declaration 

setGeneric(
  name = "simulation_series",
  def = function (time_series,outbreak_size, nsim, years){standardGeneric("simulation_series")}
  # time_series : object to complete
  # outbreak_size :int, value : 2, 3, 5 or 10
  # nsim : number of simulation per scenario 
  # years : lenght of time series in years
          )

#method definition

setMethod(
  f = "simulation_series",
  signature = "time_series",
  definition = function (time_series,outbreak_size, nsim, years){
    
    time_series@type = outbreak_size
    time_series@years = years
    time_series@nsim = nsim
    
    
    
    #parameters used for the 16 scenarii described in the research article "Joint assessment of time series segmentation and statistical algorithms for the detection of unusual events in syndromic surveillance in a one health perspective", S. Brilleaud et al , 202X
    
    theta=c(6, 0.5, 5.5, 2, 6, 1, 6, 3, 3, 5, 0.5, 9, 2, 0.05, 3, 6)
    beta=c(0, 0, 0, 0, 0, 0, 0.0001, 0, 0, 0, 0, 0, 0.0005, 0, 0, 0)
    gama1=c(0.2, 1.5, 0, 0, 0.3, 0.1, 0, 1.5, 0.2, 0.2, 0.4, 0.5, 0.8, 0.01, 0.8, 0)
    gama2=c(0.2, 1.4 ,0, 0, 2, 2, 0, 0.1, 0.1, 0.1, 0, 0.2, 0.8, 0.01, 0.6, 0)
    gama3=c(0.5, 0.5, 0.3, 0.3, 0.3, 0.05, 0.6, 0.2, 0.05, 0.05, 0.05, 0.2, 0.8, 1.8, 0.8, 0.8)
    gama4=c(0.4, 0.4, 0.25, 0.25, 0.5, 0.05, 0.9, 0.3, 0.15, 0.1, 0.15, 0.5, 0.4, 0.1, 0.4, 0.4)
    phi=c(2, 1, 1, 1, 1.5, 1, 1.5, 1, 1, 1, 1, 1, 4, 1, 4, 4)
    shift=c(29, -167, 1, 1, -50, -50, 0, -150, -200, 0, 0, 0, 57, -85, 29, 1)
    shift2=c(29, -167, 1, 1, -50, -50, 0, -150, -200, 0, 0, 0, 57, -85, 29, 1) #for seasonal outbreaks
    k=c(1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 2, 1, 1, 4, 1, 0)
    k2=c(2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2)
    
    
    
    for (i in c(1:16)) { #for each scenario
      
      
      scen <- paste("scenario", i, sep="" )

      #table to save the simulations
      
      time_series@baseline[[scen]]=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
      time_series@outbreaks[[scen]]=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
      time_series@seasoutbreaks[[scen]]=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
      time_series@total[[scen]]=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
      
      
      
      
      
      #theoretical dates :
      myDates   <- seq(ymd('2010-01-01'), ymd('2016-12-30'), by = '1 day')
      
      dropDays <- as.Date(c('2010-12-31','2011-12-31', '2012-12-31',
                            '2013-12-31', '2014-12-31', '2015-12-31', 
                            '2016-02-29,', '2012-02-29'))
      myDates <- myDates[!(myDates %in% dropDays)]
      
      time_series@baseline[[scen]]["dates"]<-myDates
      time_series@outbreaks[[scen]]["dates"]<-myDates
      time_series@seasoutbreaks[[scen]]["dates"]<-myDates
      time_series@total[[scen]]["dates"]<-myDates

      
      if (i %in% c(5:12)){ # 5 days systems
        
        N<-52*5*years
        weekend<-seq(5,5*years*52,5)
        zeros<-rep(0,2)
        
        for (sim in c(1:nsim)){#for each simulation
          paste("sim=", sim)
     
          set.seed(sim)
          yt<-negbinNoise1(N=N,k=k[i],k2=k2[i],theta=theta[i],beta=beta[i],gama1=gama1[i],
                          gama2=gama2[i],gama3=gama3[i],gama4=gama4[i],phi=phi[i],shift=shift[i],shift2=shift2[i])
          #mu=exp(h1(N=N,k=k[i],k2=k2[i],theta=theta[i],beta=beta[i], gama1=gama1[i],gama2=gama2[i],gama3=gama3[i],gama4=gama4[i],phi=phi[i],shift=shift[i]))
          
          
          if (i %in% c(5,7)){ # cases number/10 for scenarii 5 and 7
            yt <- yt/10
          }#end if
          if (i == 12){ #cases number/100 for scenario 12
              yt <- yt/100
            }#end if
       
          
          set.seed(sim)
          
          
          if (i %in% c(5,7)){ # outbreaks/10 for scenarii 5 et 7
            out2<-outbreak5(currentday=5*52*years,weeklength=52*5*years,wtime=(length(yt)-49*5+1):length(yt),interval=0.25,
                           yi=yt,k=k[i],k2=k2[i],theta=theta[i],beta=beta[i],
                           gama1=gama1[i],gama2=gama2[i],gama3=gama3[i],gama4=gama4[i],phi=phi[i],shift=shift[i], 
                           shift2=shift2[i],numoutbk=1,peakoutbk=outbreak_size*5,meanlog=0,sdlog=0.5, correction = 10)
          }#end if
          if (i == 12){ #outbreaks/100 for scenario 12
            out2<-outbreak5(currentday=5*52*years,weeklength=52*5*years,wtime=(length(yt)-49*5+1):length(yt),interval=0.25,
                           yi=yt,k=k[i],k2=k2[i],theta=theta[i],beta=beta[i],
                           gama1=gama1[i],gama2=gama2[i],gama3=gama3[i],gama4=gama4[i],phi=phi[i],shift=shift[i], 
                           shift2=shift2[i],numoutbk=1,peakoutbk=outbreak_size*5,meanlog=0,sdlog=0.5, correction = 100)
            
          }#end if
          if (i %ni% c(5, 7, 12)){
            out2<-outbreak5(currentday=5*52*years,weeklength=52*5*years,wtime=(length(yt)-49*5+1):length(yt),interval=0.25,
                           yi=yt,k=k[i],k2=k2[i],theta=theta[i],beta=beta[i],
                           gama1=gama1[i],gama2=gama2[i],gama3=gama3[i],gama4=gama4[i],phi=phi[i],shift=shift[i], 
                           shift2=shift2[i],numoutbk=1,peakoutbk=outbreak_size*5,meanlog=0,sdlog=0.5, correction = 1)
          }
          
          
          zoutbreak<-out2$outbreak
          zt<-yt
          zitot<-yt + out2$outbreak 
          
          out1 <- rep(0,N)
          if (i == 5 ){ #seasonal outbreaks
            
            
            for(j in 1:years){
              set.seed(j+years*sim)
              out <- outbreak5(currentday=5*52*years,weeklength=52*5*years,wtime=((1+(j-1)*5*52):(20+(j-1)*5*52)),yi=yt,interval=0.02,k=1,
                            k2=1,theta=6,beta=0,gama1=0.3,gama2=2,gama3=0.3,gama4=0.5,phi=1.5,shift=-50,shift2=-50,numoutbk=1,peakoutbk=3*5*80,meanlog=0,sdlog=0.5, correction = 1)
              out1 <- out1+out$outbreak
            }#end for j
            zt <- zt + out1 
            zitot <- yt + out2$outbreak + out1
            zseasoutbreak <- out2$outbreak +out1
            
            for(s in 1:length(weekend)){
              zseasoutbreak <- append(zseasoutbreak,zeros,after=2*(s-1)+weekend[s])
            }#end for s
            
            
            time_series@seasoutbreaks[[scen]][,sim] <- round(zseasoutbreak)
            
          }#end if seasonal outbreaks
          
          
          
          if (i == 6){ ##seasonal outbreaks
            
            
            for(j in 1:years){
              set.seed(j+years*sim)
              out<-outbreak5(currentday=5*52*years,weeklength=52*5*years,wtime=((1+(j-1)*5*52):(20+(j-1)*5*52)),yi=yt,interval=0.02,k=1,k2=1,theta=1,beta=0,gama1=0.1,
                            gama2=2,gama3=0.05,gama4=0.05,phi=1,shift=-50,shift2=-50,numoutbk=1,peakoutbk=3*5*50,meanlog=0,sdlog=0.5, correction = 1)
              out1<-out1+out$outbreak
            }#end for j
            zt <- zt + out1 
            zitot <- yt + out2$outbreak + out1
            zseasoutbreak <- out2$outbreak + out1
            
            for(s in 1:length(weekend)){
              zseasoutbreak<- append(zseasoutbreak,zeros,after=2*(s-1)+weekend[s])
            }#end for s
            
            
            time_series@seasoutbreaks[[scen]][,sim]<-round(zseasoutbreak)
            
          }#end if seasonal outbreaks
          
          
          
          

          #zitot[(bankhols5==1)]=0
          #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]

          
          for(b in 1:length(zitot)){
            if(bankhols5[b]==1){
              zitot[b] <- 0
              zitot[b+1] <- 1.5*zitot[b+1]
              } 
            }#end for b
          
          
          #weekend=seq(0,days5*years*52-1,days5)
          for(s in 1:length(weekend)){
            zt <- append(zt,zeros,after=2*(s-1)+weekend[s])
            zitot <- append(zitot,zeros,after=2*(s-1)+weekend[s])
            zoutbreak <- append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
          }#end for s
          
          

          
          
          time_series@baseline[[scen]][,sim] <- round(zt)
          time_series@outbreaks[[scen]][,sim] <- round(zoutbreak)
          
          time_series@total[[scen]][,sim]<- round(zitot)
          
          
          
          }#end for simulation 5d
      
        }#end if for sim 5d
       
      else { # 7 days systems
         
         N <- 52*7*years
         
         for (sim in c(1:nsim)){
           
            
           set.seed(sim)
           yt <- negbinNoise2(N=N,k=k[i],k2=k2[i],theta=theta[i],beta=beta[i],
                           gama1=gama1[i],gama2=gama2[i],gama3=gama3[i],gama4=gama4[i],phi=phi[i],shift=shift[i])
           #mu=exp(h2(N=N,k=1,k2=2,alpha=6,beta=0,gama1=0.2,gama2=0.2,gama3=0.5,gama4=0.4,shift=29))
           
           set.seed(sim)
           
           # if (i == 13){
           #   out2=outbreak7(currentday=N,weeklength=52*7*years,wtime=(length(yt)-49*7+1):length(yt),yi=yt,interval=0.25,k=1,k2=2,theta=6,beta=0,gama1=0.2,gama2=0.2,
           #                  gama3=0.5,gama4=0.4,phi=2,shift=29,numoutbk=1,peakoutbk=outbreak_size*7,meanlog=0,sdlog=0.5, correction = 1)
           #   
           # }
           # if (i == 14){
           #   out2=outbreak7(currentday=N,weeklength=52*7*years,wtime=(length(yt)-49*7+1):length(yt),yi=yt,interval=0.25,k=1,k2=2,theta=2,beta=0,gama1=0.8,
           #                  gama2=0.8,gama3=0.8,gama4=0.4,phi=4,shift=57,numoutbk=1,peakoutbk=outbreak_size*7,meanlog=0,sdlog=0.5, correction = 1)
           # }
           # if (i== 15 ){
           #   out2=outbreak7(currentday=N,weeklength=52*7*years,wtime=(length(yt)-49*7+1):length(yt),yi=yt,interval=0.25,k=1,k2=2,theta=3,beta=0,gama1=0.8,
           #                  gama2=0.6,gama3=0.8,gama4=0.4,phi=4,shift=29,numoutbk=1,peakoutbk=outbreak_size*7,meanlog=0,sdlog=0.5, correction = 1)
           # }
           # if (i %ni% c(13, 14, 15)){
             out2<-outbreak7(currentday=N,weeklength=52*7*years,wtime=(length(yt)-49*7+1):length(yt),yi=yt,
                            interval=0.25,k=k[i],k2=k2[i],theta=theta[i],beta=beta[i],gama1=gama1[i],gama2=gama2[i],
                            gama3=gama3[i],gama4=gama4[i],phi=phi[i],shift=shift[i],numoutbk=1,peakoutbk=outbreak_size*7,
                            meanlog=0,sdlog=0.5, correction = 1)
           # }

           
           zoutbreak<-out2$outbreak
           
           
           out1<-rep(0,N)
           
           if (i == 15){ #seasonal outbreaks only for scenario 15
             for(j in 1:years){
               set.seed(j+years*sim)
               out<-outbreak5(currentday=7*52*years,weeklength=52*7*years,wtime=((210+(j-1)*7*52):(230+(j-1)*7*52)),yi=yt,interval=0.02,k=1,k2=1,theta=1,beta=0,gama1=0.1,
                             gama2=2,gama3=0.05,gama4=0.05,phi=1,shift=-50,shift2=-50,numoutbk=1,peakoutbk=3*7*150,meanlog=0,sdlog=0.5, correction = 1)
             out1<-out1+out$outbreak
             } #end boucle for j
             
             zoutbreak<-out2$outbreak
             zseasoutbreak<-out2$outbreak+out1
             
             zt<- yt+out1 
             zitot<-yt+out2$outbreak+out1
             time_series@seasoutbreaks[[scen]][,sim]=round(zseasoutbreak)
            }#end if scenario==15
           else{
             
             zt<-yt
             zitot<-yt + out2$outbreak
             
           }

           for(b in 1:length(zitot)){
             if(bankhols7[b]==1){
               zitot[b]<-2*zitot[b]
             } 
           }
           
           time_series@baseline[[scen]][,sim]<-round(zt)
           time_series@outbreaks[[scen]][,sim]<-round(zoutbreak)
           
           time_series@total[[scen]][,sim]<-round(zitot)
           
         }#end for simulation 7d
       }#end else for sim 7d 
      
      
      
      
      
      
      
    } #end for scenario
    
      return(time_series)
  }#end of function
)




# -----------------
# End of file
# -----------------




















