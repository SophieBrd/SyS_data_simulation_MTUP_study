############################################################################################
#                                                                                          #
#                                                                                          #
#                       MODIFIABLE TEMPORAL UNIT PROMBLEM - MTUP                           #
#                           Aggregation of time series                                     #  
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
## Declaration of functions needed to aggregate the time series in main.R file
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
require(surveillance)
require(lubridate)
require(zoo)
library(sqldf)

#function to calculate sum on 7-moving-day counts
rolling <- function(x){
  rollapplyr(x, width=7, FUN=sum, na.rm=T, fill=NA)
} 

# From the S4 object time_series


setGeneric(name = "aggreg_time_series",
           def = function(time_series, options){standardGeneric("aggreg_time_series")}
             #time_series: object with simulated time series to aggregate (filling lists 
                      #totalRollSum and totalSegweek)
            # options : vector with time segmentation type : "rollSum" for 7-moving-day counts
                                                        # "segWeek" weekly counts
           )





setMethod(
  f = "aggreg_time_series",
  signature = "time_series",
  deendition = function (time_series,options){
    
        for (i in c(1:length(options))){
          
          #verification if options are OK  
          if(options[i] %ni% c("rollSum", "segWeek")){ 
            messageErreur = paste("Erreur - option", options[i] ,"incorrect. Options are : 'rollSum' or 'segWeek' ", sep=" ")  
            return(print(messageErreur))
          }else{
            
            for (scenario in c(1:16)){
                #print(scenario)
              scen <- paste("scenario", scenario, sep="" )
             
               if (options[i]=="rollSum"){
                  
                #print("rollSum")
                 
                 
                  time_series@totalRollSum[[scen]] = rolling(time_series@total[[scenario]][,1:100]) 
                  #time_series@totalRollSum[[scen]] = time_series@totalRollSum[[scenario]][-(1:6),] #to delete NA from first week

                  }#end rollSum
                
                if (options[i]=="segWeek"){
                 
                  time_series@totalSegweek[[scen]]=data.frame(array(rep(0,time_series@nsim*364),dim=c(364,time_series@nsim)))
                  time_series@outbreaksWeek[[scen]]=data.frame(array(rep(0,time_series@nsim*364),dim=c(364,time_series@nsim)))
                  
                  for (sim in c(1:time_series@nsim)){
                   
                  
                    
                    temp <- data.frame(cut=rep(time_series@total[[scen]]$dates[seq(from = 1, to = 2548, by = 7)], each = 7), 
                                    cas=time_series@total[[scen]][,sim], outbrk=time_series@outbreaks[[scen]][,sim])
                    
                    temp <- sqldf("select sum(cas) as cas, sum(outbrk) as outbrk, cut from temp group by cut ")
                    
                    time_series@totalSegweek[[scen]][,sim]=temp$cas
                    time_series@totalSegweek[[scen]]$dates=temp$cut
                    time_series@outbreaksWeek[[scen]][,sim]=temp$outbrk
                    time_series@outbreaksWeek[[scen]]$dates=temp$cut
                    
                    
                    rm(temp)
                    
                    
                    
                  }#end simulation segweek
                
                  
                  
                }#end if segWeek
                
                
                
            
              
            } #end for scenario
            
          }#end else
          
      }#end for i options
    
    return(time_series)
    
  }#end definition
)#end setMethod





# -----------------
# End of file
# -----------------