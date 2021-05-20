############################################################################################
#                                                                                          #
#                                                                                          #
#                       MODIFIABLE TEMPORAL UNIT PROMBLEM - MTUP                           #
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
## Simulation and event detection
## 
##
## This code was written by Sophie Brilleaud
## 
## This is the R code used in the research article "Joint assessment of time series segmentation and statistical algorithms for the detection of unusual events in syndromic surveillance in a one health perspective", S. Brilleaud et al , 202X
## 
## ############################################################################


#############################
# LOADING FUNCTIONS
#############################

source("1_data_generator.R")

source("2_time_series_aggregation.R")

source("3_event_detection.R")



######################################
# DATA SIMULATION
######################################


#very small outbreaks
outbreak2 <- new(Class = "time_series")
outbreak2 <- simulation_series(outbreak2, 2, 100, 7)

# small outbreaks
outbreak3 <- new(Class = "time_series")
outbreak3 <- simulation_series(outbreak3, 3, 100, 7)

#medium outbreaks
outbreak5 <- new(Class = "time_series")
outbreak5 <- simulation_series(outbreak5, 5, 100, 7)

#big outbreaks
outbreak10 <- new(Class = "time_series")
outbreak10 <- simulation_series(outbreak10, 10, 100, 7)



#######################
# TIME SERIES AGGREGATION
#######################



outbreak2 <- decoup_series(outbreak2, c("rollSum","segWeek"))

outbreak3 <- decoup_series(outbreak3, c("rollSum","segWeek"))

outbreak5 <- decoup_series(outbreak5, c("rollSum","segWeek"))

outbreak10 <- decoup_series(outbreak10, c("rollSum","segWeek"))



#############
# FARRINGTON
#############



outbreak2 <- detection_farrington(outbreak2)

outbreak3 <- detection_farrington(outbreak3)

outbreak5 <- detection_farrington(outbreak5)

outbreak10 <- detection_farrington(outbreak10)




###########
# EARS
###########

outbreak2 <- detection_ears.c1(outbreak2)
outbreak3 <- detection_ears.c1(outbreak3)
outbreak5 <- detection_ears.c1(outbreak5)
outbreak10 <- detection_ears.c1(outbreak10)


outbreak2 <- detection_ears.c2(outbreak2)
outbreak3 <- detection_ears.c2(outbreak3)
outbreak5 <- detection_ears.c2(outbreak5)
outbreak10 <- detection_ears.c2(outbreak10)


outbreak2 <- detection_ears.c3(outbreak2)
outbreak3 <- detection_ears.c3(outbreak3)
outbreak5 <- detection_ears.c3(outbreak5)
outbreak10 <- detection_ears.c3(outbreak10)



#############
#  RKI
#############

outbreak2 <- detection_rki1(outbreak2)
outbreak3 <- detection_rki1(outbreak3)
outbreak5 <- detection_rki1(outbreak5)
outbreak10 <- detection_rki1(outbreak10)



outbreak2 <- detection_rki2(outbreak2)
outbreak3 <- detection_rki2(outbreak3)
outbreak5 <- detection_rki2(outbreak5)
outbreak10 <- detection_rki2(outbreak10)


outbreak2 <- detection_rki3(outbreak2)
outbreak3 <- detection_rki3(outbreak3)
outbreak5 <- detection_rki3(outbreak5)
outbreak10 <- detection_rki3(outbreak10)

#############
#  BAYES
#############

outbreak2 <- detection_bayes1(outbreak2)
outbreak3 <- detection_bayes1(outbreak3)
outbreak5 <- detection_bayes1(outbreak5)
outbreak10 <- detection_bayes1(outbreak10)


outbreak2 <- detection_bayes2(outbreak2)
outbreak3 <- detection_bayes2(outbreak3)
outbreak5 <- detection_bayes2(outbreak5)
outbreak10 <- detection_bayes2(outbreak10)


outbreak2 <- detection_bayes3(outbreak2)
outbreak3 <- detection_bayes3(outbreak3)
outbreak5 <- detection_bayes3(outbreak5)
outbreak10 <- detection_bayes3(outbreak10)



##############################
# PERFORMANCE MEASURES
##############################

# In this part of the code, we show how the performance measures were calculated for the EARS NB algorithm:
# 16 Scenarii parameters used in the researc article

theta.vecteur <- c(6, 0.5, 5.5, 2, 6, 1, 6, 3, 3, 5, 0.5, 9, 2, 0.05, 3, 6)
tendance.vecteur <- as.factor(c(0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0)) 
k1.vecteur <- as.factor(c(1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 2, 1, 1, 4, 1, 0))
k2.vecteur <- as.factor(c(2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2))
phi.vecteur <- c(2, 1, 1, 1, 1.5, 1, 1.5, 1, 1, 1, 1, 1, 4, 1, 4, 4)


#### DAILY COUNTS - outbreak 2 ####

##### SCENARIO 1 - outbreak 2 
# table for 100 simulation of one scenario
dataNB <- data.frame(detect = rep(2, 100),
                     VP = rep(999, 100), 
                     VN = rep(999, 100),
                     FP = rep(999, 100), 
                     FN = rep(999, 100),
                     exces = rep(9999, 100),
                     duree= rep(999, 100), retard= rep(999, 100), retard= rep(999, 100), 
                     taille.outbk = rep("2", 100),
                     
                     algorithme = rep("ears.NB", 100),
                     
                     time_seg = rep("daily", 100), 
                     
                     scenario = as.factor(rep(1, 100)),
                     theta = rep(6, 100), 
                     tendance = as.factor(rep(0, 100)),
                     k1 = as.factor(rep(1, 100)), 
                     k2 = as.factor(rep(2, 100)), 
                     phi = rep(2, 100)
                     
)
longueur_tab_alarme <- length(outbreak2@alarmall$ears.NB$scenario1$X1)
longueur_tab_outbk <- length(outbreak2@outbreaks$scenario1[,1])
for (i in 1:100){
  
  #detected/not detected:
  
  if (sum(outbreak2@alarmall$ears.NB$scenario1[i] == 1 & outbreak2@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)){
    dataNB$detect[i] <- 1
  }else {dataNB$detect[i] <- 0}
  
  dataNB$VP[i] <-  sum(outbreak2@alarmall$ears.NB$scenario1[i] == 1 & outbreak2@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  dataNB$VN[i] <- sum(outbreak2@alarmall$ears.NB$scenario1[i] == 0 & outbreak2@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  
  dataNB$FP[i] <- sum(outbreak2@alarmall$ears.NB$scenario1[i] == 1 & outbreak2@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  dataNB$FN[i] <- sum(outbreak2@alarmall$ears.NB$scenario1[i] == 0 & outbreak2@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  
  
  
  
  dataNB$exces[i] <- sum(outbreak2@outbreaks$scenario1[i])/sum(outbreak2@total$scenario1[outbreak2@outbreaks$scenario1[i]>0,i])
  
  dataNB$duree[i] <- sum(outbreak2@outbreaks$scenario1[i]>0)
}

# all of the other scenarios



for (i in 2:16){
  
  tab <- data.frame(detect = rep(2, 100), 
                    exces = rep(9999, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("2", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("daily", 100), 
                    
                    scenario = as.factor(rep(i, 100)),
                    theta = rep(theta.vecteur[i], 100), 
                    tendance = as.factor(rep(tendance.vecteur[i], 100)),
                    k1 = as.factor(rep(k1.vecteur[i], 100)), 
                    k2 = as.factor(rep(k2.vecteur[i], 100)), 
                    phi = rep(phi.vecteur[i], 100)
                    
  )
  
  longueur_tab_alarme <- length(outbreak2@alarmall$ears.NB[[i]]$X1)
  longueur_tab_outbk <- length(outbreak2@outbreaks[[i]][,1])
  for (j in 1:100){
    
    #detected/not detected:
    
    if (sum(outbreak2@alarmall$ears.NB[[i]][j] == 1 & outbreak2@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)){
      tab$detect[j] <- 1
    }else {tab$detect[j] <- 0}
    tab$VP[j] <-  sum(outbreak2@alarmall$ears.NB[[i]][j] == 1 & outbreak2@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    tab$VN[j] <- sum(outbreak2@alarmall$ears.NB[[i]][j] == 0 & outbreak2@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    
    tab$FP[j] <- sum(outbreak2@alarmall$ears.NB[[i]][j] == 1 & outbreak2@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    tab$FN[j] <- sum(outbreak2@alarmall$ears.NB[[i]][j] == 0 & outbreak2@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    
    tab$exces[j] <- sum(outbreak2@outbreaks[[i]][j])/sum(outbreak2@total[[i]][outbreak2@outbreaks[[i]][j]>0,j])
    
    tab$duree[j] <- sum(outbreak2@outbreaks[[i]][j]>0)
  }
  
  
  
  dataNB <- rbind(dataNB, tab)
  
  rm(tab)
  
}


#### 7-MOVING-DAY COUNTS - outbreak 2 ####

##### SCENARIO 1 - outbreak 2 
# table for 100 simulation of one scenario
data2 <- data.frame(detect = rep(2, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    exces = rep(9999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("2", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("somme glissante", 100), 
                    
                    scenario = as.factor(rep(1, 100)),
                    theta = rep(6, 100), 
                    tendance = as.factor(rep(0, 100)),
                    k1 = as.factor(rep(1, 100)), 
                    k2 = as.factor(rep(2, 100)), 
                    phi = rep(2, 100)
                    
)
longueur_tab_alarme <- length(outbreak2@alarmRollSumall$ears.NB$scenario1$X1)
longueur_tab_outbk <- length(outbreak2@outbreaks$scenario1[,1])
for (i in 1:100){
  
  #detected/not detected:
  
  if (sum(outbreak2@alarmRollSumall$ears.NB$scenario1[i] == 1 & outbreak2@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)){
    data2$detect[i] <- 1
  }else {data2$detect[i] <- 0}
  
  data2$VP[i] <- sum(outbreak2@alarmRollSumall$ears.NB$scenario1[i] == 1 & outbreak2@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  data2$VN[i] <- sum(outbreak2@alarmRollSumall$ears.NB$scenario1[i] == 0 & outbreak2@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  
  data2$FP[i] <- sum(outbreak2@alarmRollSumall$ears.NB$scenario1[i] == 1 & outbreak2@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  data2$FN[i] <- sum(outbreak2@alarmRollSumall$ears.NB$scenario1[i] == 0 & outbreak2@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  
  
  
  
  data2$exces[i] <- sum(outbreak2@outbreaks$scenario1[i])/sum(outbreak2@total$scenario1[outbreak2@outbreaks$scenario1[i]>0,i])
  
  data2$duree[i] <- sum(outbreak2@outbreaks$scenario1[i]>0)
}

# all of the other scenarios



for (i in 2:16){
  
  tab <- data.frame(detect = rep(2, 100), 
                    exces = rep(9999, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("2", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("somme glissante", 100), 
                    
                    scenario = as.factor(rep(i, 100)),
                    theta = rep(theta.vecteur[i], 100), 
                    tendance = as.factor(rep(tendance.vecteur[i], 100)),
                    k1 = as.factor(rep(k1.vecteur[i], 100)), 
                    k2 = as.factor(rep(k2.vecteur[i], 100)), 
                    phi = rep(phi.vecteur[i], 100)
                    
  )
  
  longueur_tab_alarme <- length(outbreak2@alarmRollSumall$ears.NB[[i]]$X1)
  longueur_tab_outbk <- length(outbreak2@outbreaks[[i]][,1])
  for (j in 1:100){
    
    #detected/not detected:
    
    if (sum(outbreak2@alarmRollSumall$ears.NB[[i]][j] == 1 & outbreak2@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)){
      tab$detect[j] <- 1
    }else {tab$detect[j] <- 0}
    tab$VP[j] <- sum(outbreak2@alarmRollSumall$ears.NB[[i]][j] == 1 & outbreak2@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    tab$VN[j] <- sum(outbreak2@alarmRollSumall$ears.NB[[i]][j] == 0 & outbreak2@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    
    tab$FP[j] <- sum(outbreak2@alarmRollSumall$ears.NB[[i]][j] == 1 & outbreak2@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    tab$FN[j] <- sum(outbreak2@alarmRollSumall$ears.NB[[i]][j] == 0 & outbreak2@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    
    tab$exces[j] <- sum(outbreak2@outbreaks[[i]][j])/sum(outbreak2@total[[i]][outbreak2@outbreaks[[i]][j]>0,j])
    
    tab$duree[j] <- sum(outbreak2@outbreaks[[i]][j]>0)
  }
  
  
  
  data2 <- rbind(data2, tab)
  
  rm(tab)
  
}




######## WEEKLY COUNTS - outbreak 2 #######

##### SCENARIO 1 - outbreak 2 
# table for 100 simulation of one scenario
data3 <- data.frame(detect = rep(2, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    exces = rep(9999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("2", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("segmentation", 100), 
                    
                    scenario = as.factor(rep(1, 100)),
                    theta = rep(6, 100), 
                    tendance = as.factor(rep(0, 100)),
                    k1 = as.factor(rep(1, 100)), 
                    k2 = as.factor(rep(2, 100)), 
                    phi = rep(2, 100)
                    
)
longueur_tab_alarme <- length(outbreak2@alarmSegWeekall$ears.NB$scenario1$X1)
longueur_tab_outbk <- length(outbreak2@outbreaksWeek$scenario1[,1])
for (i in 1:100){
  
  #detected/not detected:
  
  if (sum(outbreak2@alarmSegWeekall$ears.NB$scenario1[i] == 1 & outbreak2@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)){
    data3$detect[i] <- 1
  }else {data3$detect[i] <- 0}
  
  data3$VP[i] <-  sum(outbreak2@alarmSegWeekall$ears.NB$scenario1[i] == 1 & outbreak2@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  data3$VN[i] <- sum(outbreak2@alarmSegWeekall$ears.NB$scenario1[i] == 0 & outbreak2@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  
  data3$FP[i] <- sum(outbreak2@alarmSegWeekall$ears.NB$scenario1[i] == 1 & outbreak2@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  data3$FN[i] <- sum(outbreak2@alarmSegWeekall$ears.NB$scenario1[i] == 0 & outbreak2@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  
  
  
  
  data3$exces[i] <- sum(outbreak2@outbreaksWeek$scenario1[i])/sum(outbreak2@totalSegweek$scenario1[outbreak2@outbreaksWeek$scenario1[i]>0,i])
  
  data3$duree[i] <- sum(outbreak2@outbreaks$scenario1[i]>0)
}

# all of the other scenarios



for (i in 2:16){
  
  tab <- data.frame(detect = rep(2, 100), 
                    exces = rep(9999, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("2", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("segmentation", 100), 
                    
                    scenario = as.factor(rep(i, 100)),
                    theta = rep(theta.vecteur[i], 100), 
                    tendance = as.factor(rep(tendance.vecteur[i], 100)),
                    k1 = as.factor(rep(k1.vecteur[i], 100)), 
                    k2 = as.factor(rep(k2.vecteur[i], 100)), 
                    phi = rep(phi.vecteur[i], 100)
                    
  )
  
  longueur_tab_alarme <- length(outbreak2@alarmSegWeekall$ears.NB[[i]]$X1)
  longueur_tab_outbk <- length(outbreak2@outbreaksWeek[[i]][,1])
  for (j in 1:100){
    
    #detected/not detected:
    
    if (sum(outbreak2@alarmSegWeekall$ears.NB[[i]][j] == 1 & outbreak2@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)){
      tab$detect[j] <- 1
    }else {tab$detect[j] <- 0}
    tab$VP[j] <- sum(outbreak2@alarmSegWeekall$ears.NB[[i]][j] == 1 & outbreak2@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    tab$VN[j] <- sum(outbreak2@alarmSegWeekall$ears.NB[[i]][j] == 0 & outbreak2@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    
    tab$FP[j] <- sum(outbreak2@alarmSegWeekall$ears.NB[[i]][j] == 1 & outbreak2@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    tab$FN[j] <- sum(outbreak2@alarmSegWeekall$ears.NB[[i]][j] == 0 & outbreak2@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    
    tab$exces[j] <- sum(outbreak2@outbreaksWeek[[i]][j])/sum(outbreak2@totalSegweek[[i]][outbreak2@outbreaksWeek[[i]][j]>0,j])
    
    tab$duree[j] <- sum(outbreak2@outbreaks[[i]][j]>0)
  }
  
  
  
  data3 <- rbind(data3, tab)
  
  rm(tab)
  
}



#### DAILY COUNTS - outbreak 3 ####

##### SCENARIO 1 - outbreak 3 
# table for 100 simulation of one scenario
data4 <- data.frame(detect = rep(2, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    exces = rep(9999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("3", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("daily", 100), 
                    
                    scenario = as.factor(rep(1, 100)),
                    theta = rep(6, 100), 
                    tendance = as.factor(rep(0, 100)),
                    k1 = as.factor(rep(1, 100)), 
                    k2 = as.factor(rep(2, 100)), 
                    phi = rep(2, 100)
                    
)
longueur_tab_alarme <- length(outbreak3@alarmall$ears.NB$scenario1$X1)
longueur_tab_outbk <- length(outbreak3@outbreaks$scenario1[,1])
for (i in 1:100){
  
  #detected/not detected:
  
  if (sum(outbreak3@alarmall$ears.NB$scenario1[i] == 1 & outbreak3@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)){
    data4$detect[i] <- 1
  }else {data4$detect[i] <- 0}
  
  data4$VP[i] <- sum(outbreak3@alarmall$ears.NB$scenario1[i] == 1 & outbreak3@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  data4$VN[i] <- sum(outbreak3@alarmall$ears.NB$scenario1[i] == 0 & outbreak3@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  
  data4$FP[i] <- sum(outbreak3@alarmall$ears.NB$scenario1[i] == 1 & outbreak3@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  data4$FN[i] <- sum(outbreak3@alarmall$ears.NB$scenario1[i] == 0 & outbreak3@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  
  
  
  
  data4$exces[i] <- sum(outbreak3@outbreaks$scenario1[i])/sum(outbreak3@total$scenario1[outbreak3@outbreaks$scenario1[i]>0,i])
  
  data4$duree[i] <- sum(outbreak3@outbreaks$scenario1[i]>0)
}

# all of the other scenarios



for (i in 2:16){
  
  tab <- data.frame(detect = rep(2, 100), 
                    exces = rep(9999, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("3", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("daily", 100), 
                    
                    scenario = as.factor(rep(i, 100)),
                    theta = rep(theta.vecteur[i], 100), 
                    tendance = as.factor(rep(tendance.vecteur[i], 100)),
                    k1 = as.factor(rep(k1.vecteur[i], 100)), 
                    k2 = as.factor(rep(k2.vecteur[i], 100)), 
                    phi = rep(phi.vecteur[i], 100)
                    
  )
  
  longueur_tab_alarme <- length(outbreak3@alarmall$ears.NB[[i]]$X1)
  longueur_tab_outbk <- length(outbreak3@outbreaks[[i]][,1])
  for (j in 1:100){
    
    #detected/not detected:
    
    if (sum(outbreak3@alarmall$ears.NB[[i]][j] == 1 & outbreak3@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)){
      tab$detect[j] <- 1
    }else {tab$detect[j] <- 0}
    tab$VP[j] <- sum(outbreak3@alarmall$ears.NB[[i]][j] == 1 & outbreak3@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    tab$VN[j] <- sum(outbreak3@alarmall$ears.NB[[i]][j] == 0 & outbreak3@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    
    tab$FP[j] <- sum(outbreak3@alarmall$ears.NB[[i]][j] == 1 & outbreak3@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    tab$FN[j] <- sum(outbreak3@alarmall$ears.NB[[i]][j] == 0 & outbreak3@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    
    tab$exces[j] <- sum(outbreak3@outbreaks[[i]][j])/sum(outbreak3@total[[i]][outbreak3@outbreaks[[i]][j]>0,j])
    
    tab$duree[j] <- sum(outbreak3@outbreaks[[i]][j]>0)
  }
  
  
  
  data4 <- rbind(data4, tab)
  
  rm(tab)
  
}


#### 7-MOVING-DAY COUNTS - outbreak 3 ####

##### SCENARIO 1 - outbreak 3 
# table for 100 simulation of one scenario
data5 <- data.frame(detect = rep(2, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    exces = rep(9999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("3", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("somme glissante", 100), 
                    
                    scenario = as.factor(rep(1, 100)),
                    theta = rep(6, 100), 
                    tendance = as.factor(rep(0, 100)),
                    k1 = as.factor(rep(1, 100)), 
                    k2 = as.factor(rep(2, 100)), 
                    phi = rep(2, 100)
                    
)
longueur_tab_alarme <- length(outbreak3@alarmRollSumall$ears.NB$scenario1$X1)
longueur_tab_outbk <- length(outbreak3@outbreaks$scenario1[,1])
for (i in 1:100){
  
  #detected/not detected:
  
  if (sum(outbreak3@alarmRollSumall$ears.NB$scenario1[i] == 1 & outbreak3@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)){
    data5$detect[i] <- 1
  }else {data5$detect[i] <- 0}
  
  data5$VP[i] <- sum(outbreak3@alarmRollSumall$ears.NB$scenario1[i] == 1 & outbreak3@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  data5$VN[i] <- sum(outbreak3@alarmRollSumall$ears.NB$scenario1[i] == 0 & outbreak3@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  
  data5$FP[i] <- sum(outbreak3@alarmRollSumall$ears.NB$scenario1[i] == 1 & outbreak3@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  data5$FN[i] <- sum(outbreak3@alarmRollSumall$ears.NB$scenario1[i] == 0 & outbreak3@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  
  
  
  
  data5$exces[i] <- sum(outbreak3@outbreaks$scenario1[i])/sum(outbreak3@total$scenario1[outbreak3@outbreaks$scenario[i]>0,i])
  if(is.infinite(data5$exces[i])){data5$exces[i] <- 0}
  
  
  
  data5$duree[i] <- sum(outbreak3@outbreaks$scenario1[i]>0)
}

# all of the other scenarios



for (i in 2:16){
  
  tab <- data.frame(detect = rep(2, 100), 
                    exces = rep(9999, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("3", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("somme glissante", 100), 
                    
                    scenario = as.factor(rep(i, 100)),
                    theta = rep(theta.vecteur[i], 100), 
                    tendance = as.factor(rep(tendance.vecteur[i], 100)),
                    k1 = as.factor(rep(k1.vecteur[i], 100)), 
                    k2 = as.factor(rep(k2.vecteur[i], 100)), 
                    phi = rep(phi.vecteur[i], 100)
                    
  )
  
  longueur_tab_alarme <- length(outbreak3@alarmRollSumall$ears.NB[[i]]$X1)
  longueur_tab_outbk <- length(outbreak3@outbreaks[[i]][,1])
  for (j in 1:100){
    
    #detected/not detected:
    
    if (sum(outbreak3@alarmRollSumall$ears.NB[[i]][j] == 1 & outbreak3@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)){
      tab$detect[j] <- 1
    }else {tab$detect[j] <- 0}
    tab$VP[j] <- sum(outbreak3@alarmRollSumall$ears.NB[[i]][j] == 1 & outbreak3@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    tab$VN[j] <- sum(outbreak3@alarmRollSumall$ears.NB[[i]][j] == 0 & outbreak3@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    
    tab$FP[j] <- sum(outbreak3@alarmRollSumall$ears.NB[[i]][j] == 1 & outbreak3@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    tab$FN[j] <- sum(outbreak3@alarmRollSumall$ears.NB[[i]][j] == 0 & outbreak3@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    
    tab$exces[j] <- sum(outbreak3@outbreaks[[i]][j])/sum(outbreak3@total[[i]][outbreak3@outbreaks[[i]][j]>0,j])
    
    tab$duree[j] <- sum(outbreak3@outbreaks[[i]][j]>0)
  }
  
  
  
  data5 <- rbind(data5, tab)
  
  rm(tab)
  
}




######## WEEKLY COUNTS - outbreak 3 #######

##### SCENARIO 1 - outbreak 3 
# table for 100 simulation of one scenario
data6 <- data.frame(detect = rep(2, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    exces = rep(9999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("3", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("segmentation", 100), 
                    
                    scenario = as.factor(rep(1, 100)),
                    theta = rep(6, 100), 
                    tendance = as.factor(rep(0, 100)),
                    k1 = as.factor(rep(1, 100)), 
                    k2 = as.factor(rep(2, 100)), 
                    phi = rep(2, 100)
                    
)
longueur_tab_alarme <- length(outbreak3@alarmSegWeekall$ears.NB$scenario1$X1)
longueur_tab_outbk <- length(outbreak3@outbreaksWeek$scenario1[,1])
for (i in 1:100){
  
  #detected/not detected:
  
  if (sum(outbreak3@alarmSegWeekall$ears.NB$scenario1[i] == 1 & outbreak3@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)){
    data6$detect[i] <- 1
  }else {data6$detect[i] <- 0}
  
  data6$VP[i] <- sum(outbreak3@alarmSegWeekall$ears.NB$scenario1[i] == 1 & outbreak3@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  data6$VN[i] <- sum(outbreak3@alarmSegWeekall$ears.NB$scenario1[i] == 0 & outbreak3@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  
  data6$FP[i] <- sum(outbreak3@alarmSegWeekall$ears.NB$scenario1[i] == 1 & outbreak3@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  data6$FN[i] <- sum(outbreak3@alarmSegWeekall$ears.NB$scenario1[i] == 0 & outbreak3@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  
  
  
  
  data6$exces[i] <- sum(outbreak3@outbreaksWeek$scenario1[i])/sum(outbreak3@totalSegweek$scenario1[outbreak3@outbreaksWeek$scenario1[i]>0,i])
  
  data6$duree[i] <- sum(outbreak3@outbreaks$scenario1[i]>0)
}

# all of the other scenarios



for (i in 2:16){
  
  tab <- data.frame(detect = rep(2, 100), 
                    exces = rep(9999, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("3", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("segmentation", 100), 
                    
                    scenario = as.factor(rep(i, 100)),
                    theta = rep(theta.vecteur[i], 100), 
                    tendance = as.factor(rep(tendance.vecteur[i], 100)),
                    k1 = as.factor(rep(k1.vecteur[i], 100)), 
                    k2 = as.factor(rep(k2.vecteur[i], 100)), 
                    phi = rep(phi.vecteur[i], 100)
                    
  )
  
  longueur_tab_alarme <- length(outbreak3@alarmSegWeekall$ears.NB[[i]]$X1)
  longueur_tab_outbk <- length(outbreak3@outbreaksWeek[[i]][,1])
  for (j in 1:100){
    
    #detected/not detected:
    
    if (sum(outbreak3@alarmSegWeekall$ears.NB[[i]][j] == 1 & outbreak3@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)){
      tab$detect[j] <- 1
    }else {tab$detect[j] <- 0}
    tab$VP[j] <- sum(outbreak3@alarmSegWeekall$ears.NB[[i]][j] == 1 & outbreak3@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    tab$VN[j] <- sum(outbreak3@alarmSegWeekall$ears.NB[[i]][j] == 0 & outbreak3@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    
    tab$FP[j] <- sum(outbreak3@alarmSegWeekall$ears.NB[[i]][j] == 1 & outbreak3@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    tab$FN[j] <- sum(outbreak3@alarmSegWeekall$ears.NB[[i]][j] == 0 & outbreak3@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    
    tab$exces[j] <- sum(outbreak3@outbreaksWeek[[i]][j])/sum(outbreak3@totalSegweek[[i]][outbreak3@outbreaksWeek[[i]][j]>0,j])
    
    tab$duree[j] <- sum(outbreak3@outbreaks[[i]][j]>0)
  }
  
  
  
  data6 <- rbind(data6, tab)
  
  rm(tab)
  
}


#### DAILY COUNTS - outbreak 5 ####

##### SCENARIO 1 - outbreak 5 
# table for 100 simulation of one scenario
data7 <- data.frame(detect = rep(2, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    exces = rep(9999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("5", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("daily", 100), 
                    
                    scenario = as.factor(rep(1, 100)),
                    theta = rep(6, 100), 
                    tendance = as.factor(rep(0, 100)),
                    k1 = as.factor(rep(1, 100)), 
                    k2 = as.factor(rep(2, 100)), 
                    phi = rep(2, 100)
                    
)
longueur_tab_alarme <- length(outbreak5@alarmall$ears.NB$scenario1$X1)
longueur_tab_outbk <- length(outbreak5@outbreaks$scenario1[,1])
for (i in 1:100){
  
  #detected/not detected:
  
  if (sum(outbreak5@alarmall$ears.NB$scenario1[i] == 1 & outbreak5@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)){
    data7$detect[i] <- 1
  }else {data7$detect[i] <- 0}
  
  data7$VP[i] <- sum(outbreak5@alarmall$ears.NB$scenario1[i] == 1 & outbreak5@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  data7$VN[i] <- sum(outbreak5@alarmall$ears.NB$scenario1[i] == 0 & outbreak5@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  
  data7$FP[i] <- sum(outbreak5@alarmall$ears.NB$scenario1[i] == 1 & outbreak5@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  data7$FN[i] <- sum(outbreak5@alarmall$ears.NB$scenario1[i] == 0 & outbreak5@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  
  
  
  
  data7$exces[i] <- sum(outbreak5@outbreaks$scenario1[i])/sum(outbreak5@total$scenario1[outbreak5@outbreaks$scenario1[i]>0,i])
  
  data7$duree[i] <- sum(outbreak5@outbreaks$scenario1[i]>0)
}

# all of the other scenarios



for (i in 2:16){
  
  tab <- data.frame(detect = rep(2, 100), 
                    exces = rep(9999, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("5", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("daily", 100), 
                    
                    scenario = as.factor(rep(i, 100)),
                    theta = rep(theta.vecteur[i], 100), 
                    tendance = as.factor(rep(tendance.vecteur[i], 100)),
                    k1 = as.factor(rep(k1.vecteur[i], 100)), 
                    k2 = as.factor(rep(k2.vecteur[i], 100)), 
                    phi = rep(phi.vecteur[i], 100)
                    
  )
  
  longueur_tab_alarme <- length(outbreak5@alarmall$ears.NB[[i]]$X1)
  longueur_tab_outbk <- length(outbreak5@outbreaks[[i]][,1])
  for (j in 1:100){
    
    #detected/not detected:
    
    if (sum(outbreak5@alarmall$ears.NB[[i]][j] == 1 & outbreak5@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)){
      tab$detect[j] <- 1
    }else {tab$detect[j] <- 0}
    tab$VP[j] <- sum(outbreak5@alarmall$ears.NB[[i]][j] == 1 & outbreak5@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    tab$VN[j] <- sum(outbreak5@alarmall$ears.NB[[i]][j] == 0 & outbreak5@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    
    tab$FP[j] <- sum(outbreak5@alarmall$ears.NB[[i]][j] == 1 & outbreak5@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    tab$FN[j] <- sum(outbreak5@alarmall$ears.NB[[i]][j] == 0 & outbreak5@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    
    tab$exces[j] <- sum(outbreak5@outbreaks[[i]][j])/sum(outbreak5@total[[i]][outbreak5@outbreaks[[i]][j]>0,j])
    
    tab$duree[j] <- sum(outbreak5@outbreaks[[i]][j]>0)
  }
  
  
  
  data7 <- rbind(data7, tab)
  
  rm(tab)
  
}


#### 7-MOVING-DAY COUNTS - outbreak 5 ####

##### SCENARIO 1 - outbreak 5 
# table for 100 simulation of one scenario
data8 <- data.frame(detect = rep(2, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    exces = rep(9999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("5", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("somme glissante", 100), 
                    
                    scenario = as.factor(rep(1, 100)),
                    theta = rep(6, 100), 
                    tendance = as.factor(rep(0, 100)),
                    k1 = as.factor(rep(1, 100)), 
                    k2 = as.factor(rep(2, 100)), 
                    phi = rep(2, 100)
                    
)
longueur_tab_alarme <- length(outbreak5@alarmRollSumall$ears.NB$scenario1$X1)
longueur_tab_outbk <- length(outbreak5@outbreaks$scenario1[,1])
for (i in 1:100){
  
  #detected/not detected:
  
  if (sum(outbreak5@alarmRollSumall$ears.NB$scenario1[i] == 1 & outbreak5@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)){
    data8$detect[i] <- 1
  }else {data8$detect[i] <- 0}
  
  data8$VP[i] <- sum(outbreak5@alarmRollSumall$ears.NB$scenario1[i] == 1 & outbreak5@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  data8$VN[i] <- sum(outbreak5@alarmRollSumall$ears.NB$scenario1[i] == 0 & outbreak5@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  
  data8$FP[i] <- sum(outbreak5@alarmRollSumall$ears.NB$scenario1[i] == 1 & outbreak5@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  data8$FN[i] <- sum(outbreak5@alarmRollSumall$ears.NB$scenario1[i] == 0 & outbreak5@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  
  
  
  
  data8$exces[i] <- sum(outbreak5@outbreaks$scenario1[i])/sum(outbreak5@total$scenario1[outbreak5@outbreaks$scenario1[i]>0,i])
  
  data8$duree[i] <- sum(outbreak5@outbreaks$scenario1[i]>0)
}

# all of the other scenarios



for (i in 2:16){
  
  tab <- data.frame(detect = rep(2, 100), 
                    exces = rep(9999, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("5", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("somme glissante", 100), 
                    
                    scenario = as.factor(rep(i, 100)),
                    theta = rep(theta.vecteur[i], 100), 
                    tendance = as.factor(rep(tendance.vecteur[i], 100)),
                    k1 = as.factor(rep(k1.vecteur[i], 100)), 
                    k2 = as.factor(rep(k2.vecteur[i], 100)), 
                    phi = rep(phi.vecteur[i], 100)
                    
  )
  
  longueur_tab_alarme <- length(outbreak5@alarmRollSumall$ears.NB[[i]]$X1)
  longueur_tab_outbk <- length(outbreak5@outbreaks[[i]][,1])
  for (j in 1:100){
    
    #detected/not detected:
    
    if (sum(outbreak5@alarmRollSumall$ears.NB[[i]][j] == 1 & outbreak5@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)){
      tab$detect[j] <- 1
    }else {tab$detect[j] <- 0}
    tab$VP[j] <- sum(outbreak5@alarmRollSumall$ears.NB[[i]][j] == 1 & outbreak5@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    tab$VN[j] <- sum(outbreak5@alarmRollSumall$ears.NB[[i]][j] == 0 & outbreak5@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    
    tab$FP[j] <- sum(outbreak5@alarmRollSumall$ears.NB[[i]][j] == 1 & outbreak5@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    tab$FN[j] <- sum(outbreak5@alarmRollSumall$ears.NB[[i]][j] == 0 & outbreak5@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    
    tab$exces[j] <- sum(outbreak5@outbreaks[[i]][j])/sum(outbreak5@total[[i]][outbreak5@outbreaks[[i]][j]>0,j])
    
    tab$duree[j] <- sum(outbreak5@outbreaks[[i]][j]>0)
  }
  
  
  
  data8 <- rbind(data8, tab)
  
  rm(tab)
  
}




######## WEEKLY COUNTS - outbreak 5 #######

##### SCENARIO 1 - outbreak 5 
# table for 100 simulation of one scenario
data9 <- data.frame(detect = rep(2, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    exces = rep(9999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("5", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("segmentation", 100), 
                    
                    scenario = as.factor(rep(1, 100)),
                    theta = rep(6, 100), 
                    tendance = as.factor(rep(0, 100)),
                    k1 = as.factor(rep(1, 100)), 
                    k2 = as.factor(rep(2, 100)), 
                    phi = rep(2, 100)
                    
)
longueur_tab_alarme <- length(outbreak5@alarmSegWeekall$ears.NB$scenario1$X1)
longueur_tab_outbk <- length(outbreak5@outbreaksWeek$scenario1[,1])
for (i in 1:100){
  
  #detected/not detected:
  
  if (sum(outbreak5@alarmSegWeekall$ears.NB$scenario1[i] == 1 & outbreak5@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)){
    data9$detect[i] <- 1
  }else {data9$detect[i] <- 0}
  
  data9$VP[i] <- sum(outbreak5@alarmSegWeekall$ears.NB$scenario1[i] == 1 & outbreak5@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  data9$VN[i] <- sum(outbreak5@alarmSegWeekall$ears.NB$scenario1[i] == 0 & outbreak5@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  
  data9$FP[i] <- sum(outbreak5@alarmSegWeekall$ears.NB$scenario1[i] == 1 & outbreak5@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  data9$FN[i] <- sum(outbreak5@alarmSegWeekall$ears.NB$scenario1[i] == 0 & outbreak5@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  
  
  
  
  data9$exces[i] <- sum(outbreak5@outbreaksWeek$scenario1[i])/sum(outbreak5@totalSegweek$scenario1[outbreak5@outbreaksWeek$scenario1[i]>0,i])
  
  data9$duree[i] <- sum(outbreak5@outbreaks$scenario1[i]>0)
}

# all of the other scenarios



for (i in 2:16){
  
  tab <- data.frame(detect = rep(2, 100), 
                    exces = rep(9999, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("5", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("segmentation", 100), 
                    
                    scenario = as.factor(rep(i, 100)),
                    theta = rep(theta.vecteur[i], 100), 
                    tendance = as.factor(rep(tendance.vecteur[i], 100)),
                    k1 = as.factor(rep(k1.vecteur[i], 100)), 
                    k2 = as.factor(rep(k2.vecteur[i], 100)), 
                    phi = rep(phi.vecteur[i], 100)
                    
  )
  
  longueur_tab_alarme <- length(outbreak5@alarmSegWeekall$ears.NB[[i]]$X1)
  longueur_tab_outbk <- length(outbreak5@outbreaksWeek[[i]][,1])
  for (j in 1:100){
    
    #detected/not detected:
    
    if (sum(outbreak5@alarmSegWeekall$ears.NB[[i]][j] == 1 & outbreak5@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)){
      tab$detect[j] <- 1
    }else {tab$detect[j] <- 0}
    tab$VP[j] <- sum(outbreak5@alarmSegWeekall$ears.NB[[i]][j] == 1 & outbreak5@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    tab$VN[j] <- sum(outbreak5@alarmSegWeekall$ears.NB[[i]][j] == 0 & outbreak5@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    
    tab$FP[j] <- sum(outbreak5@alarmSegWeekall$ears.NB[[i]][j] == 1 & outbreak5@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    tab$FN[j] <- sum(outbreak5@alarmSegWeekall$ears.NB[[i]][j] == 0 & outbreak5@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    
    tab$exces[j] <- sum(outbreak5@outbreaksWeek[[i]][j])/sum(outbreak5@totalSegweek[[i]][outbreak5@outbreaksWeek[[i]][j]>0,j])
    
    tab$duree[j] <- sum(outbreak5@outbreaks[[i]][j]>0)
  }
  
  
  
  data9 <- rbind(data9, tab)
  
  rm(tab)
  
}






#### DAILY COUNTS - outbreak 10 ####

##### SCENARIO 1 - outbreak 10 
# table for 100 simulation of one scenario
data10 <- data.frame(detect = rep(2, 100),
                     VP = rep(999, 100), 
                     VN = rep(999, 100),
                     FP = rep(999, 100), 
                     FN = rep(999, 100),
                     exces = rep(9999, 100),
                     duree= rep(999, 100), retard= rep(999, 100), 
                     taille.outbk = rep("10", 100),
                     
                     algorithme = rep("ears.NB", 100),
                     
                     time_seg = rep("daily", 100), 
                     
                     scenario = as.factor(rep(1, 100)),
                     theta = rep(6, 100), 
                     tendance = as.factor(rep(0, 100)),
                     k1 = as.factor(rep(1, 100)), 
                     k2 = as.factor(rep(2, 100)), 
                     phi = rep(2, 100)
                     
)
longueur_tab_alarme <- length(outbreak10@alarmall$ears.NB$scenario1$X1)
longueur_tab_outbk <- length(outbreak10@outbreaks$scenario1[,1])
for (i in 1:100){
  
  #detected/not detected:
  
  if (sum(outbreak10@alarmall$ears.NB$scenario1[i] == 1 & outbreak10@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)){
    data10$detect[i] <- 1
  }else {data10$detect[i] <- 0}
  
  data10$VP[i] <- sum(outbreak10@alarmall$ears.NB$scenario1[i] == 1 & outbreak10@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  data10$VN[i] <- sum(outbreak10@alarmall$ears.NB$scenario1[i] == 0 & outbreak10@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  
  data10$FP[i] <- sum(outbreak10@alarmall$ears.NB$scenario1[i] == 1 & outbreak10@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  data10$FN[i] <- sum(outbreak10@alarmall$ears.NB$scenario1[i] == 0 & outbreak10@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  
  
  
  
  data10$exces[i] <- sum(outbreak10@outbreaks$scenario1[i])/sum(outbreak10@total$scenario1[outbreak10@outbreaks$scenario1[i]>0,i])
  
  data10$duree[i] <- sum(outbreak10@outbreaks$scenario1[i]>0)
}

# all of the other scenarios



for (i in 2:16){
  
  tab <- data.frame(detect = rep(2, 100), 
                    exces = rep(9999, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("10", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("daily", 100), 
                    
                    scenario = as.factor(rep(i, 100)),
                    theta = rep(theta.vecteur[i], 100), 
                    tendance = as.factor(rep(tendance.vecteur[i], 100)),
                    k1 = as.factor(rep(k1.vecteur[i], 100)), 
                    k2 = as.factor(rep(k2.vecteur[i], 100)), 
                    phi = rep(phi.vecteur[i], 100)
                    
  )
  
  longueur_tab_alarme <- length(outbreak10@alarmall$ears.NB[[i]]$X1)
  longueur_tab_outbk <- length(outbreak10@outbreaks[[i]][,1])
  for (j in 1:100){
    
    #detected/not detected:
    
    if (sum(outbreak10@alarmall$ears.NB[[i]][j] == 1 & outbreak10@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)){
      tab$detect[j] <- 1
    }else {tab$detect[j] <- 0}
    tab$VP[j] <- sum(outbreak10@alarmall$ears.NB[[i]][j] == 1 & outbreak10@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    tab$VN[j] <- sum(outbreak10@alarmall$ears.NB[[i]][j] == 0 & outbreak10@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    
    tab$FP[j] <- sum(outbreak10@alarmall$ears.NB[[i]][j] == 1 & outbreak10@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    tab$FN[j] <- sum(outbreak10@alarmall$ears.NB[[i]][j] == 0 & outbreak10@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    
    tab$exces[j] <- sum(outbreak10@outbreaks[[i]][j])/sum(outbreak10@total[[i]][outbreak10@outbreaks[[i]][j]>0,j])
    
    tab$duree[j] <- sum(outbreak10@outbreaks[[i]][j]>0)
  }
  
  
  
  data10 <- rbind(data10, tab)
  
  rm(tab)
  
}


#### 7-MOVING-DAY COUNTS - outbreak 10 ####

##### SCENARIO 1 - outbreak 10 
# table for 100 simulation of one scenario
data11 <- data.frame(detect = rep(2, 100),
                     VP = rep(999, 100), 
                     VN = rep(999, 100),
                     FP = rep(999, 100), 
                     FN = rep(999, 100),
                     exces = rep(9999, 100),
                     duree= rep(999, 100), retard= rep(999, 100), 
                     taille.outbk = rep("10", 100),
                     
                     algorithme = rep("ears.NB", 100),
                     
                     time_seg = rep("somme glissante", 100), 
                     
                     scenario = as.factor(rep(1, 100)),
                     theta = rep(6, 100), 
                     tendance = as.factor(rep(0, 100)),
                     k1 = as.factor(rep(1, 100)), 
                     k2 = as.factor(rep(2, 100)), 
                     phi = rep(2, 100)
                     
)
longueur_tab_alarme <- length(outbreak10@alarmRollSumall$ears.NB$scenario1$X1)
longueur_tab_outbk <- length(outbreak10@outbreaks$scenario1[,1])
for (i in 1:100){
  
  #detected/not detected:
  
  if (sum(outbreak10@alarmRollSumall$ears.NB$scenario1[i] == 1 & outbreak10@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)){
    data11$detect[i] <- 1
  }else {data11$detect[i] <- 0}
  
  data11$VP[i] <- sum(outbreak10@alarmRollSumall$ears.NB$scenario1[i] == 1 & outbreak10@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  data11$VN[i] <- sum(outbreak10@alarmRollSumall$ears.NB$scenario1[i] == 0 & outbreak10@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  
  data11$FP[i] <- sum(outbreak10@alarmRollSumall$ears.NB$scenario1[i] == 1 & outbreak10@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  data11$FN[i] <- sum(outbreak10@alarmRollSumall$ears.NB$scenario1[i] == 0 & outbreak10@outbreaks$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  
  
  
  
  data11$exces[i] <- sum(outbreak10@outbreaks$scenario1[i])/sum(outbreak10@total$scenario1[outbreak10@outbreaks$scenario1[i]>0,i])
  
  data11$duree[i] <- sum(outbreak10@outbreaks$scenario1[i]>0)
}

# all of the other scenarios



for (i in 2:16){
  
  tab <- data.frame(detect = rep(2, 100), 
                    exces = rep(9999, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("10", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("somme glissante", 100), 
                    
                    scenario = as.factor(rep(i, 100)),
                    theta = rep(theta.vecteur[i], 100), 
                    tendance = as.factor(rep(tendance.vecteur[i], 100)),
                    k1 = as.factor(rep(k1.vecteur[i], 100)), 
                    k2 = as.factor(rep(k2.vecteur[i], 100)), 
                    phi = rep(phi.vecteur[i], 100)
                    
  )
  
  longueur_tab_alarme <- length(outbreak10@alarmRollSumall$ears.NB[[i]]$X1)
  longueur_tab_outbk <- length(outbreak10@outbreaks[[i]][,1])
  for (j in 1:100){
    
    #detected/not detected:
    
    if (sum(outbreak10@alarmRollSumall$ears.NB[[i]][j] == 1 & outbreak10@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)){
      tab$detect[j] <- 1
    }else {tab$detect[j] <- 0}
    tab$VP[j] <- sum(outbreak10@alarmRollSumall$ears.NB[[i]][j] == 1 & outbreak10@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    tab$VN[j] <- sum(outbreak10@alarmRollSumall$ears.NB[[i]][j] == 0 & outbreak10@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    
    tab$FP[j] <- sum(outbreak10@alarmRollSumall$ears.NB[[i]][j] == 1 & outbreak10@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    tab$FN[j] <- sum(outbreak10@alarmRollSumall$ears.NB[[i]][j] == 0 & outbreak10@outbreaks[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    
    tab$exces[j] <- sum(outbreak10@outbreaks[[i]][j])/sum(outbreak10@total[[i]][outbreak10@outbreaks[[i]][j]>0,j])
    
    tab$duree[j] <- sum(outbreak10@outbreaks[[i]][j]>0)
  }
  
  
  
  data11 <- rbind(data11, tab)
  
  rm(tab)
  
}




######## WEEKLY COUNTS - outbreak 10 #######

##### SCENARIO 1 - outbreak 10 
# table for 100 simulation of one scenario
data12 <- data.frame(detect = rep(2, 100),
                     VP = rep(999, 100), 
                     VN = rep(999, 100),
                     FP = rep(999, 100), 
                     FN = rep(999, 100),
                     exces = rep(9999, 100),
                     duree= rep(999, 100), retard= rep(999, 100), 
                     taille.outbk = rep("10", 100),
                     
                     algorithme = rep("ears.NB", 100),
                     
                     time_seg = rep("segmentation", 100), 
                     
                     scenario = as.factor(rep(1, 100)),
                     theta = rep(6, 100), 
                     tendance = as.factor(rep(0, 100)),
                     k1 = as.factor(rep(1, 100)), 
                     k2 = as.factor(rep(2, 100)), 
                     phi = rep(2, 100)
                     
)
longueur_tab_alarme <- length(outbreak10@alarmSegWeekall$ears.NB$scenario1$X1)
longueur_tab_outbk <- length(outbreak10@outbreaksWeek$scenario1[,1])
for (i in 1:100){
  
  #detected/not detected:
  
  if (sum(outbreak10@alarmSegWeekall$ears.NB$scenario1[i] == 1 & outbreak10@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)){
    data12$detect[i] <- 1
  }else {data12$detect[i] <- 0}
  
  data12$VP[i] <- sum(outbreak10@alarmSegWeekall$ears.NB$scenario1[i] == 1 & outbreak10@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  data12$VN[i] <- sum(outbreak10@alarmSegWeekall$ears.NB$scenario1[i] == 0 & outbreak10@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  
  data12$FP[i] <- sum(outbreak10@alarmSegWeekall$ears.NB$scenario1[i] == 1 & outbreak10@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] == 0)
  data12$FN[i] <- sum(outbreak10@alarmSegWeekall$ears.NB$scenario1[i] == 0 & outbreak10@outbreaksWeek$scenario1[(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,i] >0)
  
  
  
  
  data12$exces[i] <- sum(outbreak10@outbreaksWeek$scenario1[i])/sum(outbreak10@totalSegweek$scenario1[outbreak10@outbreaksWeek$scenario1[i]>0,i])
  
  data12$duree[i] <- sum(outbreak10@outbreaks$scenario1[i]>0)
}

# all of the other scenarios



for (i in 2:16){
  
  tab <- data.frame(detect = rep(2, 100), 
                    exces = rep(9999, 100),
                    VP = rep(999, 100), 
                    VN = rep(999, 100),
                    FP = rep(999, 100), 
                    FN = rep(999, 100),
                    duree= rep(999, 100), retard= rep(999, 100), 
                    taille.outbk = rep("10", 100),
                    
                    algorithme = rep("ears.NB", 100),
                    
                    time_seg = rep("segmentation", 100), 
                    
                    scenario = as.factor(rep(i, 100)),
                    theta = rep(theta.vecteur[i], 100), 
                    tendance = as.factor(rep(tendance.vecteur[i], 100)),
                    k1 = as.factor(rep(k1.vecteur[i], 100)), 
                    k2 = as.factor(rep(k2.vecteur[i], 100)), 
                    phi = rep(phi.vecteur[i], 100)
                    
  )
  
  longueur_tab_alarme <- length(outbreak10@alarmSegWeekall$ears.NB[[i]]$X1)
  longueur_tab_outbk <- length(outbreak10@outbreaksWeek[[i]][,1])
  for (j in 1:100){
    
    #detected/not detected:
    
    if (sum(outbreak10@alarmSegWeekall$ears.NB[[i]][j] == 1 & outbreak10@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)){
      tab$detect[j] <- 1
    }else {tab$detect[j] <- 0}
    tab$VP[j] <- sum(outbreak10@alarmSegWeekall$ears.NB[[i]][j] == 1 & outbreak10@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    tab$VN[j] <- sum(outbreak10@alarmSegWeekall$ears.NB[[i]][j] == 0 & outbreak10@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    
    tab$FP[j] <- sum(outbreak10@alarmSegWeekall$ears.NB[[i]][j] == 1 & outbreak10@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] == 0)
    tab$FN[j] <- sum(outbreak10@alarmSegWeekall$ears.NB[[i]][j] == 0 & outbreak10@outbreaksWeek[[i]][(longueur_tab_outbk-longueur_tab_alarme+1):longueur_tab_outbk,j] >0)
    
    tab$exces[j] <- sum(outbreak10@outbreaksWeek[[i]][j])/sum(outbreak10@totalSegweek[[i]][outbreak10@outbreaksWeek[[i]][j]>0,j])
    
    tab$duree[j] <- sum(outbreak10@outbreaks[[i]][j]>0)
  }
  
  
  
  data12 <- rbind(data12, tab)
  
  rm(tab)
  
}




dataNB <- rbind(dataNB, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11, data12)




dataNB$exces[is.na(dataNB$exces)] <- 0
dataNB$VP.FN <- dataNB$VP + dataNB$FN
dataNB$VN.FP <- dataNB$VN + dataNB$FP
dataNB$VP.FP <- dataNB$VP + dataNB$FP


data <- rbind(data, dataNB)

write.table(data, file="data_regression.csv", row.names = FALSE, sep="\t", dec=".", col.names = TRUE)
#write.table(data, file="data_regression.txt", row.names = FALSE, sep="\t", dec=".", col.names = TRUE)


# End of File