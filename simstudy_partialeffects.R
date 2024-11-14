########################################################--
########################################################--
#Partial Effects Simulation Study
#A.Keil, M.Kamenetsky
#Created: 2022-12-21
#Last Updated: 2023-04-28

########################################################--
########################################################--
rm(list=ls())
#Load Libraries ----
library(qgcomp)
library(future)
library(gWQS)
library(future.apply)
library(ggplot2)
library(dplyr)
library(here)
library(tidyr)
library(glmnet)
library(forcats)
library(caret)
library(MuMIn)
library(purrr);library(broom);library(matrixStats);library(plyr)
#plan(strategy=multicore)
`%ni%` <- Negate(`%in%`)
source("simstudy_partialeffectsfunctions.R")

args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])
set.seed(i)

resultsdir <- "../Results/r1/"
########################################################-
########################################################-
#SIMULATION STUDY ----
########################################################-
########################################################-

############################################################################################################################################################--
############################################################################################################################################################--
#BASE SIMULATIONS ----
############################################################################################################################################################--
############################################################################################################################################################--
#

#1) Base Case ----
ptm <- proc.time()
resl0 <- sim2(n=500, coef=c(0.1, -0.1, 0.1, -0.1, 0, 0, 0.0, 0.0, 0.0, 0.0),
                                      cor=c(0),
                                      V=10, future.seed = TRUE)

out.time=proc.time() - ptm
print(out.time)
res0 <- as.data.frame(resl0) %>%
 dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 0.2,
         trueneg=-0.2) %>%
  mutate(N=500,
         corrx1x2 = 0,
         corrx1x3 = 0,
         corrx1x4 = 0,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 4,
         Nexposures0 = 6,
         jointpos = 0.2,
         jointneg = -0.2,
         type="base0") #%>%
res <- cleandat(res0, truepos=0.2, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_clean_iter",i, ".csv"), row.names = FALSE)



#2) Correlated Exposures ----
resl0 = sim2(n=500,
                     coef=c(0.1, -0.1, 0.1, -0.1, 0, 0, 0.0, 0.0, 0.0, 0.0),
                     cor=c(0.75,0.5,0.8),
                     V=10, future.seed = TRUE)
res0 <- as.data.frame(resl0) %>%
  dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 0.2,
         trueneg=-0.2) %>%
  mutate(N=500,
         corrx1x2 = 0.75,
         corrx1x3 = 0.5,
         corrx1x4 = 0.8,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 4,
         Nexposures0 = 6,
         jointpos = 0.2,
         jointneg = -0.2,
         type="base0 correlated Exp")
res <- cleandat(res0, truepos=0.2, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_corr_clean_iter",i, ".csv"), row.names = FALSE)


#3) Big N ----
resl0 = sim2(n=1000,
                     coef=c(0.1, -0.1, 0.1, -0.1, 0, 0, 0.0, 0.0, 0.0, 0.0),
                     cor=c(0),
                     V=10, future.seed = TRUE)

res0 <- as.data.frame(resl0) %>%
  dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 0.2,
         trueneg=-0.2) %>%
  mutate(N=1000,
         corrx1x2 = 0,
         corrx1x3 = 0,
         corrx1x4 = 0,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 4,
         Nexposures0 = 6,
         jointpos = 0.2,
         jointneg = -0.2,
         type="base0 bigN")
res <- cleandat(res0, truepos=0.2, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_bigN_clean_iter",i, ".csv"), row.names = FALSE)

# #3aa) Bigger N----
resl0 = sim2(n=10000,
                     coef=c(0.1, -0.1, 0.1, -0.1, 0, 0, 0.0, 0.0, 0.0, 0.0),
                     cor=c(0),
                     V=10, future.seed = TRUE)


res0 <- as.data.frame(resl0) %>%
  dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 0.2,
         trueneg=-0.2) %>%
  mutate(N=10000,
         corrx1x2 = 0,
         corrx1x3 = 0,
         corrx1x4 = 0,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 4,
         Nexposures0 = 6,
         jointpos = 0.2,
         jointneg = -0.2,
         type="base0 biggerN")
res <- cleandat(res0, truepos=0.2, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_biggerN_clean_iter",i, ".csv"), row.names = FALSE)



#
#
#3a) Small N----
resl0 = sim2(n=300,
                     coef=c(0.1, -0.1, 0.1, -0.1, 0, 0, 0.0, 0.0, 0.0, 0.0),
                     cor=c(0),
                     V=10, future.seed = TRUE)
res0 <- as.data.frame(resl0) %>%
  dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 0.2,
         trueneg=-0.2) %>%
  mutate(N=300,
         corrx1x2 = 0,
         corrx1x3 = 0,
         corrx1x4 = 0,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 4,
         Nexposures0 = 6,
         jointpos = 0.2,
         jointneg = -0.2,
         type="base0 smallN")
res <- cleandat(res0, truepos=0.2, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_smallN_clean_iter",i, ".csv"), row.names = FALSE)


#
#4) Partial Effect Spread Out Among 4 Exposures----
resl0 = sim2(n=500,
                     coef=c(0.1, -0.05, 0.1, -0.05, -0.05, -0.05, 0.0, 0.0, 0.0, 0.0),
                     cor=c(0),
                     V=10, future.seed = TRUE)

res0 <- as.data.frame(resl0) %>%
  dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 0.2,
         trueneg=-0.2) %>%
  mutate(N=500,
         corrx1x2 = 0,
         corrx1x3 = 0,
         corrx1x4 = 0,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 6,
         Nexposures0 = 4,
         jointpos = 0.2,
         jointneg = -0.2,
         type="base0 spread4")
res <- cleandat(res0, truepos=0.2, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_spread4_clean_iter",i, ".csv"), row.names = FALSE)



#4)a) Partial Effect Spread Out Among 1 Exposures (1 Exposures) ----
resl0 = sim2(n=500,
                     coef=c(0.1, -0.2, 0.1, 0, 0.0, 0.0, 0.0, 0.0, 0.0),
                     cor=c(0),
                     V=10, future.seed = TRUE)
res0 <- as.data.frame(resl0) %>%
  dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 0.2,
         trueneg=-0.2) %>%
  mutate(N=500,
         corrx1x2 = 0,
         corrx1x3 = 0,
         corrx1x4 = 0,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 3,
         Nexposures0 = 7,
         jointpos = 0.2,
         jointneg = -0.2,
         type="base0 spread1")
res <- cleandat(res0, truepos=0.2, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_spread1_clean_iter",i, ".csv"), row.names = FALSE)

#5) Partial Effect Spread Out Among 8 Exposures ----

resl0 = sim2(n=500,
                      coef=c(0.1, -0.025, 0.1, -0.025, -0.025, -0.025,
                             -0.025, -0.025, -0.025, -0.025),
                      cor=c(0),
                      V=10, future.seed = TRUE)

res0 <- as.data.frame(resl0) %>%
  dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 0.2,
         trueneg=-0.2) %>%
  mutate(N=500,
         corrx1x2 = 0,
         corrx1x3 = 0,
         corrx1x4 = 0,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 10,
         Nexposures0 = 10,
         jointpos = 0.2,
         jointneg = -0.2,
         type="base0 spread8")
res <- cleandat(res0, truepos=0.2, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_spread8_clean_iter",i, ".csv"), row.names = FALSE)





#5) Partial Effect Spread Out Among More Exposures + Noise (10 Exposures)

resl0 = sim2(n=500,
                      coef=c(0.1, -0.025, 0.1, -0.025, -0.025, -0.025,
                             -0.025, -0.025, -0.025, -0.025, rep(0, 10)),
                      cor=c(0),
                      V=10, future.seed = TRUE)

res0 <- as.data.frame(resl0) %>%
  dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 0.2,
         trueneg=-0.2) %>%
  mutate(N=500,
         corrx1x2 = 0,
         corrx1x3 = 0,
         corrx1x4 = 0,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 10,
         Nexposures0 = 10,
         jointpos = 0.2,
         jointneg = -0.2,
         type="base0 spread+noise")
res <- cleandat(res0, truepos=0.2, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_spreadwithnoise_clean_iter",i, ".csv"), row.names = FALSE)




#6) Imbalance in Partial Effects ----
##a) Small Imbalance ----
## Make positive effect be 0.3 and negative effect 0.2 (positive effect is 0.1 > negative effect)
resl0 = sim2(n=500,
                     coef=c(0.15, -0.1, 0.15, -0.1, 0, 0, 0.0, 0.0, 0.0, 0.0),
                     cor=c(0),
                     V=10, future.seed = TRUE)#V isn't used in the simiter function, what is this?
res0 <- as.data.frame(resl0) %>%
  dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 0.3,
         trueneg=-0.2) %>%
  mutate(N=500,
         corrx1x2 = 0,
         corrx1x3 = 0,
         corrx1x4 = 0,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 4,
         Nexposures0 = 6,
         jointpos = 0.3,
         jointneg = -0.2,
         type="base0 small imbalance")
res <- cleandat(res0, truepos=0.3, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_imbalancepartialsmall_clean_iter",i, ".csv"), row.names = FALSE)


#
#
##b) Medium Imbalance ----
## Make positive effect be 0.6 and negative effect 0.2 (positive effect is 0.1 > negative effect)
resl0 = sim2(n=500,
                     coef=c(0.3, -0.1, 0.3, -0.1, 0, 0, 0.0, 0.0, 0.0, 0.0),
                     cor=c(0),
                     V=10, future.seed = TRUE)#V isn't used in the simiter function, what is this?
res0 <- as.data.frame(resl0) %>%
  dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 0.6,
         trueneg=-0.2) %>%
  mutate(N=500,
         corrx1x2 = 0,
         corrx1x3 = 0,
         corrx1x4 = 0,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 4,
         Nexposures0 = 6,
         jointpos = 0.6,
         jointneg = -0.2,
         type="base0 med imbalance")
res <- cleandat(res0, truepos=0.6, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_imbalancepartialmed_clean_iter",i, ".csv"), row.names = FALSE)


##c) Large Imbalance ----
## Make positive effect be 1 and negative effect 0.2 (positive effect is 0.1 > negative effect)
resl0 = sim2(n=500,
                     coef=c(0.5, -0.1, 0.5, -0.1, 0, 0, 0.0, 0.0, 0.0, 0.0),
                     cor=c(0),
                     V=10, future.seed = TRUE)#V isn't used in the simiter function, what is this?
res0 <- as.data.frame(resl0) %>%
  dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 1,
         trueneg=-0.2) %>%
  mutate(N=500,
         corrx1x2 = 0,
         corrx1x3 = 0,
         corrx1x4 = 0,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 4,
         Nexposures0 =6,
         jointpos = 1,
         jointneg = -0.2,
         type="base0 large imbalance") %>%
  filter(variable %ni% c("X", "trueoverall", "trueneg","truepos"))
res <- cleandat(res0, truepos=1, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_imbalancepartiallarge_clean_iter",i, ".csv"), row.names = FALSE)


#9) High, medium, low correlation among exposures ----
###a) Hi Correlation ----
resl0 = sim2(n=500,
                     coef=c(0.1, -0.1, 0.1, -0.1, 0, 0, 0.0, 0.0, 0.0, 0.0),
                     cor=c(.9, .9, .9),
                     V=10, future.seed = TRUE)#V isn't used in the simiter function, what is this?
res0 <- as.data.frame(resl0) %>%
  dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 0.2,
         trueneg=-0.2) %>%
  mutate(N=500,
         corrx1x2 = 0.9,
         corrx1x3 = 0.9,
         corrx1x4 = 0.9,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 4,
         Nexposures0 = 6,
         jointpos = 0.2,
         jointneg = -0.2,
         type="base0 hicorr Exp")
res <- cleandat(res0, truepos=0.2, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_hicorr_clean_iter",i, ".csv"), row.names = FALSE)


#
#
##b) Med Correlation ----
resl0 = sim2(n=500,
                     coef=c(0.1, -0.1, 0.1, -0.1, 0, 0, 0.0, 0.0, 0.0, 0.0),
                     #???there are 2 coefs here, one for each exposure, right?
                     cor=c(.6, .6, .6),#correlations between x1 and x2, x1 and x3,
                     ##and x1 and x4, right?
                     V=10, future.seed = TRUE)#V isn't used in the simiter function, what is this?

res0 <- as.data.frame(resl0) %>%
  dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 0.2,
         trueneg=-0.2) %>%
  mutate(N=500,
         corrx1x2 = 0.6,
         corrx1x3 = 0.6,
         corrx1x4 = 0.6,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 4,
         Nexposures0 = 6,
         jointpos = 0.2,
         jointneg = -0.2,
         type="base0 medcorr Exp")
res <- cleandat(res0, truepos=0.2, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_medcorr_clean_iter",i, ".csv"), row.names = FALSE)




##c) Small Correlation ----
resl0 = sim2(n=500,
                     coef=c(0.1, -0.1, 0.1, -0.1, 0, 0, 0.0, 0.0, 0.0, 0.0),
                     cor=c(.1, .1, .1),
                     V=10, future.seed = TRUE)
res0 <- as.data.frame(resl0) %>%
  dplyr::rename(value=resl0) %>%
  mutate(variable=attributes(resl0)$names,
         iter=i,
         truepos= 0.2,
         trueneg=-0.2) %>%
  mutate(N=500,
         corrx1x2 = 0.1,
         corrx1x3 = 0.1,
         corrx1x4 = 0.1,
         corrx1x5 = 0,
         corrx1x6 = 0,
         Nexposures = 4,
         Nexposures0 = 6,
         jointpos = 0.2,
         jointneg = -0.2,
         type="base0 smallcorr Exp")
res <- cleandat(res0, truepos=0.2, trueneg=-0.2)
write.csv(res, file=paste0(resultsdir, "partialsims_base0_smallcorr_clean_iter",i, ".csv"), row.names = FALSE)


