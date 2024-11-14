#############################################--
#############################################--
#############################################--
#Analysis of isoprostane-8 (LIFECODES Study)
#MKamenetsky
#Last modified: 2024-11-1
#############################################--
#############################################--
#############################################--
sink("LIFECODES_analysis.txt")

#Load packages ----
library(tidyverse)
library(ggcorrplot)
library(kableExtra)
library(qgcomp)
library(gWQS)
library(MetBrewer)
# table one
library(tableone)
library(gt)
library(glmnet)
library(data.table)
library(RColorBrewer)
`%ni%` <- Negate(`%in%`)
library(forcats)
library(patchwork)
library(stringr)

#Prep data ----
dat1  <- read.csv("../../Data/keil_lifec_20230208.csv") %>%
  mutate(iso = X_15F2T_ISOP_T1) %>%
  filter(iso !="N") %>%
  mutate(iso = as.numeric(iso),
         concept_yr = as.factor(concept_yr),
         mat_edu = as.factor(mat_edu),
         alc_cat = as.factor(alc_cat),
         smpreg = as.factor(smpreg)) %>%
  mutate(sumDEHP = ((MEHP_T1/278) + (MEHHP_T1/294) + (MEOHP_T1/292) + (MECPP_T1/308))*308) 
dat2 <- dat1 %>%
  dplyr::select(c(Subject_ID:TCS_T1,iso, sumDEHP)) %>%
  drop_na(MCNP_T1:sumDEHP, age:concept_yr) #drop NA based on analytes AND covariates
sg_median = median(dat2$sg_T1, na.rm=TRUE)
dat <- dat2 %>%
  pivot_longer(cols = c(MCNP_T1:TCS_T1,sumDEHP), names_to="exposure", values_to="expvalue") %>%
  mutate(lnexpvalue = log(expvalue)) %>% #natural log transform and specific gravity adjustment
  mutate(sgexpvalue = expvalue * ((sg_median-1)/(sg_T1 - 1))) %>%
  mutate(sglnexpvalue = log(sgexpvalue),
         sgiso = iso * ((sg_median-1)/(sg_T1-1)),
         lnsgiso = log(sgiso))


# Pre-quantize data and test/train split----
dat2 <- dat %>%
  filter(exposure %ni% c("MEHP_T1", "MEHHP_T1", "MEOHP_T1", "MECPP_T1"))%>% #omit the analytes that go into sumDEHP
  pivot_wider(id_cols = Subject_ID,
              names_from=exposure,
              values_from=c(sglnexpvalue, sgiso)) %>%
  mutate(sgiso = sgiso_MCNP_T1) %>% #take the first one because of repetiations
  dplyr::select(-c(sgiso_MCNP_T1:sgiso_sumDEHP)) %>%
  mutate_if(is.factor, na_if, y='')
#merge in covariates back
dat2 <- merge(dat2, dat1[, c("Subject_ID","age", "BMI_prepreg", "mat_edu", "alc_cat",
                             "smpreg", "ga_t1", "concept_yr")], by="Subject_ID") %>%
  mutate_if(is.factor, funs(factor(replace(., .=="", NA)))) %>% #convert to factor
  drop_na()

#set exposure names to specific-gravity corrected log tansformed analytes
expnames = grep("sglnexpvalue", names(dat2), value = TRUE)
covars <- c("age", "BMI_prepreg", "mat_edu", "alc_cat",
            "smpreg", "ga_t1", "concept_yr")
nexp = length(expnames) #number of exposures considered

#quantized dataset
qdat <- quantize(dat2, expnames)$data %>%
  dplyr::select(-c(Subject_ID))

##########################################################################################--
##########################################################################################--
##########################################################################################--
#Make a Table1 ----
##########################################################################################--
##########################################################################################--
##########################################################################################--
tab0 <- dat %>%
  filter(exposure %ni% c("MEHP_T1", "MEHHP_T1", "MEOHP_T1", "MECPP_T1"))%>%
  pivot_wider(id_cols = Subject_ID,
              names_from=exposure,
              values_from=c(expvalue, iso)) %>%
  mutate(iso = iso_MCNP_T1) %>% #take the first one
  dplyr::select(-c(iso_MCNP_T1:iso_sumDEHP))
#merge in covariates back
dat2 <- merge(tab0, dat1[, c("Subject_ID","age", "BMI_prepreg", "mat_edu", "alc_cat",
                             "smpreg", "ga_t1", "concept_yr")], by="Subject_ID") %>%
  rename_all(~stringr::str_replace(.,"^expvalue_","")) %>%
  mutate_if(is.factor, funs(factor(replace(., .=="", NA)))) %>%
  dplyr::select(-Subject_ID) %>%
  relocate(age, BMI_prepreg, ga_t1, concept_yr,smpreg, mat_edu, alc_cat,
           iso)
out_lifecodes_tab <- print(CreateTableOne(data=dat2, includeNA = TRUE), missing=TRUE)
write.csv(out_lifecodes_tab, file="../../Manuscript/EHP/Resubmission/tables/table1_lifecodes.csv")

##########################################################################################--
##########################################################################################--
##########################################################################################--
#Run ANALYSES ----
##########################################################################################--
##########################################################################################--
##########################################################################################--

################################################################--
# WQSSS ----
################################################################--
## Neg -----
try(suppressWarnings(wqsss_neg<- gWQS::gwqs(
  sgiso ~ wqs + age + BMI_prepreg+ mat_edu + alc_cat + smpreg + ga_t1 + concept_yr, b1_pos = FALSE,
  na.action=na.exclude, b=100, mix_name=expnames, q=NULL,
  data = qdat[, c(expnames,covars,"sgiso")], family = "gaussian", seed=2)))


psi_wqsss_neg =  ifelse(!exists("wqsss_neg"), NA, wqsss_neg$fit$coefficients[2])
if(!exists("wqsss_neg")){
  print("WQS:Neg doesn't exist")
  wqsss_negCI <- matrix(c(NA,NA,NA,NA), ncol=2)
  wqsnegpsi <- 0
}else{
  summary(wqsss_neg)
  gwqs_barplot(wqsss_neg)
  gwqs_scatterplot(wqsss_neg)
  gwqs_fitted_vs_resid(wqsss_neg)
  wqsss_negCI <- confint(wqsss_neg)
}
ID_wqsss_neg <- ifelse(exists("wqsss_neg"),
                       wqsss_neg$final_weights$mix_name[which(wqsss_neg$final_weights$mean_weight>= 1/nexp)],
                       NA)


## Pos -----
try(suppressWarnings(
  wqsss_pos <- gWQS::gwqs(
    sgiso ~ wqs + age + BMI_prepreg+ mat_edu + alc_cat + smpreg + ga_t1 + concept_yr, b1_pos = TRUE,
    na.action=na.exclude, b=100, mix_name=expnames, q=NULL,
    data =  qdat[, c(expnames,covars,"sgiso")], family = "gaussian", seed=2)
))
psi_wqsss_pos =  ifelse(!exists("wqsss_pos"), NA, wqsss_pos$fit$coefficients[2])
if(!exists("wqsss_pos")){
  print("WQS:Pos doesn't exist")
  wqsss_posCI <- matrix(c(NA,NA,NA,NA), ncol=2)
}else{
  summary(wqsss_pos)
  gwqs_barplot(wqsss_pos)
  gwqs_scatterplot(wqsss_pos)
  gwqs_fitted_vs_resid(wqsss_pos)
  
  #gwqs_weights_tab(splitres5)
  wqsss_posCI <- confint(wqsss_pos)
}
ID_wqsss_pos <- wqsss_pos$final_weights$mix_name[which(wqsss_pos$final_weights$mean_weight>= 1/nexp)]


################################################################--
# WQSNS ----
################################################################--
## Neg -----
try(suppressWarnings(
  wqsns_neg <- gWQS::gwqs(
    sgiso ~ wqs+ age+ BMI_prepreg + mat_edu + alc_cat + smpreg + ga_t1 + concept_yr, b1_pos = FALSE,
    na.action=na.exclude, b=100, mix_name=expnames, q=NULL, validation = 0.0,
    data = qdat[, c(expnames,covars,"sgiso")], family = "gaussian", seed=2)
))
psi_wqsns_neg =  ifelse(!exists("wqsns_neg"), NA, wqsns_neg$fit$coefficients[2])


if(!exists("wqsns_neg")){
  print("WQSNS:Neg doesn't exist")
  wqsns_negCI <- matrix(c(NA,NA,NA,NA), ncol=2)
  wqsnegNSpsi <-0
}else{
  summary(wqsns_neg)
  gwqs_barplot(wqsns_neg)
  gwqs_scatterplot(wqsns_neg)
  gwqs_fitted_vs_resid(wqsns_neg)
  wqsns_negCI <- confint(wqsns_neg)
}
ID_wqsns_neg <- ifelse(exists("wqsns_neg"),
                       wqsns_neg$final_weights$mix_name[which(wqsns_neg$final_weights$mean_weight>= 1/nexp)],
                       NA)

## Pos -----
try(suppressWarnings(
  wqsns_pos<- gWQS::gwqs(
    sgiso ~ wqs+ age+ BMI_prepreg + mat_edu + alc_cat + smpreg + ga_t1 + concept_yr, b1_pos = TRUE,
    na.action=na.exclude, b=100, mix_name=expnames, q=NULL, validation = 0.0,
    data = qdat[, c(expnames,covars,"sgiso")], family = "gaussian", seed=2)
))
psi_wqsns_pos=  ifelse(!exists("wqsns_pos"), NA, wqsns_pos$fit$coefficients[2]) 
if(!exists("wqsns_pos")){
  print("WQSNS:Pos doesn't exist")
  wqsns_posCI <- matrix(c(NA,NA,NA,NA), ncol=2)
  wqsposNSpsi <-0
}else{
  summary(wqsns_pos)
  gwqs_barplot(wqsns_pos)
  gwqs_scatterplot(wqsns_pos)
  gwqs_fitted_vs_resid(wqsns_pos)
  wqsns_posCI <- confint(wqsns_pos)
}
ID_wqsns_pos <- wqsns_pos$final_weights$mix_name[which(wqsns_pos$final_weights$mean_weight>= 1/nexp)]

################################################################--
#QGCAP ----
################################################################--
expnms_pos <- names(qdat)[1:21]
set.seed(2)
qgcap_pos <- qgcomp.noboot(f=sgiso~ ., q=NULL, 
                             expnms=expnms_pos, data=qdat)

psi_qgcap_pos  = qgcap_pos$psi[[1]]
qgcap_posLB = qgcap_pos$ci[1]
qgcap_posUB = qgcap_pos$ci[2]

#we assume a priori that no metabolites are in the negative partial effect
psi_qgcap_neg= NA
qgcap_negLB = NA
qgcap_negUB = NA



################################################################--
#QGCSS ----
################################################################--
nseeds <- 100

#set empties
qgcss_neg=   rep(NA, length(nseeds))
qgcss_negLB = rep(NA, length(nseeds))
qgcss_negUB = rep(NA, length(nseeds))
qgcss_pos =   rep(NA, length(nseeds))
qgcss_posLB = rep(NA, length(nseeds))
qgcss_posUB = rep(NA, length(nseeds))
psi_qgcss <- rep(NA, length(nseeds))

for(i in 1:nseeds){
  set.seed(i)
  trainidx <- sample(1:nrow(qdat), round(nrow(qdat)*0.4))
  valididx <- setdiff(1:nrow(qdat),trainidx)
  traindata <- qdat[trainidx,]
  validdata <- qdat[valididx,]
  nrow(traindata);nrow(validdata)
  
  qgcss <- qgcomp.partials(
    fun="qgcomp.noboot", f=sgiso~., q=NULL,
    traindata=traindata[, c(expnames,covars,"sgiso")],
    validdata=validdata[, c(expnames, covars,"sgiso")],
    expnms=expnames,
    .fixbreaks=TRUE
  )
  qgcss_neg[i] =   ifelse(is.null(qgcss$neg.fit$coef[2]), 0, qgcss$neg.fit$coef[2])
  qgcss_negLB[i] = ifelse(is.null(qgcss$neg.fit$coef[2]),NA, qgcss$neg.fit$ci[1])
  qgcss_negUB[i] = ifelse(is.null(qgcss$neg.fit$coef[2]),NA,qgcss$neg.fit$ci[2])
  qgcss_pos[i] =   ifelse(is.null(qgcss$pos.fit$coef[2]), 0, qgcss$pos.fit$coef[2])
  qgcss_posLB[i] = ifelse(is.null(qgcss$pos.fit$coef[2]),NA, qgcss$pos.fit$ci[1])
  qgcss_posUB[i] = ifelse(is.null(qgcss$pos.fit$coef[2]),NA,qgcss$pos.fit$ci[2])
  psi_qgcss[i] <- qgcss_neg[i] + qgcss_pos[i] #overall effect from summation
}
ID_qgcss_neg <-  attributes(qgcss$neg.fit$neg.weights)$names
ID_qgcss_pos <-  attributes(qgcss$pos.fit$pos.weights)$names


################################################################--
#WQSAP ----
################################################################--

try(suppressWarnings(
  wqsap_pos<- gWQS::gwqs(
    sgiso ~ wqs+ age+ BMI_prepreg + mat_edu + alc_cat + smpreg + ga_t1 + concept_yr, b1_pos = TRUE,
    na.action=na.exclude, b=100, mix_name=expnames, q=NULL, validation = 0.0,
    data = qdat, family = "gaussian", seed=3)
))
psi_wqsap_pos =  ifelse(!exists("wqsap_pos"), NA, wqsap_pos$fit$coefficients[2])
if(!exists("wqsap_pos")){
  print("WQS A Priori:Pos doesn't exist")
  wqsap_posCI  <- matrix(c(NA,NA,NA,NA), ncol=2)
}else{
  summary(wqsap_pos )
  gwqs_barplot(wqsap_pos )
  gwqs_scatterplot(wqsap_pos )
  gwqs_fitted_vs_resid(wqsap_pos )
  wqsap_posCI <- confint(wqsap_pos )
}
#No WQSAP-Neg because assuming all mixture are positive effect
psi_wqsap_neg <- NA 
psi_wqsap_negCI <- c(NA, NA)


################################################################--
#WQS2i ----
################################################################--
wqs2i<- gWQS::gwqs(
  sgiso ~ pwqs+ nwqs+   age+ BMI_prepreg + mat_edu + alc_cat + smpreg + ga_t1 + concept_yr, 
  na.action=na.exclude, mix_name=expnames, q=NULL, 
  data = qdat, family = "gaussian", seed=4,rh=1, future.seed = TRUE,b=1)
wqs2i_CI<- confint(wqs2i)

ID_wqs2i_pos <- wqs2i$final_weights$mix_name[wqs2i$final_weights$mean_weight_p>=1/nexp]
ID_wqs2i_neg <- wqs2i$final_weights$mix_name[wqs2i$final_weights$mean_weight_n>=1/nexp]

################################################################--
#what metabolites are identified uniquely?
positives <- unique(c(ID_qgcss_pos, as.character(ID_wqs2i_pos),
                      as.character(ID_wqsns_pos), as.character(ID_wqsss_pos)))
negatives <- unique(c(ID_qgcss_neg, as.character(ID_wqs2i_neg),
                     as.character(ID_wqsss_neg)))
intersect(positives, negatives)


##########################################################################################--
##########################################################################################--
#Partial Effects Plot ----
##########################################################################################--
##########################################################################################--
mycols <- c( "darkorchid4",  "plum4", "#FDBB84",
             "rosybrown1","#E34A33", "brown4")
neffect <- c(psi_qgcap_neg , qgcss$neg.fit$coef[2] , psi_wqsss_neg , psi_wqsns_neg, psi_wqsap_neg, wqs2i$fit$coefficients[3] )
peffect <- c(psi_qgcap_pos, qgcss$pos.fit$coef[2], psi_wqsss_pos , psi_wqsns_pos, psi_wqsap_pos ,wqs2i$fit$coefficients[2])
nLB <- c(qgcap_negLB, qgcss$neg.fit$ci[1], wqsss_negCI[2,1], wqsns_negCI[2,1], NA, wqs2i_CI[3,1] )
nUB <- c(qgcap_negUB, qgcss$neg.fit$ci[2] , wqsss_negCI[2,2],wqsns_negCI[2,2] , NA, wqs2i_CI[3,2])
pLB <- c(qgcap_posLB, qgcss$pos.fit$ci[1], wqsss_posCI[2,1], wqsns_posCI[2,1], wqsap_posCI[2,1], wqs2i_CI[2,1] )
pUB <- c(qgcap_posUB, qgcss$pos.fit$ci[2], wqsss_posCI[2,2], wqsns_posCI[2,2] , wqsap_posCI[2,2] , wqs2i_CI[2,2] )
method <- c("QGCAP", "QGCSS", "WQSSS", "WQSNS", "WQSAP", "WQS2i")

dat <- cbind.data.frame(method, neffect, nLB,nUB, peffect, pLB, pUB,mycols) %>%
  mutate(noeffect = ifelse(neffect==0,1,0)) %>%
  mutate(method=as.factor(method)) %>%
  mutate(method = fct_relevel(method, "WQSNS", after=3)) 
dat

p1 <-ggplot()+
  geom_segment(data=dat,aes(y=nLB, yend=nUB, x=method, color=method, xend=method), size=5)+
  geom_point(data=dat,aes(y=neffect, x=method, shape=method), size=6) +
  scale_x_discrete(limits = rev(levels(dat$method)))+
  theme_classic() +
  coord_flip() +
  geom_hline(yintercept = 0, linetype="dashed") +
  xlab("")+ylab(expression(paste(psi^"-")))+
  scale_color_manual(values = mycols)+
  guides(color="none") +
  ggtitle("Negative Partial Effect")+
  ylim(-0.3, 0.6)+
  theme(legend.position = "bottom", legend.title = element_blank())



p2 <-ggplot(data=dat)+
  geom_segment(aes(y=pLB, yend=pUB, x=method, color=method), xend=method, size=5)+
  geom_point(aes(y=peffect, x=method, shape=method), size=6) +
  theme_classic() +
  scale_x_discrete(limits=rev(levels(dat$method)))+
  coord_flip() +
  geom_hline(yintercept = 0, linetype="dashed") +
  xlab("")+ylab(expression(paste(psi^"+")))+
  scale_color_manual(values = mycols)+
  guides(color="none") +
  ggtitle("Positive Partial Effect")+
  ylim(-0.3, 0.6)+
  theme(legend.position = "bottom", legend.title = element_blank())


p1+p2+plot_annotation(title = "Isoprostane: Partial Effects and 95% CIs",
                      theme = theme(plot.title=element_text(size=16,hjust = 0.5)))
ggsave("../../Manuscript/EHP/Resubmission/figures/lifecodes_partialeffectsresults.pdf", height=4, width=8, units="in", dpi=300)
ggsave("../../Manuscript/EHP/Resubmission/figures/lifecodes_partialeffectsresults.eps", height=4, width=8, units="in", dpi=300)
ggsave("../../Manuscript/EHP/Resubmission/figures/Figure3.pdf", height=4, width=8, units="in", dpi=300)
ggsave("../../Manuscript/EHP/Resubmission/figures/Figure3.eps", height=4, width=8, units="in", dpi=300)

################################################################--
################################################################--
##Make Tables----
################################################################--
################################################################--
### Partial Effects data example table ----

datmem <- dat %>%
  mutate(membership_pos = case_when(
    method=="QGCAP" ~ paste(expnames, collapse=", "),
    method=="QGCSS" ~ paste(ID_qgcss_pos, collapse=", "),
    method=="WQSSS" ~ paste0(ID_wqsss_pos, collapse = ", "),
    method=="WQSNS" ~ paste(ID_wqsns_pos, collapse = ", "),
    method=="WQSAP" ~ paste(expnames, collapse=", "),
    method=="WQS2i" ~ paste(ID_wqs2i_pos, collapse=", "))) %>%
  mutate(membership_neg = case_when(
    method=="QGCAP" ~ "",
    method=="QGCSS" ~ paste(ID_qgcss_neg, collapse=", "),
    method=="WQSSS" ~ paste(ID_wqsss_neg, collapse=", "),
    method=="WQSNS" ~ paste(ID_wqsns_neg, collapse=", "),
    method=="WQSAP" ~ "",
    method=="WQS2i" ~ paste(ID_wqs2i_neg, collapse=", ")
  ))
datmem

datmem %>%
  mutate(psi_neg = paste0(round(neffect, 3), " (", round(nLB,3), ", ", round(nUB,3), ")")) %>%
  mutate(psi_pos = paste0(round(peffect, 3), " (", round(pLB,3), ", ", round(pUB,3), ")")) %>%
  dplyr::select(-c(noeffect,mycols,neffect, peffect, nLB, pLB, nUB, pUB)) %>%
  relocate(method,psi_neg, membership_neg, psi_pos, membership_pos) %>%
  gt() %>%
  tab_header(title="LifeCodes Data Analysis") %>%
  gt::gtsave(., "../../Manuscript/EHP/Resubmission/tables/lifecodes_datanalysis.rtf")



################################################################--
################################################################--
#Joint Effects plot ----
################################################################--
################################################################--

set.seed(2)
overall <- qgcomp.noboot(f=sgiso~., q=NULL, expnms=expnames, data=qdat[, c(expnames,covars,"sgiso")])
qgcoverall <- coef(overall)[2]

load("../../Results/wqs_100seeds_isoprostanedata.RData") #these results were run on biowulf server


psiest <- round(c(psi_qgcss , 
                  psi_qgcap_pos, 
                  (wqspospsi+ wqsnegpsi), 
                  sum(psi_wqsns_neg, psi_wqsns_pos, na.rm = TRUE),
                  0+ psi_wqsap_pos,
                  sum(wqs2i$fit$coefficients[2:3])),3)


method <- factor(c(rep("QGCSS",100),
                   "QGCAP",
                   rep("WQSSS",100), 
                   "WQSNS",
                   "WQSAP", 
                   "WQS2i"),
                 levels=unique(c(rep("QGCSS",100),"QGCAP",rep("WQSSS", 100), "WQSNS", "WQSAP", "WQS2i")),
                 ordered = TRUE)
overallpsibyqgc <- rep(qgcoverall,204)

psidat <- cbind.data.frame(method,psiest,overallpsibyqgc) %>%
  group_by(method) %>%
  dplyr::mutate(meanpsiest = mean(psiest,na.rm=TRUE)) %>%
  dplyr::mutate(id = case_when( #ID the corresponding single shot estimate
    round(psiest,3) == round(psi_qgcap_pos,3) & method == "QGCAP" ~ 1,
    round(psiest,3) == round(psi_qgcss[100],3) & method== "QGCSS" ~ 1,
    round(psiest,3) == round((wqspospsi+ wqsnegpsi)[100],3) & method == "WQSSS" ~ 1,
    round(psiest,3) == round(sum(psi_wqsns_neg, psi_wqsns_pos, na.rm = TRUE),3) & method == "WQSNS" ~ 1,
    round(psiest,3) == round(psi_wqsap_pos,3) & method == "WQSAP" ~ 1,
    round(psiest,3) == round( sum(wqs2i$fit$coefficients[2:3]),3) & method == "WQS2i" ~1
  )) %>%
  dplyr::mutate(method = as.factor(method)) %>%
  dplyr::mutate(method = fct_relevel(method, "WQSNS", after=3)) %>%
  dplyr::mutate(method = fct_relevel(method, "QGCSS", after=1)) 

psidat %>%
  group_by(method) %>%
  dplyr::summarize(q90 = round(quantile(psiest, probs=0.9, na.rm=TRUE),3),
            q10 = round(quantile(psiest, probs=0.1, na.rm = TRUE),3), 
            avg = round(quantile(psiest, probs=0.5, na.rm = TRUE),3))
with(psidat, tapply(psiest, method, summary))
psidat



ggplot(data=psidat)+
  geom_point(aes(x=method, y=psiest, fill=method, color=method, shape=method), size=5, aes=0.5,stat="identity") +
  geom_point(aes(x=method, y=meanpsiest, fill=method, shape=method),color="black", size=15, pch=21) +
  geom_point(data=subset(psidat, id==1), 
             aes(x=method, y=psiest, fill=method, shape=method),color="black", size=5,  stroke=2) +
  
  geom_hline(aes(yintercept = overallpsibyqgc), linetype="dashed", size=1.5) +
  scale_x_discrete(limits=rev(levels(psidat$method)))+
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = mycols)+
  scale_fill_manual(values = mycols)+
  guides(color="none")+xlab("") +ylab(expression(paste(psi))) +
  ggtitle("Isoprostane Example:\nJoint Effect Estimate vs. Sum of Partial Effects")+
  theme(legend.position = "bottom",legend.title=element_blank())




ggsave(filename = "../../Manuscript/EHP/Resubmission/sumpartialsvsjoint_isoprostane.pdf", height = 5, width = 8, units = "in", dpi=300)
ggsave(filename = "../../Manuscript/figures/sumpartialsvsjoint_isoprostane.eps", height = 5, width = 8, units = "in", dpi=300)

ggsave(filename = "../../Manuscript/EHP/Resubmission/FigureS2.pdf", height = 5, width = 8, units = "in", dpi=300)
ggsave(filename = "../../Manuscript/EHP/Resubmission/FigureS2.eps", height = 5, width = 8, units = "in", dpi=300)

################################################################--
##BW Version ----
################################################################--
mycolsbw <- c( "black",  "grey40", "grey50","grey60","grey70", "grey80")

ggplot(data=psidat)+
  geom_point(aes(x=method, y=psiest, fill=method, color=method, shape=method), size=5, aes=0.5,stat="identity") +
  geom_point(aes(x=method, y=meanpsiest, fill=method, shape=method),color="black", size=15, pch=21) +
  geom_point(data=subset(psidat, id==1), 
             aes(x=method, y=psiest, fill=method, shape=method),color="black", size=5,  stroke=2) +
  
  geom_hline(aes(yintercept = overallpsibyqgc), linetype="dashed", size=1.5) +
  scale_x_discrete(limits=rev(levels(psidat$method)))+
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = mycolsbw)+
  scale_fill_manual(values = mycolsbw)+
  guides(color="none")+xlab("") +ylab(expression(paste(psi))) +
  ggtitle("Isoprostane Example:\nJoint Effect Estimate vs. Sum of Partial Effects")+
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(filename = "../../Manuscript/EHP/Resubmission/figures/sumpartialsvsjoint_isoprostane_BW.pdf", height = 5, width = 8, units = "in", dpi=300)
ggsave(filename = "../../Manuscript/EHP/Resubmission/figures/sumpartialsvsjoint_isoprostane_BW.eps", height = 5, width = 8, units = "in", dpi=300)

ggsave(filename = "../../Manuscript/EHP/Resubmission/FigureS2_BW.pdf", height = 5, width = 8, units = "in", dpi=300)
ggsave(filename = "../../Manuscript/EHP/Resubmission/FigureS2_BW.eps", height = 5, width = 8, units = "in", dpi=300)

################################################################--
################################################################--
##Make Tables----
################################################################--
################################################################--
### Partial Effects data example table ----
###Joint vs partial effects tables ----
psidat %>%
  group_by(method) %>%
  dplyr::select(method, psiest)%>%
  dplyr::summarise(amean = mean(psiest, na.rm=TRUE),
                   asd = sd(psiest, na.rm=TRUE)) %>%
  mutate(out = paste0(round(amean,3), " (", ifelse(is.na(asd),"",round(asd,3)), ")")) %>%
  dplyr::select(method, out) %>%
  gt() %>%
  #tab_spanner(label="(Log) MSE", columns=3:4, id=1) %>%
  #tab_spanner(label="Bias", columns=5:6, id=2) %>%
  tab_header(title="LifeCodes") %>%
  cols_label(method="Method",
             out = "Psi Est (SD)") %>%
  gt::gtsave(., "../../Manuscript/EHP/Resubmission/tables/lifecodes_results.rtf")

sink()

