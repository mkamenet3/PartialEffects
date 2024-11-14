#############################################--
#############################################--
#############################################--
#NHANES Leukocytic Telomere Length data
#MKamenetsky
#Last modified: 2024-11-3
#############################################--
#############################################--
#############################################--
sink("LTLNHANES_analysis.txt")
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
library(MetBrewer)
library(RColorBrewer)
library(patchwork)
library(tableone)
#source("partialeffectsfuncs_20230428.R")
`%ni%` <- Negate(`%in%`)


mitro0 <- read.csv("../../Data/telomere_data/data/mitro_data_ccgibson_covars.csv") %>%
  dplyr::select(wtspo2yr, lbx074la:lbxf09la,
                lbxcot,lbxwbcsi, #added 06-26-23
                "ridageyr", "ridreth1", "dmdeduc2", "riagendr",
                "bmxbmi","telomean") %>%
  mutate(raceth = case_when( #recode race/ethnicity
    ridreth1 == 1 ~ "Mexican-American",
    ridreth1 == 2 ~ "Other Hispanic",
    ridreth1 == 3 ~ "Non-Hispanic White",
    ridreth1 == 4 ~ "Non-Hispanic Black",
    ridreth1 == 5 ~ "Other, Including Mixed Race"
  )) %>%
  mutate(educ = case_when( #recode education
    dmdeduc2 == 1 ~ "< 9th grade",
    dmdeduc2 == 2 ~ "9-11th grade",
    dmdeduc2 == 3 ~ "HS/GED",
    dmdeduc2 == 4 ~ "Some college/AA",
    dmdeduc2 == 5 ~ "College or above",
    dmdeduc2 == 7 ~ "Refused",
    dmdeduc2 == 9 ~ "Don't Know"
  )) %>%
  mutate(sex = case_when( #recode sex
    riagendr == 1 ~ "Male",
    riagendr == 2 ~ "Female"
  )) %>%
  mutate(age = ridageyr) %>%
  dplyr::select(-c(ridageyr, ridreth1, dmdeduc2, riagendr)) %>%
  dplyr::select(c(bmxbmi, telomean, raceth,
                  educ, sex, age, 
                  lbx074la, lbx099la,
                  lbx138la, lbx153la,
                  lbx170la, lbx180la,
                  lbx187la, lbx194la,
                  lbxpcbla,
                  lbxhxcla,
                  lbx118la,
                  lbxd03la,
                  lbxd05la,
                  lbxd07la,
                  lbxf03la,
                  lbxf04la,
                  lbxf05la,
                  lbxf08la,lbxpcbla,lbxhxcla,wtspo2yr,
                  lbxcot,lbxwbcsi)) %>%
  mutate(PCB74 = lbx074la,
         PCB99 = lbx099la,
         PCB138 = lbx138la,
         PCB153 = lbx153la,
         PCB170 = lbx170la,
         PCB180 = lbx180la,
         PCB187 = lbx187la,
         PCB194 = lbx194la,
         PCB126 = lbxpcbla,
         PCB118 = lbx118la,
         PCB169= lbxhxcla,
         `1,2,3,6,7,8-HxCDD` = lbxd03la,
         `1,2,3,4,6,7,8-HpCDD` = lbxd05la,
         `1,2,3,4,6,7,8,9-OCDD` = lbxd07la,
         `2,3,4,7,8-PeCDF` = lbxf03la,
         `1,2,3,4,7,8-HxCDF` = lbxf04la,
         `1,2,3,6,7,8-HxCDF` = lbxf05la,
         `1,2,3,4,6,7,8-HpCDF` =lbxf08la,
         cotinine = lbxcot,
         wbc = lbxwbcsi) %>%
  mutate(raceth= as.factor(raceth),
         educ=as.factor(educ),
         sex=as.factor(sex)) %>%
  mutate_at(vars(matches("lbx")),log) %>%
  rename_at(vars(matches("lbx")), .funs = funs(paste0("ln",.)))%>%
  mutate(age2 = age^2) %>%
  dplyr::select(-wtspo2yr) %>%
  mutate(
    ltelomean = log(telomean))

# Pre-quantize data and test/train split----
covars <- c("raceth", "educ","sex","age", "age2", "bmxbmi", "wbc")
expnms_neg <- c("lnlbx074la", "lnlbx099la",
                "lnlbx138la", "lnlbx153la",
                "lnlbx170la", "lnlbx180la",
                "lnlbx187la", "lnlbx194la")
expnms_pos <- c("lnlbxpcbla",
                "lnlbxhxcla",
                "lnlbx118la",
                "lnlbxd03la",
                "lnlbxd05la",
                "lnlbxd07la",
                "lnlbxf03la",
                "lnlbxf04la",
                "lnlbxf05la",
                "lnlbxf08la")
expnms2 = c(expnms_pos, expnms_neg)
nexp = length(expnms2) #number of exposures considered

mitro <- mitro0 %>% #select only the variables we want (the log-transfomred versions)
  dplyr::select(-c(cotinine, telomean, PCB74:`1,2,3,4,6,7,8-HpCDF`))
qdat <- quantize(mitro, expnms2)$data %>% # %>% log-transformed wbc (covariate)
  dplyr::select(c(-lnlbxwbcsi))
##########################################################################################--
##########################################################################################--
##########################################################################################--
#Make a Table1 ----
##########################################################################################--
##########################################################################################--
##########################################################################################--
tab0 <- mitro0 %>%
  dplyr::select(c(1:6, 27:46)) %>%
  dplyr::rename(Age = age,
                BMI = bmxbmi,
                Cotinine = cotinine,
                `Telomere Length` = telomean,
                `WBC` = wbc) %>%
  relocate(Age, BMI, Cotinine, `Telomere Length`, `WBC`) 

out_mitro2_tab <- print(CreateTableOne(data=tab0, includeNA = TRUE), missing=TRUE)
write.csv(out_mitro2_tab, file="../../Manuscript/EHP/Resubmission/tables/table1_ltl.csv")

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
  ltelomean ~ wqs+ bmxbmi + raceth + educ + sex + age + age2+ lnlbxcot +wbc, b1_pos = FALSE,
  na.action=na.exclude, b=100, mix_name=expnms2, q=NULL,
  data = qdat, family = "gaussian", seed = 2)))

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
ID_wqsss_neg <- wqsss_neg$final_weights$mix_name[which(wqsss_neg$final_weights$mean_weight>= 1/nexp)]


## Pos -----
try(suppressWarnings(
  wqsss_pos  <- gWQS::gwqs(
    ltelomean ~ wqs+ bmxbmi + raceth + educ + sex + age + age2+ lnlbxcot+wbc, b1_pos = TRUE,
    na.action=na.exclude, b=100, mix_name=expnms2, q=NULL,
    data =  qdat, family = "gaussian", seed = 3)
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
  wqsss_posCI <- confint(wqsss_pos)
}
ID_wqsss_pos <- wqsss_pos$final_weights$mix_name[which(wqsss_pos$final_weights$mean_weight>= 1/nexp)]

################################################################--
# WQSNS ----
################################################################--
## Neg -----
try(suppressWarnings(
  wqsns_neg <- gWQS::gwqs(
    ltelomean ~ wqs+ bmxbmi + raceth + educ + sex + age + age2+ lnlbxcot +wbc, b1_pos = FALSE,
    na.action=na.exclude, b=100, mix_name=expnms2, q=NULL, validation = 0.0,
    data = qdat, family = "gaussian", seed = 3)
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
  wqsns_pos <- gWQS::gwqs(
    ltelomean ~ wqs+ bmxbmi + raceth + educ + sex + age + age2+ lnlbxcot+wbc, b1_pos = TRUE,
    na.action=na.exclude, b=100, mix_name=expnms2, q=NULL, validation = 0.0,
    data = qdat, family = "gaussian", seed=3)
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
## Pos-----
qgcap_pos <- qgcomp.noboot(f=ltelomean~.,
                           expnms=expnms_pos, data=qdat, q=NULL)

psi_qgcap_pos  = qgcap_pos$psi[[1]]
qgcap_posLB = qgcap_pos$ci[1]
qgcap_posUB = qgcap_pos$ci[2]


## Neg-----
qgcap_neg <- qgcomp.noboot(f=ltelomean~.,
                           expnms=expnms_neg, data=qdat, q=NULL)

psi_qgcap_neg  = qgcap_neg$psi[[1]]
qgcap_negLB = qgcap_neg$ci[1]
qgcap_negUB = qgcap_neg$ci[2]

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
  set.seed(i+100) 
  trainidx <- sample(1:nrow(qdat), round(nrow(qdat)*0.4))
  valididx <- setdiff(1:nrow(qdat),trainidx)
  traindata <- qdat[trainidx,]
  validdata <- qdat[valididx,]
  nrow(traindata);nrow(validdata)
  
  qgcss <- qgcomp.partials(
    fun="qgcomp.noboot", f=ltelomean~., q=NULL,
    traindata=traindata[, c(expnms2,covars,"ltelomean")],
    validdata=validdata[, c(expnms2, covars,"ltelomean")],
    expnms=expnms2,
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
## Pos-----
try(suppressWarnings(
  wqsap_pos<- gWQS::gwqs(
    ltelomean ~ wqs+ bmxbmi + raceth + educ + sex + age + age2+ lnlbxcot+wbc + lnlbx074la +
      lnlbx099la + lnlbx138la + lnlbx153la + lnlbx170la +
      lnlbx180la + lnlbx187la + lnlbx194la, b1_pos = TRUE,
    na.action=na.exclude, b=100, mix_name=expnms_pos, q=NULL, validation = 0.0,
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

## Neg -----
try(suppressWarnings(
  wqsap_neg<- gWQS::gwqs(
    ltelomean ~ wqs+ bmxbmi + raceth + educ + sex + age + age2+ lnlbxcot+wbc + lnlbxpcbla +
      lnlbxhxcla +  lnlbx118la +
      lnlbxd03la + lnlbxd05la +
      lnlbxd07la + lnlbxf03la +
      lnlbxf04la + lnlbxf05la +
      lnlbxf08la , b1_pos = TRUE,
    na.action=na.exclude, b=100, mix_name=expnms_neg, q=NULL, validation = 0.0,
    data = qdat, family = "gaussian", seed=3)
))

psi_wqsap_neg =  ifelse(!exists("wqsap_neg"), NA, wqsap_neg$fit$coefficients[2])
if(!exists("wqsap_neg")){
  print("WQS A Priori:Neg doesn't exist")
  wqsap_negCI  <- matrix(c(NA,NA,NA,NA), ncol=2)
}else{
  summary(wqsap_neg )
  gwqs_barplot(wqsap_neg )
  gwqs_scatterplot(wqsap_neg )
  gwqs_fitted_vs_resid(wqsap_neg )
  wqsap_negCI <- confint(wqsap_neg )
}


################################################################--
#WQS2i ----
################################################################--
################################################################--


wqs2i<- gWQS::gwqs(
  ltelomean ~ pwqs+ nwqs+  bmxbmi + raceth + educ + sex + age + age2+ lnlbxcot +wbc, 
  na.action=na.exclude, mix_name=expnms2, q=NULL, 
  data = qdat, family = "gaussian", seed=3)
wqs2i_CI<- confint(wqs2i)

ID_wqs2i_pos <- wqs2i$final_weights$mix_name[wqs2i$final_weights$mean_weight_p>=1/18]
ID_wqs2i_neg <- wqs2i$final_weights$mix_name[wqs2i$final_weights$mean_weight_n>=1/18]



################################################################--
positives <- unique(c(ID_qgcss_pos, as.character(ID_wqs2i_pos),
                      as.character(ID_wqsns_pos), as.character(ID_wqsss_pos)))
negatives <- unique(c(ID_qgcss_neg, as.character(ID_wqs2i_neg),
                     as.character(ID_wqsss_neg)))
intersect(positives, negatives)


mycols <- c( "darkorchid4",  "plum4", "#FDBB84",
             "rosybrown1","#E34A33", "brown4")
neffect <- c(psi_qgcap_neg , qgcss$neg.fit$coef[2] , psi_wqsss_neg , psi_wqsns_neg, psi_wqsap_neg, wqs2i$fit$coefficients[3] )
peffect <- c(psi_qgcap_pos, qgcss$pos.fit$coef[2], psi_wqsss_pos , psi_wqsns_pos, psi_wqsap_pos ,wqs2i$fit$coefficients[2])
nLB <- c(qgcap_negLB, qgcss$neg.fit$ci[1], wqsss_negCI[2,1], wqsns_negCI[2,1], wqsap_negCI[2,1], wqs2i_CI[3,1] )
nUB <- c(qgcap_negUB, qgcss$neg.fit$ci[2] , wqsss_negCI[2,2],wqsns_negCI[2,2] , wqsap_negCI[2,2], wqs2i_CI[3,2])
pLB <- c(qgcap_posLB, qgcss$pos.fit$ci[1], wqsss_posCI[2,1], wqsns_posCI[2,1], wqsap_posCI[2,1], wqs2i_CI[2,1] )
pUB <- c(qgcap_posUB, qgcss$pos.fit$ci[2], wqsss_posCI[2,2], wqsns_posCI[2,2] , wqsap_posCI[2,2] , wqs2i_CI[2,2] )
method <- c("QGCAP", "QGCSS", "WQSSS", "WQSNS", "WQSAP", "WQS2i")

dat <- cbind.data.frame(method, neffect, nLB,nUB, peffect, pLB, pUB,mycols) %>%
  mutate(noeffect = ifelse(neffect==0,1,0)) %>%
  mutate(method=as.factor(method)) %>%
  mutate(method = fct_relevel(method, "WQSNS", after=3)) 


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
  ylim(-0.11, 0.11)+
  theme(legend.position = "bottom", legend.title = element_blank())

p2 <- ggplot()+
  geom_segment(data=dat,aes(y=pLB, yend=pUB, x=method, color=method, xend=method), size=5)+
  geom_point(data=dat,aes(y=peffect, x=method, shape=method), size=6) +
  theme_classic() +
  scale_x_discrete(limits=rev(levels(dat$method)))+
  coord_flip() +
  geom_hline(yintercept = 0, linetype="dashed") +
  xlab("")+ylab(expression(paste(psi^"+")))+
  scale_color_manual(values = mycols)+
  guides(color="none") +
  ggtitle("Positive Partial Effect")+
  ylim(-0.11, 0.11)+
  theme(legend.position = "bottom", legend.title = element_blank())

p1+p2+plot_annotation(title = "Telomere Length: Partial Effect Estimates and 95% CIs",
                      theme = theme(plot.title=element_text(size=16,hjust = 0.5)))
ggsave("../../Manuscript/EHP/Resubmission/figures/ltl_partialeffectsresults.pdf", height=4, width=8, units="in")
ggsave("../../Manuscript/EHP/Resubmission/figures/ltl_partialeffectsresults.eps", height=4, width=8, units="in")

ggsave("../../Manuscript/EHP/Resubmission/figures/Figure4.pdf", height=4, width=8, units="in")
ggsave("../../Manuscript/EHP/Resubmission/figures/Figure4.eps", height=4, width=8, units="in")

################################################################--
################################################################--
##Make Tables----
################################################################--
################################################################--
### Partial Effects data example table ----

datmem <- dat %>%
  mutate(membership_pos = case_when(
    method=="QGCAP" ~ paste(expnms_pos, collapse=", "),
    method=="QGCSS" ~ paste(ID_qgcss_pos, collapse=", "),
    method=="WQSSS" ~ paste0(ID_wqsss_pos, collapse = ", "),
    method=="WQSNS" ~ paste(ID_wqsns_pos, collapse = ", "),
    method=="WQSAP" ~ paste(expnms_pos, collapse=", "),
    method=="WQS2i" ~ paste(ID_wqs2i_pos, collapse=", "))) %>%
  mutate(membership_neg = case_when(
    method=="QGCAP" ~ paste(expnms_neg, collapse=", "),
    method=="QGCSS" ~ paste(ID_qgcss_neg, collapse=", "),
    method=="WQSSS" ~ paste(ID_wqsss_neg, collapse=", "),
    method=="WQSNS" ~ paste(ID_wqsns_neg, collapse=", "),
    method=="WQSAP" ~ paste(expnms_neg, collapse=", "),
    method=="WQS2i" ~ paste(ID_wqs2i_neg, collapse=", ")
  ))



datmem %>%
  mutate(psi_neg = paste0(round(neffect, 3), " (", round(nLB,3), ", ", round(nUB,3), ")")) %>%
  mutate(psi_pos = paste0(round(peffect, 3), " (", round(pLB,3), ", ", round(pUB,3), ")")) %>%
  dplyr::select(-c(noeffect,mycols,neffect, peffect, nLB, pLB, nUB, pUB)) %>%
  relocate(method,psi_neg, membership_neg, psi_pos, membership_pos) %>%
  gt() %>%
  tab_header(title="Telomere Length Data Analysis") %>%
  gt::gtsave(., "../../Manuscript/EHP/Resubmission/tables/ltl_datanalysis.rtf")


################################################################--
################################################################--
#Joint Effects plot ----
################################################################--
################################################################--
set.seed(2)

overall <- qgcomp.noboot(f=ltelomean~.,
                         expnms=expnms2, data=qdat, q=NULL)
qgcoverall <- coef(overall)[2]

load("../../Results/wqs_100seeds_telomeredata.RData")


psiest <- round(c(psi_qgcss , 
                  psi_qgcap_pos + psi_qgcap_neg, 
                  (wqspospsi+ wqsnegpsi), 
                  sum(psi_wqsns_neg, psi_wqsns_pos, na.rm = TRUE),
                  psi_wqsap_neg+ psi_wqsap_pos,
                  sum(wqs2i$fit$coefficients[2:3])),3)

method <- factor(c(rep("QGCSS",100),
                   "QGCAP",
                   rep("WQSSS",100), 
                   "WQSNS",
                   "WQSAP", 
                   "WQS2i"),
                 levels=unique(c(rep("QGCSS",100),"QGCAP","WQSSS", "WQSNS", "WQSAP", "WQS2i")),
                 ordered = TRUE)
overallpsibyqgc <- rep(qgcoverall,204)

psidat <- cbind.data.frame(method,psiest,overallpsibyqgc) %>%
  group_by(method) %>%
  mutate(meanpsiest = mean(psiest,na.rm=TRUE)) %>%
  mutate(id = case_when( #ID the corresponding single shot estimate
    round(psiest,3) == round(psi_qgcap_neg+psi_qgcap_pos,3) & method == "QGCAP" ~ 1,
    round(psiest,3) == round(psi_qgcss[1],3) & method== "QGCSS" ~ 1,
    round(psiest,3) == round((psi_wqsss_neg+ psi_wqsss_pos)[1],3) & method == "WQSSS" ~ 1,
    round(psiest,3) == round(sum(psi_wqsns_neg, psi_wqsns_pos, na.rm = TRUE),3) & method == "WQSNS" ~ 1,
    round(psiest,3) == round(psi_wqsap_pos+psi_wqsap_neg,3) & method == "WQSAP" ~ 1,
    round(psiest,3) == round( sum(wqs2i$fit$coefficients[2:3]),3) & method == "WQS2i" ~1
  )) %>%
  mutate(method = as.factor(method)) %>%
  mutate(method = fct_relevel(method, "WQSNS", after=3)) %>%
  mutate(method = fct_relevel(method, "QGCSS", after=1))

psidat %>%
  group_by(method) %>%
  summarize(q90 = round(quantile(psiest, probs=0.9, na.rm=TRUE),3),
            q10 = round(quantile(psiest, probs=0.1, na.rm = TRUE),3))

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
  ggtitle("Telomere Length Example:\nJoint Effect Estimate vs. Sum of Partial Effects")+
  theme(legend.position = "bottom", legend.title = element_blank())


ggsave(filename = "../../Manuscript/EHP/Resubmission/figures/sumpartialsvsjoint_telomere.pdf", height = 5, width = 8, units = "in", dpi=300)
ggsave(filename = "../../Manuscript/EHP/Resubmission/figures/sumpartialsvsjoint_telomere.eps", height = 5, width = 8, units = "in", dpi=300)

ggsave(filename = "../../Manuscript/EHP/Resubmission/figures/FigureS4.pdf", height = 5, width = 8, units = "in", dpi=300)
ggsave(filename = "../../Manuscript/EHP/Resubmission/figures/FigureS4.eps", height = 5, width = 8, units = "in", dpi=300)


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
  ggtitle("Telomere Length Example:\nJoint Effect Estimate vs. Sum of Partial Effects")+
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(filename = "../../Manuscript/EHP/Resubmission/figures/sumpartialsvsjoint_telomere_BW.pdf", height = 5, width = 8, units = "in", dpi=300)
ggsave(filename = "../../Manuscript/EHP/Resubmission/figures/sumpartialsvsjoint_telomere_BW.eps", height = 5, width = 8, units = "in", dpi=300)

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
  #mutate(asd = ifelse(is.na))
  mutate(out = paste0(round(amean,3), " (", ifelse(is.na(asd),"",round(asd,3)), ")")) %>%
  dplyr::select(method, out) %>%
  gt() %>%
  #tab_spanner(label="(Log) MSE", columns=3:4, id=1) %>%
  #tab_spanner(label="Bias", columns=5:6, id=2) %>%
  tab_header(title="Telomere Length Eestimates") %>%
  cols_label(method="Method",
             out = "Psi Est (SD (Psi Est)") %>%
  gt::gtsave(., "../../Manuscript/EHP/Resubmission/tables/telomerelength_results.rtf")


sink()