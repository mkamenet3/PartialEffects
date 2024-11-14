#####################################################################--
#####################################################################--
#Functions used for partial effects simulation study
#####################################################################--
#####################################################################--
#2023-04-28: this version sets elastic net alpha to 0.5 so no need to
#do cv to determine optimal alpha, just cv.glmnet for lambda


sim2<- function(n=1000, coef=c(0.25,-0.25,0,0), cor=c(.75, 0), V=10, future.seed=TRUE){
  #create quantized data
  dat = qgcomp::simdata_quantized(n=n, outcomtype="continuous", cor=cor, 
                                  b0=0, coef=coef, q=4)
  #extract all and +/- exposures
  expnms = grep("x[0-9]+", names(dat), value = TRUE)
  expnms_pos <- which(coef>0)
  expnms_neg <- which(coef<0)
  
  #QGC Overall
  overall <- qgcomp.noboot(f=y~., q=NULL, expnms=expnms, data=dat)
  print("Overall joint effect estimated")
  
  #QGC AP
  qgcap_pos <- qgcomp.noboot(f=y~., q=NULL, expnms=expnms[expnms_pos], data=dat)
  qgcap_neg <- qgcomp.noboot(f=y~., q=NULL, expnms=expnms[expnms_neg], data=dat)
  print("QGCAP estimated")
  
  
  
  # partial effects using 40% training/60% validation split
  splitdat = split_data(dat, prop.train = .4)
  
  #estimate via model-averaging on train/test split
  fullmod <- lm(y~., data=splitdat$traindata,na.action = "na.fail")
  dd <- dredge(fullmod, rank="AIC")
  ma_out <- model.avg(dd)
  ma_coefs <- ma_out$coefficients[1,]# extract full model
  ma_exneg <- grep("x[0-9]+",names(which(ma_coefs<0)), value=TRUE)
  ma_expos <- grep("x[0-9]+",names(which(ma_coefs>0)), value=TRUE)
  #run qgc on coefs determined by MA
  qgc_ma_pos <-qgcomp.noboot(f=y~., q=NULL, expnms=expnms[which(expnms%in%ma_expos)], data=splitdat$validdata)
  qgc_ma_neg <-qgcomp.noboot(f=y~., q=NULL, expnms=expnms[which(expnms%in%ma_exneg)], data=splitdat$validdata)
  print("QGCMA estimated")
  
  #QGC-Enet
  qgc_enet <- qgcomp.partials2(
    fun="qgcomp.noboot", f=y~., q=NULL, 
    selectfun = "qgcomp_enet2",
    traindata=splitdat$traindata,
    validdata=splitdat$validdata, 
    expnms=expnms,
    .fixbreaks=TRUE
  ) 
  print("QGC-enet estimated")

  #QGC_SS
  qgc_ss <- qgcomp.partials(
    fun="qgcomp.noboot", f=y~., q=NULL, 
    traindata=splitdat$traindata,
    validdata=splitdat$validdata, 
    expnms=expnms,
    .fixbreaks=TRUE
  ) 
  print("QGCSS estimated")
  
  #WQS: Neg
  try(suppressWarnings(
    wqs_neg <- gWQS::gwqs(
      y ~ wqs, b1_pos = FALSE,
      na.action=na.exclude, b=100, mix_name=expnms, q=NULL,
      data = dat, family = "gaussian", future.seed=TRUE)
  ))
  print("WQS-Neg estimated")
  #WQS: Pos
  try(suppressWarnings(
    wqs_pos <- gWQS::gwqs(
      y ~ wqs, b1_pos = TRUE,
      na.action=na.exclude, b=100, mix_name=expnms, q=NULL,
      data = dat, family = "gaussian", future.seed=TRUE)
  ))
  print("WQS-Pos estimated")
  #WQS NS: Neg
  try(suppressWarnings(
    wqsns_neg <- gWQS::gwqs(
      y ~ wqs, b1_pos = FALSE,
      na.action=na.exclude, b=100, mix_name=expnms, q=NULL, validation = 0.0,
      data = dat, family = "gaussian", future.seed=TRUE)
  ))
  print("WQSNS-Neg estimated")
  #WQS NS: Pos
  try(suppressWarnings(
    wqsns_pos <- gWQS::gwqs(
      y ~ wqs, b1_pos = TRUE,
      na.action=na.exclude, b=100, mix_name=expnms, q=NULL, validation = 0.0,
      data = dat, family = "gaussian", future.seed=TRUE)
  ))
  print("WQSNS-Pos estimated")

  ###Add in new WQS methods
  #RH-WQS
  try(suppressWarnings(wqsrh_neg <- gWQS::gwqs(y ~ wqs, 
                     mix_name = expnms, 
                     data = dat, 
                     q = NULL, 
                     validation = 0.6, # 40% trainging data, 60% testing data
                     b = 100, # 100 bootstraps recommended, but gWQS runs faster with fewer
                     rh = 100,
                     b1_pos = FALSE, future.seed=TRUE)))
  print("RH-WQS-Neg estimated")
  
  try(suppressWarnings(wqsrh_pos <- gWQS::gwqs(y ~ wqs, 
                            mix_name = expnms, 
                            data = dat, 
                            q = NULL, 
                            validation = 0.6, # 40% trainging data, 60% testing data
                            b = 100, # 100 bootstraps recommended, but gWQS runs faster with fewer.
                            rh=100,
                            b1_pos = TRUE,  future.seed=TRUE)))
  print("RH-WQS-Pos estimated")
  #WQSAP
  #browser()
  try(suppressWarnings(
    wqsap_pos <- gWQS::gwqs(
    y ~ wqs + x2 +x4, b1_pos = TRUE,
    na.action=na.exclude, b=100, mix_name=expnms[expnms_pos], q=NULL, validation = 0.0,
    data = dat, family = "gaussian",future.seed=TRUE)
    ))
  print("WQSAP-Pos estimated")
  try(suppressWarnings(
    wqsap_neg <- gWQS::gwqs(
    y ~ wqs + x1 +x3, b1_pos = FALSE,
    na.action=na.exclude, b=100, mix_name=expnms[expnms_neg], q=NULL, validation = 0.0,
    data = dat, family = "gaussian", future.seed=TRUE)
    ))
  print("WQSAP-Neg estimated")
  
  #WQS2i double index model
  #browser()
  try(suppressWarnings(
    wqs_2i<- gWQS::gwqs(
      y ~ pwqs + nwqs, b1_pos = TRUE,
      na.action=na.exclude, b=1, mix_name=expnms, q=NULL,
      data = dat, family = "gaussian", future.seed=TRUE)
  ))
  print("WQS-double index estimated")
  
  
  
  
  
  
  c(
    trueoverall = sum(coef),
    res_overall = overall$psi[1], 
    res_qgcap_pos = qgcap_pos$psi[1],
    res_qgcap_neg = qgcap_neg$psi[1],
    trueneg = sum(coef[coef<0]),
    
    res_qgcma_pos = qgc_ma_pos$psi[1],
    res_qgcma_neg = qgc_ma_neg$psi[1],
    

    res_qgcenet_neg = ifelse(is.null(qgc_enet$neg.fit$coef[2]), 0, qgc_enet$neg.fit$coef[2]),
    res_qgcss_neg =   ifelse(is.null(qgc_ss$neg.fit$coef[2]), 0, qgc_ss$neg.fit$coef[2]), 
    res_wqs_neg =  ifelse(!exists("wqs_neg"), 0, wqs_neg$fit$coefficients[2]), 
    res_wqsns_neg =  ifelse(!exists("wqsns_neg"), 0, wqsns_neg$fit$coefficients[2]), 
    truepos= sum(coef[coef>0]),
    res_qgcenet_pos = ifelse(is.null(qgc_enet$pos.fit$coef[2]), 0, qgc_enet$pos.fit$coef[2]),
    res_qgcss_pos =   ifelse(is.null(qgc_ss$pos.fit$coef[2]), 0, qgc_ss$pos.fit$coef[2]), 
    res_wqs_pos =  ifelse(!exists("wqs_pos"), 0, wqs_pos$fit$coefficients[2]),
    res_wqsns_pos =  ifelse(!exists("wqsns_pos"), 0, wqsns_pos$fit$coefficients[2]) ,
    #additional methods
    res_wqsrh_neg = ifelse(!exists("wqsrh_neg"),0, wqsrh_neg$fit$coefficients[2,1]),
    res_wqsrh_pos =  ifelse(!exists("wqsrh_pos"),0, wqsrh_pos$fit$coefficients[2,1]),
    res_wqsap_pos = ifelse(!exists("wqsap_pos"),0, wqsap_pos$fit$coefficients[2]),
    res_wqsap_neg = ifelse(!exists("wqsap_neg"),0, wqsap_neg$fit$coefficients[2]),
    res_wqs2i_pos = wqs_2i$fit$coefficients[2],
    res_wqs2i_neg = wqs_2i$fit$coefficients[3]
    
  )
}





cleandat <- function(res0, truepos, trueneg){
 res <- res0 %>%
    mutate(
      name = case_when(
        variable %in% c("res_qgcap_pos.psi1", "res_qgcap_neg.psi1") ~ "QGC: A Priori",
        variable %in% c("res_qgcenet_neg", "res_qgcenet_pos") ~ "QGC: enet",
        variable %in% c("res_qgcss_neg", "res_qgcss_pos") ~ "QGCSS",
        variable %in% c("res_wqs_neg", "res_wqs_pos") ~ "WQSSS",
        variable %in% c("res_wqsns_neg", "res_wqsns_pos") ~ "WQSNS",
        variable %in% c("res_wqsrh_neg", "res_wqsrh_pos") ~ "RH-WQS",
        variable %in% c("res_wqsap_pos", "res_wqsap_neg") ~ "WQSAP",
        variable %in% c("res_qgcma_pos.psi1", "res_qgcma_neg.psi1") ~ "QGCMA",
        variable == "res_overall.psi1" ~ "Overall (Joint Effect)",
        variable %in% c("res_wqs2i_pos.pwqs", "res_wqs2i_neg.nwqs") ~ "WQS2i"
      )
    ) %>%
    filter(variable %ni% c("X", "trueoverall", "trueneg","truepos")) %>%
    mutate(posneg = ifelse(variable %in% c("res_qgcap_neg.psi1",
                                           "res_qgcenet_neg",
                                           "res_qgcss_neg",
                                           "res_qgcma_neg.psi1",
                                           "res_wqs_neg",
                                           "res_wqsns_neg",
                                           "res_wqsrh_neg",
                                           "res_wqsap_neg",
                                           "res_qgcma_neg",
                                           "res_wqs2i_neg.nwqs"), "Negative Effect", "Positive Effect")) %>%
    mutate(posneg = ifelse(name=="Overall (Joint Effect)", "Joint Effect", posneg)) %>%
    mutate(bias = case_when(
      posneg == "Positive Effect" ~ value - truepos,
      posneg == "Negative Effect" ~ value - trueneg,
      posneg == "Joint Effect" ~ value - 0
    ))
 return(res)
}








qgcomp_enet2 = function(f, expnms, data, q, breaks=NULL, ...){
  
  newform <- terms(f, data = data)
  if(!(is.null(q))){
    ql <- qgcomp::quantize(data, expnms, q, breaks)
    qdata <- ql$data
    br <- ql$breaks
  } else{
    qdata=data
    br = NULL
  }
  
  X = model.matrix(newform, data = qdata)
  y = qdata[,all.vars(f)[1],drop=TRUE]
  fit <- cv.glmnet(X[,-1],y, alpha=0.5, intercept=TRUE, ...) # drop constant term from model matrix
  whichind = which(fit$lambda == fit$lambda.min)
  coefs = coef(fit$glmnet.fit)[,whichind]
  epsilon = .Machine$double.eps
  pc = coefs[(coefs > epsilon) & (names(coefs) %in% expnms)]
  nc = coefs[(coefs < -epsilon) & (names(coefs) %in% expnms)]
  nullc = coefs[(coefs > -epsilon) & (coefs < epsilon) & (names(coefs) %in% expnms)]
  list(
    coefs = coefs,
    pos.weights  = pc/sum(pc),
    neg.weights  = nc/sum(nc)
  )
}



qgcomp.partials2 = function (fun = c("qgcomp.noboot", "qgcomp.cox.noboot", "qgcomp.zi.noboot"), 
                             selectfun = c("qgcomp_enet2"),
                             traindata = NULL, validdata = NULL, expnms = NULL, .fixbreaks = TRUE, 
                             .globalbreaks = FALSE, ...) 
{
  if (is.null(traindata) | is.null(validdata)) 
    stop("traindata and validdata must both be specified")
  traincall <- validcall <- match.call(expand.dots = TRUE)
  #  start modification
  #droppers <- match(c("traindata", "validdata", ".fixbreaks", 
  #                    ".globalbreaks", "fun"), names(traincall), 0L)
  droppers <- match(c("traindata", "validdata", ".fixbreaks", 
                      ".globalbreaks", "fun", "selectfun"), names(traincall), 0L)
  #  end modification
  traincall[["data"]] <- eval(traincall[["traindata"]], parent.frame())
  validcall[["data"]] <- eval(validcall[["validdata"]], parent.frame())
  #browser()
  
  traincall <- traincall[-c(droppers)]
  validcall <- validcall[-c(droppers)]
  hasbreaks = ifelse("breaks" %in% names(traincall), TRUE, 
                     FALSE)
  hasq = ifelse("q" %in% names(traincall), TRUE, FALSE)
  qnull = ifelse(is.null(traincall$q), TRUE, FALSE)
  if (hasbreaks && .fixbreaks) 
    .fixbreaks = FALSE
  if (hasq) {
    if (!hasbreaks && qnull && .fixbreaks) 
      .fixbreaks = FALSE
    if (!hasbreaks && qnull && .globalbreaks) 
      .globalbreaks = FALSE
  }
  #  start modification
  #if (is.function(fun)) {
  #  traincall[[1L]] <- validcall[[1L]] <- fun
  #}
  #else {
  #  traincall[[1L]] <- validcall[[1L]] <- as.name(fun[1])
  #}
  if (is.function(fun)) {
    validcall[[1L]] <- fun
  } else
  {
    validcall[[1L]] <- as.name(fun[1])
  }
  if (is.function(selectfun)) {
    traincall[[1L]] <- selectfun
  }else
  {
    traincall[[1L]]  <- as.name(selectfun[1])
  }
  # end modification
  if (.globalbreaks) {
    globalcall <- traincall
    globalcall[["data"]] <- eval(rbind(traincall[["data"]], 
                                       validcall[["data"]]), parent.frame())
    global.fit = eval(globalcall, parent.frame())
    if (!is.null(global.fit$breaks)) {
      validcall[["breaks"]] = global.fit$breaks
      validcall[["q"]] = NULL
      traincall[["breaks"]] = global.fit$breaks
      traincall[["q"]] = NULL
    }
  }
  train.fit = eval(traincall, parent.frame())
  if (.fixbreaks && !.globalbreaks) {
    validcall$breaks = train.fit$breaks
    validcall$q = NULL
  }
  posnms = expnms[is.element(expnms, names(train.fit$pos.weights))]
  negnms = expnms[is.element(expnms, names(train.fit$neg.weights))]
  if (length(posnms) == 1 && all(posnms == c("count", "zero"))) {
    posnms = names(train.fit$pos.weights$count)
    negnms = names(train.fit$neg.weights$count)
  }
  res = list(train.fit = train.fit)
  res$negmix <- res$posmix <- "none"
  if (length(posnms) > 0) {
    res$posmix = posnms
    poscall <- validcall
    if (!is.null(poscall$breaks) && (.fixbreaks || .globalbreaks)) {
      posidx <- which(expnms %in% posnms)
      poscall$breaks <- poscall$breaks[posidx]
    }
    vc = as.list(poscall)
    vc$expnms = c(posnms)
    res$pos.fit <- eval(as.call(vc), parent.frame())
  }
  if (length(negnms) > 0) {
    res$negmix = negnms
    negcall <- validcall
    if (!is.null(negcall$breaks) && (.fixbreaks || .globalbreaks)) {
      negidx <- which(expnms %in% negnms)
      negcall$breaks <- negcall$breaks[negidx]
    }
    vc = as.list(negcall)
    vc$expnms = c(negnms)
    res$neg.fit <- eval(as.call(vc), parent.frame())
  }
  class(res) <- "qgcompmultifit"
  res
}


