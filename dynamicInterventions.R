## -----------------------------------------------------------------------------
## intervention info without age limitation
## ----
library(data.table)
prob_selection <- 0.8
datX <- data.table(x = 1:5)
datX[,`:=`(x1=min(1.2,x), x2=pmin(1.2,x))]
datX[x==1,`:=`(y1=x)]
datX[x!=1,`:=`(x=y1)]
datX

# library(vcd)
# datPP0 <- datPP[t0==0]
# datPP0[
#     ,`:=`(`Diabetes mellitus`=diabetesRecdStatusLabel,
#           `LDL-C group`=ifelse(LDLC<1.8,"<1.8",
#                          ifelse(LDLC>=1.8 & LDLC<2.6,"<2.6",
#                                 ifelse(LDLC>=2.6 & LDLC<3.4,"<3.4","3.4+"))),
#           `Risk stratification`=riskStatusLabel)]
# int.1a <- xtabs(~`Risk stratification`+`LDL-C group`+`Diabetes mellitus`,
#                 data=datPP0) 
# int.1a <- structable(
#     ~`Risk stratification`+`LDL-C group`+`Diabetes mellitus`,
#     data=datPP0)
# pushViewport(viewport(layout = grid.layout(ncol=2)))
# mosaic(int.1a,split_vertical=TRUE,cex=0.8)
# mNames <- getNames()
# high.DM <- mNames[grepl("=High",mNames)] 
# high.DM <- high.DM[grepl("=Yes",high.DM)]
# high.DM <- high.DM[!grepl("<1.8",high.DM)]
# # high.DM <- gsub("cell")
# grid.edit(high.DM[1],gp=gpar(fill="cyan"))
# grid.edit(high.DM[2],gp=gpar(fill="cyan"))
# grid.edit(high.DM[3],gp=gpar(fill="cyan"))
# high <- c(mNames[grepl("=High",mNames)],
#           mNames[grepl("=Intermediate",mNames)])

# survfit(Surv(fAnyEvent_y,fAnyEvent)~riskStatusLabel,data=datOutcome)
ldlc_int <- function(newdf, pool, intvar, intvals, time_name, t) {
    newdf[, `:=`(LDLC = pmax(get(intvar)-1,0),
                 lipidDrugStatus = 1)]
}
ldlc_int.high <- function(newdf, pool, intvar, intvals, time_name, t) {
    riskStatusVal <- intvals[[1]]

    # reduce 1 mmol/L for those at high risk
    newdf[riskStatus==riskStatusVal[3], 
          `:=`(LDLC = pmax(get(intvar)-1,0),
               lipidDrugStatus = 1)]
}
## 2020 CSC guideline-based intervention
# No.1
csc2020_int.1a <- function(newdf, pool, intvar, intvals, time_name, t) {
    # p1014: Lipid-lowering target
    # intvals = list(riskStatus=c(1,2,3),diabeteRecdStatus=1)
    riskStatusVal <- intvals[[1]]  
    diabeteStatus <- intvals[[2]]
    ageVal <- intvals[[3]]
    
    # Eligible age and ldlc
    # set.seed(1)
    newdf[,`:=`(cond_met=((age1992+t0)<ageVal),
                ldlcLow=runif(.N,0.1,3.4),
                ldlcInt=runif(.N,0.1,2.6),
                ldlcHigh=runif(.N,0.1,1.8))]
    
    # low-risk: LDL-C < 3.4 (IIa,B)
    newdf[riskStatus==riskStatusVal[1] #& cond_met==TRUE
          ,`:=`(LDLC = pmin(ldlcLow, get(intvar)),
                lipidDrugStatus = 1)]
    # moderate-risk: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[2] #& cond_met==TRUE
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                lipidDrugStatus = 1)]
    # high-risk w/o diabetes: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus!=diabeteStatus & cond_met==TRUE
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                lipidDrugStatus = 1)]
    # high-risk w diabetes: LDLC < 1.8 or 50%+ reduction
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus==diabeteStatus #& cond_met==TRUE
          ,`:=`(LDLC = pmax(ldlcHigh, get(intvar)*0.5),
                lipidDrugStatus = 1)]
}
csc2020_int.1b <- function(newdf, pool, intvar, intvals, time_name, t) {
    # p1014: Lipid-lowering target
    # intvals = list(riskStatus=c(1,2,3),diabeteRecdStatus=1)
    riskStatusVal <- intvals[[1]]  
    diabeteStatus <- intvals[[2]]  
    ageVal <- intvals[[3]]
    
    # 80% chance to adhere the dynamic intervention 
    # set.seed(1)
    newdf[,`:=`(cond_sel=rbinom(.N, 1, 0.8),
                cond_met=((age1992+t0)<ageVal),
                ldlcLow=runif(.N,0.1,3.4),
                ldlcInt=runif(.N,0.1,2.6),
                ldlcHigh=runif(.N,0.1,1.8),
                lipidDrugStatus = 1)]
    
    # low-risk: LDL-C < 3.4 (IIa,B)
    newdf[riskStatus==riskStatusVal[1] & cond_sel==1 #& cond_met==TRUE 
          ,`:=`(LDLC = pmin(ldlcLow,get(intvar)),
                lipidDrugStatus = 1)]
    # moderate-risk: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[2] & cond_met==1 #& cond_met==TRUE
          ,`:=`(LDLC = pmin(ldlcInt,get(intvar)),
                lipidDrugStatus = 1)]
    # high-risk w/o diabetes: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus!=diabeteStatus &
              cond_sel==1 #& cond_met==TRUE 
          ,`:=`(LDLC = pmin(ldlcInt,get(intvar)),
                lipidDrugStatus = 1)]
    # high-risk w diabetes: LDLC < 1.8 or 50%+ reduction
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus==diabeteStatus &
              cond_sel==1 #& cond_met==TRUE
          ,`:=`(LDLC = pmax(ldlcHigh,get(intvar)*0.5),
                lipidDrugStatus = 1)]
}

csc2020_int.1b20 <- function(newdf, pool, intvar, intvals, time_name, t) {
    # p1014: Lipid-lowering target
    # intvals = list(riskStatus=c(1,2,3),diabeteRecdStatus=1)
    riskStatusVal <- intvals[[1]]  
    diabeteStatus <- intvals[[2]]  
    ageVal <- intvals[[3]]
    
    # 80% chance to adhere the dynamic intervention 
    # set.seed(1)
    newdf[,`:=`(cond_sel=rbinom(.N, 1, 0.2),
                cond_met=((age1992+t0)<ageVal),
                ldlcLow=runif(.N,0.1,3.4),
                ldlcInt=runif(.N,0.1,2.6),
                ldlcHigh=runif(.N,0.1,1.8),
                lipidDrugStatus = 1)]
    
    # low-risk: LDL-C < 3.4 (IIa,B)
    newdf[riskStatus==riskStatusVal[1] & cond_sel==1 #& cond_met==TRUE 
          ,`:=`(LDLC = pmin(ldlcLow,get(intvar)),
                lipidDrugStatus = 1)]
    # moderate-risk: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[2] & cond_met==1 #& cond_met==TRUE
          ,`:=`(LDLC = pmin(ldlcInt,get(intvar)),
                lipidDrugStatus = 1)]
    # high-risk w/o diabetes: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus!=diabeteStatus &
              cond_sel==1 #& cond_met==TRUE 
          ,`:=`(LDLC = pmin(ldlcInt,get(intvar)),
                lipidDrugStatus = 1)]
    # high-risk w diabetes: LDLC < 1.8 or 50%+ reduction
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus==diabeteStatus &
              cond_sel==1 #& cond_met==TRUE
          ,`:=`(LDLC = pmax(ldlcHigh,get(intvar)*0.5),
                lipidDrugStatus = 1)]
}
csc2020_int.1b40 <- function(newdf, pool, intvar, intvals, time_name, t) {
    # p1014: Lipid-lowering target
    # intvals = list(riskStatus=c(1,2,3),diabeteRecdStatus=1)
    riskStatusVal <- intvals[[1]]  
    diabeteStatus <- intvals[[2]]  
    ageVal <- intvals[[3]]
    
    # 80% chance to adhere the dynamic intervention 
    # set.seed(1)
    newdf[,`:=`(cond_sel=rbinom(.N, 1, 0.4),
                cond_met=((age1992+t0)<ageVal),
                ldlcLow=runif(.N,0.1,3.4),
                ldlcInt=runif(.N,0.1,2.6),
                ldlcHigh=runif(.N,0.1,1.8),
                lipidDrugStatus = 1)]
    
    # low-risk: LDL-C < 3.4 (IIa,B)
    newdf[riskStatus==riskStatusVal[1] & cond_sel==1 #& cond_met==TRUE 
          ,`:=`(LDLC = pmin(ldlcLow,get(intvar)),
                lipidDrugStatus = 1)]
    # moderate-risk: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[2] & cond_met==1 #& cond_met==TRUE
          ,`:=`(LDLC = pmin(ldlcInt,get(intvar)),
                lipidDrugStatus = 1)]
    # high-risk w/o diabetes: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus!=diabeteStatus &
              cond_sel==1 #& cond_met==TRUE 
          ,`:=`(LDLC = pmin(ldlcInt,get(intvar)),
                lipidDrugStatus = 1)]
    # high-risk w diabetes: LDLC < 1.8 or 50%+ reduction
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus==diabeteStatus &
              cond_sel==1 #& cond_met==TRUE
          ,`:=`(LDLC = pmax(ldlcHigh,get(intvar)*0.5),
                lipidDrugStatus = 1)]
}
csc2020_int.1b60 <- function(newdf, pool, intvar, intvals, time_name, t) {
    # p1014: Lipid-lowering target
    # intvals = list(riskStatus=c(1,2,3),diabeteRecdStatus=1)
    riskStatusVal <- intvals[[1]]  
    diabeteStatus <- intvals[[2]]  
    ageVal <- intvals[[3]]
    
    # 80% chance to adhere the dynamic intervention 
    # set.seed(1)
    newdf[,`:=`(cond_sel=rbinom(.N, 1, 0.6),
                cond_met=((age1992+t0)<ageVal),
                ldlcLow=runif(.N,0.1,3.4),
                ldlcInt=runif(.N,0.1,2.6),
                ldlcHigh=runif(.N,0.1,1.8),
                lipidDrugStatus = 1)]
    
    # low-risk: LDL-C < 3.4 (IIa,B)
    newdf[riskStatus==riskStatusVal[1] & cond_sel==1 #& cond_met==TRUE 
          ,`:=`(LDLC = pmin(ldlcLow,get(intvar)),
                lipidDrugStatus = 1)]
    # moderate-risk: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[2] & cond_met==1 #& cond_met==TRUE
          ,`:=`(LDLC = pmin(ldlcInt,get(intvar)),
                lipidDrugStatus = 1)]
    # high-risk w/o diabetes: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus!=diabeteStatus &
              cond_sel==1 #& cond_met==TRUE 
          ,`:=`(LDLC = pmin(ldlcInt,get(intvar)),
                lipidDrugStatus = 1)]
    # high-risk w diabetes: LDLC < 1.8 or 50%+ reduction
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus==diabeteStatus &
              cond_sel==1 #& cond_met==TRUE
          ,`:=`(LDLC = pmax(ldlcHigh,get(intvar)*0.5),
                lipidDrugStatus = 1)]
}
# No.2 
# with(datPP[t0==0],table(riskStatus,LDLC>1.8))
# p.int.6 <- ggplot(datPP[t0==0],aes(x=LDLC,fill=riskStatus)) +
#     geom_histogram(position = position_dodge(0.3)) +
#     theme_bw()
# p.int.6
# datX <- datPP
# datPP.tmp <- datX[riskStatus==2
#                    ,`:=`(LDLC2 = ifelse(get("LDLC")<1.8,
#                                         get("LDLC"),
#                                         pmin(1.8,get("LDLC")*0.5))),
#                    by=.(ID)]
# datPP.tmp[,`:=`(LDLC=LDLC2)]
# ggplot(datPP.tmp[t0==0],aes(x=LDLC2,fill=riskStatus)) +
#     geom_histogram(position = position_dodge(0.3)) +
#     theme_bw()

esc2021_int.2a <- function(newdf, pool, intvar, intvals, time_name, t) {
    riskStatusVal <- intvals[[1]]
    ageVal <- intvals[[2]]
    
    # set.seed(1)
    newdf[,`:=`(cond_met=((age1992+t0)<ageVal),
                ldlcHigh=runif(.N,0.1,1.8),
                ldlcProp=runif(.N,0.5,0.99))]
    # high-risk with age < 70 years of age
    newdf[riskStatus %in% riskStatusVal[3] #& cond_met==TRUE 
          ,`:=`(#LDLC = pmin(ldlcHigh,get(intvar)*(1-ldlcProp))
                LDLC = get(intvar)*(1-ldlcProp),
                lipidDrugStatus = 1)]
}
esc2021_int.2b <- function(newdf, pool, intvar, intvals, time_name, t) {
    riskStatusVal <- intvals[[1]]
    ageVal <- intvals[[2]]
    
    # set.seed(1)
    newdf[,`:=`(cond_sel=rbinom(.N,1,0.8),
                cond_met=((age1992+t0)<ageVal),
                ldlcHigh=runif(.N,0.1,1.8),
                ldlcProp=runif(.N,0.5,0.99))]
    # high-risk with age < 70 years of age
    newdf[riskStatus %in% riskStatusVal[3] & cond_sel==1 # & cond_met==TRUE
          ,`:=`(#LDLC = pmin(ldlcHigh,get(intvar)*(1-ldlcProp))
              LDLC = get(intvar)*(1-ldlcProp),
              lipidDrugStatus = 1)]
}

# No.7
# p.int.7 <- ggplot(datPP[t0==0],aes(x=LDLC,fill=riskStatus)) +
#     geom_histogram(position=position_dodge(0.3)) +
#     theme_bw() +
#     facet_grid(~diabetesRecdStatus)
# p.int.7
aha2018_int.3a <- function(newdf, pool, intvar, intvals, time_name, t) {
    riskStatusVal <- intvals[[1]]
    diabeteVal <- intvals[[2]]
    ageVal <- intvals[[3]]
    
    # set.seed(1)
    newdf[,`:=`(ldlcPropLow=runif(.N,0.3,0.49),
                ldlcPropHigh=runif(.N,0.5,0.99),
                cond_ldlc4.9=(LDLC>=4.9),
                cond_ldlc1.74.8=(LDLC>=1.7 & LDLC<4.8),
                cond_age75=((age1992+t0)<ageVal),
                cond_age4075=((age1992+t0)>=40 & (age1992+t0)<=75))]
    # LDL-C >= 4.9 mmol/L
    newdf[cond_ldlc4.9==TRUE
          ,`:=`(LDLC =get(intvar)*(1-ldlcPropHigh),
                lipidDrugStatus = 1)]
    # Diabetes with age 40-75 y
    newdf[diabetesRecdStatus %in% diabeteVal[1] & cond_age4075==TRUE
          ,`:=`(LDLC = get(intvar)*(1-ldlcPropLow),
                lipidDrugStatus = 1)]
    # High-risk with age 40-75 y + 1.7 <= LDL-C <= 4.8
    newdf[riskStatus %in% riskStatusVal[3] & 
              cond_age4075==TRUE &
              cond_ldlc1.74.8==TRUE
          ,`:=`(LDLC = get(intvar)*(1-ldlcPropHigh),
                lipidDrugStatus = 1)]
    # Intermediate-risk with age 40-75 y
    newdf[riskStatus %in% riskStatusVal[2] & 
              cond_age4075==TRUE &
              cond_ldlc1.74.8==TRUE
          ,`:=`(LDLC = get(intvar)*(1-ldlcPropLow),
                lipidDrugStatus = 1)]
}

aha2018_int.3b <- function(newdf, pool, intvar, intvals, time_name, t) {
    riskStatusVal <- intvals[[1]]
    diabeteVal <- intvals[[2]]
    ageVal <- intvals[[3]]
    
    # set.seed(1)
    newdf[,`:=`(ldlcPropLow=runif(.N,0.3,0.49),
                ldlcPropHigh=runif(.N,0.5,0.99),
                cond_ldlc4.9=(LDLC>=4.9),
                cond_ldlc1.74.8=(LDLC>=1.7 & LDLC<4.8),
                cond_age75=((age1992+t0)<ageVal),
                cond_age4075=((age1992+t0)>=40 & (age1992+t0)<=75),
                cond_sel=rbinom(.N,1,0.8))]
    # LDL-C >= 4.9 mmol/L
    newdf[cond_ldlc4.9==TRUE & cond_sel==1
          ,`:=`(LDLC =get(intvar)*(1-ldlcPropHigh),
                lipidDrugStatus = 1)]
    # Diabetes with age 40-75 y
    newdf[diabetesRecdStatus %in% diabeteVal[1] & 
              cond_age4075==TRUE &
              cond_sel==1
          ,`:=`(LDLC = get(intvar)*(1-ldlcPropLow),
                lipidDrugStatus = 1)]
    # High-risk with age 40-75 y + 1.7 <= LDL-C <= 4.8
    newdf[riskStatus %in% riskStatusVal[3] & 
              cond_age4075==TRUE &
              cond_ldlc1.74.8==TRUE &
              cond_sel==1
          ,`:=`(LDLC = get(intvar)*(1-ldlcPropHigh),
                lipidDrugStatus = 1)]
    # Intermediate-risk with age 40-75 y
    newdf[riskStatus %in% riskStatusVal[2] & 
              cond_age4075==TRUE &
              cond_ldlc1.74.8==TRUE &
              cond_sel==1
          ,`:=`(LDLC = get(intvar)*(1-ldlcProphLow),
                lipidDrugStatus = 1)]
}
# No.8
# p.int.8 <- ggplot(datPP[t0==0],aes(x=nonHDL,fill=riskStatus)) +
#     geom_histogram(position = position_dodge(0.3)) +
#     theme_bw() +
#     facet_grid(~diabeteRecdStatus)
csc2020_int.4a <- function(newdf, pool, intvar, intvals, time_name, t) {
    # p1014: Lipid-lowering target
    # intvals = list(riskStatus=c(1,2,3),diabeteRecdStatus=1)
    riskStatusVal <- intvals[[1]]  
    diabeteStatus <- intvals[[2]]
    ageVal <- intvals[[3]]
    
    # set.seed(1)
    # Eligible age and ldlc
    newdf[,`:=`(cond_met=((age1992+t0)<ageVal),
                ldlcLow=runif(.N,0.1,3.4),
                ldlcInt=runif(.N,0.1,2.6),
                ldlcHigh=runif(.N,0.1,1.8),
                ldlcProp=runif(.N,0.5,0.99),
                nonHDLLow=runif(.N,0.1,4.2),
                nonHDLInt=runif(.N,0.1,3.4),
                nonHDLHigh=runif(.N,0.1,2.6))]
    
    # low-risk: LDL-C < 3.4 (IIa,B)
    newdf[riskStatus==riskStatusVal[1] # & cond_met==TRUE
          ,`:=`(LDLC = pmin(ldlcLow, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLLow),
                lipidDrugStatus = 1)]
    # moderate-risk: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[2] # & cond_met==TRUE
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w/o diabetes: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus!=diabeteStatus # & cond_met==TRUE
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w diabetes: LDLC < 1.8 or 50%+ reduction
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus==diabeteStatus # & cond_met==TRUE
          ,`:=`(LDLC = pmax(ldlcHigh, get(intvar)*(1-ldlcProp)),
                nonHDL = pmin(get("nonHDL"),nonHDLHigh),
                lipidDrugStatus = 1)]
}
csc2020_int.4b <- function(newdf, pool, intvar, intvals, time_name, t) {
    # p1014: Lipid-lowering target
    # intvals = list(riskStatus=c(1,2,3),diabeteRecdStatus=1)
    riskStatusVal <- intvals[[1]]  
    diabeteStatus <- intvals[[2]]
    ageVal <- intvals[[3]]
    
    # set.seed(1)
    # Eligible age and ldlc
    newdf[,`:=`(cond_met=((age1992+t0)<ageVal),
                ldlcLow=runif(.N,0.1,3.4),
                ldlcInt=runif(.N,0.1,2.6),
                ldlcHigh=runif(.N,0.1,1.8),
                ldlcProp=runif(.N,0.5,0.99),
                nonHDLLow=runif(.N,0.1,4.2),
                nonHDLInt=runif(.N,0.1,3.4),
                nonHDLHigh=runif(.N,0.1,2.6),
                cond_sel=rbinom(.N,1,0.8))]
    
    # low-risk: LDL-C < 3.4 (IIa,B)
    newdf[riskStatus==riskStatusVal[1] & 
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcLow, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLLow),
                lipidDrugStatus = 1)]
    # moderate-risk: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[2] & 
              # cond_met==TRUE & 
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w/o diabetes: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus!=diabeteStatus &
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w diabetes: LDLC < 1.8 or 50%+ reduction
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus==diabeteStatus & 
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmax(ldlcHigh, get(intvar)*(1-ldlcProp)),
                nonHDL = pmin(get("nonHDL"),nonHDLHigh),
                lipidDrugStatus = 1)]
}
csc2020_int.4b20 <- function(newdf, pool, intvar, intvals, time_name, t) {
    # p1014: Lipid-lowering target
    # intvals = list(riskStatus=c(1,2,3),diabeteRecdStatus=1)
    riskStatusVal <- intvals[[1]]  
    diabeteStatus <- intvals[[2]]
    ageVal <- intvals[[3]]
    
    # set.seed(1)
    # Eligible age and ldlc
    newdf[,`:=`(cond_met=((age1992+t0)<ageVal),
                ldlcLow=runif(.N,0.1,3.4),
                ldlcInt=runif(.N,0.1,2.6),
                ldlcHigh=runif(.N,0.1,1.8),
                ldlcProp=runif(.N,0.5,0.99),
                nonHDLLow=runif(.N,0.1,4.2),
                nonHDLInt=runif(.N,0.1,3.4),
                nonHDLHigh=runif(.N,0.1,2.6),
                cond_sel=rbinom(.N,1,0.2))]
    
    # low-risk: LDL-C < 3.4 (IIa,B)
    newdf[riskStatus==riskStatusVal[1] & 
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcLow, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLLow),
                lipidDrugStatus = 1)]
    # moderate-risk: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[2] & 
              # cond_met==TRUE & 
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w/o diabetes: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus!=diabeteStatus &
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w diabetes: LDLC < 1.8 or 50%+ reduction
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus==diabeteStatus & 
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmax(ldlcHigh, get(intvar)*(1-ldlcProp)),
                nonHDL = pmin(get("nonHDL"),nonHDLHigh),
                lipidDrugStatus = 1)]
}
csc2020_int.4b30 <- function(newdf, pool, intvar, intvals, time_name, t) {
    # p1014: Lipid-lowering target
    # intvals = list(riskStatus=c(1,2,3),diabeteRecdStatus=1)
    riskStatusVal <- intvals[[1]]  
    diabeteStatus <- intvals[[2]]
    ageVal <- intvals[[3]]
    
    # set.seed(1)
    # Eligible age and ldlc
    newdf[,`:=`(cond_met=((age1992+t0)<ageVal),
                ldlcLow=runif(.N,0.1,3.4),
                ldlcInt=runif(.N,0.1,2.6),
                ldlcHigh=runif(.N,0.1,1.8),
                ldlcProp=runif(.N,0.5,0.99),
                nonHDLLow=runif(.N,0.1,4.2),
                nonHDLInt=runif(.N,0.1,3.4),
                nonHDLHigh=runif(.N,0.1,2.6),
                cond_sel=rbinom(.N,1,0.3))]
    
    # low-risk: LDL-C < 3.4 (IIa,B)
    newdf[riskStatus==riskStatusVal[1] & 
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcLow, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLLow),
                lipidDrugStatus = 1)]
    # moderate-risk: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[2] & 
              # cond_met==TRUE & 
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w/o diabetes: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus!=diabeteStatus &
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w diabetes: LDLC < 1.8 or 50%+ reduction
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus==diabeteStatus & 
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmax(ldlcHigh, get(intvar)*(1-ldlcProp)),
                nonHDL = pmin(get("nonHDL"),nonHDLHigh),
                lipidDrugStatus = 1)]
}
csc2020_int.4b40 <- function(newdf, pool, intvar, intvals, time_name, t) {
    # p1014: Lipid-lowering target
    # intvals = list(riskStatus=c(1,2,3),diabeteRecdStatus=1)
    riskStatusVal <- intvals[[1]]  
    diabeteStatus <- intvals[[2]]
    ageVal <- intvals[[3]]
    
    # set.seed(1)
    # Eligible age and ldlc
    newdf[,`:=`(cond_met=((age1992+t0)<ageVal),
                ldlcLow=runif(.N,0.1,3.4),
                ldlcInt=runif(.N,0.1,2.6),
                ldlcHigh=runif(.N,0.1,1.8),
                ldlcProp=runif(.N,0.5,0.99),
                nonHDLLow=runif(.N,0.1,4.2),
                nonHDLInt=runif(.N,0.1,3.4),
                nonHDLHigh=runif(.N,0.1,2.6),
                cond_sel=rbinom(.N,1,0.4))]
    
    # low-risk: LDL-C < 3.4 (IIa,B)
    newdf[riskStatus==riskStatusVal[1] & 
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcLow, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLLow),
                lipidDrugStatus = 1)]
    # moderate-risk: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[2] & 
              # cond_met==TRUE & 
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w/o diabetes: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus!=diabeteStatus &
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w diabetes: LDLC < 1.8 or 50%+ reduction
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus==diabeteStatus & 
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmax(ldlcHigh, get(intvar)*(1-ldlcProp)),
                nonHDL = pmin(get("nonHDL"),nonHDLHigh),
                lipidDrugStatus = 1)]
}
csc2020_int.4b50 <- function(newdf, pool, intvar, intvals, time_name, t) {
    # p1014: Lipid-lowering target
    # intvals = list(riskStatus=c(1,2,3),diabeteRecdStatus=1)
    riskStatusVal <- intvals[[1]]  
    diabeteStatus <- intvals[[2]]
    ageVal <- intvals[[3]]
    
    # set.seed(1)
    # Eligible age and ldlc
    newdf[,`:=`(cond_met=((age1992+t0)<ageVal),
                ldlcLow=runif(.N,0.1,3.4),
                ldlcInt=runif(.N,0.1,2.6),
                ldlcHigh=runif(.N,0.1,1.8),
                ldlcProp=runif(.N,0.5,0.99),
                nonHDLLow=runif(.N,0.1,4.2),
                nonHDLInt=runif(.N,0.1,3.4),
                nonHDLHigh=runif(.N,0.1,2.6),
                cond_sel=rbinom(.N,1,0.5))]
    
    # low-risk: LDL-C < 3.4 (IIa,B)
    newdf[riskStatus==riskStatusVal[1] & 
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcLow, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLLow),
                lipidDrugStatus = 1)]
    # moderate-risk: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[2] & 
              # cond_met==TRUE & 
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w/o diabetes: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus!=diabeteStatus &
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w diabetes: LDLC < 1.8 or 50%+ reduction
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus==diabeteStatus & 
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmax(ldlcHigh, get(intvar)*(1-ldlcProp)),
                nonHDL = pmin(get("nonHDL"),nonHDLHigh),
                lipidDrugStatus = 1)]
}
csc2020_int.4b60 <- function(newdf, pool, intvar, intvals, time_name, t) {
    # p1014: Lipid-lowering target
    # intvals = list(riskStatus=c(1,2,3),diabeteRecdStatus=1)
    riskStatusVal <- intvals[[1]]  
    diabeteStatus <- intvals[[2]]
    ageVal <- intvals[[3]]
    
    # set.seed(1)
    # Eligible age and ldlc
    newdf[,`:=`(cond_met=((age1992+t0)<ageVal),
                ldlcLow=runif(.N,0.1,3.4),
                ldlcInt=runif(.N,0.1,2.6),
                ldlcHigh=runif(.N,0.1,1.8),
                ldlcProp=runif(.N,0.5,0.99),
                nonHDLLow=runif(.N,0.1,4.2),
                nonHDLInt=runif(.N,0.1,3.4),
                nonHDLHigh=runif(.N,0.1,2.6),
                cond_sel=rbinom(.N,1,0.6))]
    
    # low-risk: LDL-C < 3.4 (IIa,B)
    newdf[riskStatus==riskStatusVal[1] & 
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcLow, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLLow),
                lipidDrugStatus = 1)]
    # moderate-risk: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[2] & 
              # cond_met==TRUE & 
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w/o diabetes: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus!=diabeteStatus &
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w diabetes: LDLC < 1.8 or 50%+ reduction
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus==diabeteStatus & 
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmax(ldlcHigh, get(intvar)*(1-ldlcProp)),
                nonHDL = pmin(get("nonHDL"),nonHDLHigh),
                lipidDrugStatus = 1)]
}
csc2020_int.4b70 <- function(newdf, pool, intvar, intvals, time_name, t) {
    # p1014: Lipid-lowering target
    # intvals = list(riskStatus=c(1,2,3),diabeteRecdStatus=1)
    riskStatusVal <- intvals[[1]]  
    diabeteStatus <- intvals[[2]]
    ageVal <- intvals[[3]]
    
    # set.seed(1)
    # Eligible age and ldlc
    newdf[,`:=`(cond_met=((age1992+t0)<ageVal),
                ldlcLow=runif(.N,0.1,3.4),
                ldlcInt=runif(.N,0.1,2.6),
                ldlcHigh=runif(.N,0.1,1.8),
                ldlcProp=runif(.N,0.5,0.99),
                nonHDLLow=runif(.N,0.1,4.2),
                nonHDLInt=runif(.N,0.1,3.4),
                nonHDLHigh=runif(.N,0.1,2.6),
                cond_sel=rbinom(.N,1,0.7))]
    
    # low-risk: LDL-C < 3.4 (IIa,B)
    newdf[riskStatus==riskStatusVal[1] & 
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcLow, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLLow),
                lipidDrugStatus = 1)]
    # moderate-risk: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[2] & 
              # cond_met==TRUE & 
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w/o diabetes: LDLC < 2.6 (I,A)
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus!=diabeteStatus &
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmin(ldlcInt, get(intvar)),
                nonHDL = pmin(get("nonHDL"),nonHDLInt),
                lipidDrugStatus = 1)]
    # high-risk w diabetes: LDLC < 1.8 or 50%+ reduction
    newdf[riskStatus==riskStatusVal[3] & 
              diabetesRecdStatus==diabeteStatus & 
              # cond_met==TRUE &
              cond_sel==1
          ,`:=`(LDLC = pmax(ldlcHigh, get(intvar)*(1-ldlcProp)),
                nonHDL = pmin(get("nonHDL"),nonHDLHigh),
                lipidDrugStatus = 1)]
}
# No.9
esc2021_int.5a <- function(newdf, pool, intvar, intvals, time_name, t) {
    riskStatusVal <- intvals[[1]]
    ageVal <- intvals[[2]]
    
    # set.seed(1)
    newdf[,`:=`(cond_met=((age1992+t0)<ageVal),
                ldlcHigh=runif(.N,0.1,1.8),
                ldlcProb=runif(.N,0.5,0.99),
                nonHDLHigh=runif(.N,0.1,2.6))]
    # high-risk with age < 70 years of age
    newdf[riskStatus %in% riskStatusVal[3] # & cond_met==TRUE 
          ,`:=`(LDLC = pmin(ldlcHigh,get(intvar)*(1-ldlcProb)),
                nonHDL = pmin(get("nonHDL"),nonHDLHigh))]
}
esc2021_int.5b <- function(newdf, pool, intvar, intvals, time_name, t) {
    riskStatusVal <- intvals[[1]]
    ageVal <- intvals[[2]]
    
    # set.seed(1)
    newdf[,`:=`(cond_met=((age1992+t0)<ageVal),
                ldlcHigh=runif(.N,0.1,1.8),
                ldlcProb=runif(.N,0.5,0.99),
                nonHDLHigh=runif(.N,0.1,2.6),
                cond_sel=rbinom(.N,1,0.8))]
    # high-risk with age < 70 years of age
    newdf[riskStatus %in% riskStatusVal[3] & 
              # cond_met==TRUE & 
              cond_sel==1 
          ,`:=`(LDLC = pmin(ldlcHigh,get(intvar)*(1-ldlcProb)),
                nonHDL = pmin(get("nonHDL"),nonHDLHigh))]
}

## Version 1
## Old version
interventions.list <- 
    list(
        list(c(ldlc_int,0)),
        list(c(csc2020_int.1a,c(0,1,2),1,75)),
        list(c(csc2020_int.1b,c(0,1,2),1,75)),
        # list(c(esc2021_int.2a,c(0,1,2),70)),
        # list(c(esc2021_int.2b,c(0,1,2),70)),
        list(c(aha2018_int.3a,c(0,1,2),1,75)),
        list(c(aha2018_int.3b,c(0,1,2),1,75)),
        list(c(csc2020_int.4a,c(0,1,2),1,75)),
        list(c(csc2020_int.4b,c(0,1,2),1,75)))#,
        # list(c(esc2021_int.5a,c(0,1,2),70)),
        # list(c(esc2021_int.5b,c(0,1,2),70)))
interventions.descript <- c(
    "per 1 mmol/L reduction",
    "1a. 2020 CSC Guideline",
    "1b. 2020 CSC Guideline",
    # "2a. 2021 ESC Guideline",
    # "2b. 2021 ESC Guideline",
    "3a. 2018 AHA Guideline",
    "3b. 2018 AHA Guideline",
    "4a. 2020 CSC Guideline",
    "4b. 2020 CSC Guideline")#,
    # "5a. 2021 ESC Guideline",
    # "5b. 2021 ESC Guideline")
intvars.list <- as.list(rep("LDLC",length(interventions.list)))

## Version 2
## Update: 20231015
interventions.list <- 
    list(
        list(c(ldlc_int,0)),
        # list(c(csc2020_int.1a,c(0,1,2),1,75)),
        # list(c(csc2020_int.1b,c(0,1,2),1,75)),
        # list(c(esc2021_int.2a,c(0,1,2),70)),
        # list(c(esc2021_int.2b,c(0,1,2),70)),
        list(c(aha2018_int.3a,c(0,1,2),1,75)),
        list(c(aha2018_int.3b,c(0,1,2),1,75)),
        list(c(csc2020_int.4a,c(0,1,2),1,75)),
        list(c(csc2020_int.4b,c(0,1,2),1,75)))#,
# list(c(esc2021_int.5a,c(0,1,2),70)),
# list(c(esc2021_int.5b,c(0,1,2),70)))
interventions.descript <- c(
    "per 1 mmol/L reduction",
    # "1a. 2020 CSC Guideline",
    # "1b. 2020 CSC Guideline",
    # "2a. 2021 ESC Guideline",
    # "2b. 2021 ESC Guideline",
    "3a. 2018 AHA Guideline",
    "3b. 2018 AHA Guideline",
    "4a. 2020 CSC Guideline",
    "4b. 2020 CSC Guideline")#,
# "5a. 2021 ESC Guideline",
# "5b. 2021 ESC Guideline")
intvars.list <- as.list(rep("LDLC",length(interventions.list)))

## Version 3
## Update: 20231114
interventions.list.SA <- 
    list(list(c(csc2020_int.4b20,c(0,1,2),1,75)),
         list(c(csc2020_int.4b30,c(0,1,2),1,75)),
         list(c(csc2020_int.4b40,c(0,1,2),1,75)),
         list(c(csc2020_int.4b50,c(0,1,2),1,75)),
         list(c(csc2020_int.4b60,c(0,1,2),1,75)),
         list(c(csc2020_int.4b70,c(0,1,2),1,75)))

interventions.descript.SA <- c(
    "2020 CSC Guideline (20%)",
    "2020 CSC Guideline (30%)",
    "2020 CSC Guideline (40%)",
    "2020 CSC Guideline (50%)",
    "2020 CSC Guideline (60%)",
    "2020 CSC Guideline (70%)")#,
# "5a. 2021 ESC Guideline",
# "5b. 2021 ESC Guideline")
intvars.list.SA <- as.list(rep("LDLC",length(interventions.list.SA)))
