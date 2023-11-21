## Functions 
area <- function(age,ldlc) sum(diff(age)*(ldlc[-length(ldlc)] + ldlc[-1])*0.5)
## load data files from SPSS
loadSPSS <- function(name,path) {
    data.table::as.data.table(foreign::read.spss(
        paste(path,name,sep=""),
        trim_values=TRUE,
        stringsAsFactors=FALSE,
        add.undeclared.levels="no"))
}


## get CIF from gformula objective
gformCIF <- function(allCauseMortality,CVD,ASCVD,intDescript=NULL) {
    ## Combined the cumulative risk
    datRes <- rbind(
        allCauseMortality$result[,`:=`(Outcome="All-cause mortality")],
        CVD$result[,`:=`(Outcome="Cardiovascular disease")],
        ASCVD$result[,`:=`(Outcome="Atherosclerostic cardiovascular disease")])
    if ('Risk lower 95% CI' %in% names(datRes)) {
        datResNew <- datRes[
            ,`:=`('Cumulative risk'=get('g-form risk'),
                  'Lower 95% CI'=get('Risk lower 95% CI'),
                  'Upper 95% CI'=get('Risk upper 95% CI'))]
        # datResNew <- melt(datRes,
        #                id.vars=c("k","Outcome","Interv."),
        #                measure.vars=c("NP Risk","g-form risk","Risk lower 95% CI","Risk upper 95% CI"),
        #                variable.name="Method",
        #                value.name="Cumulative risk")
        # datResNew <- datResNew[
        #     ,`:=`(Model=ifelse(Method=="NP Risk",
        #                        "Non-parametric estimate",
        #                 ifelse(Method=="g-form risk",
        #                        "Parametric g-formula estimate",
        #                 ifelse(Method=="Risk lower 95% CI",
                               # "Lower 95% CI","Upper 95%"))))]
        datResNew$Intervention <- with(
            datResNew,
            factor(`Interv.`,
                   labels = c("Natural course",intDescript),
                   levels = unique(`Interv.`)))
        datRes0 <- datResNew[k==0][,`:=`(`Cumulative risk`=0,t0=0,
                                         `Lower 95% CI`=0,
                                         `Upper 95% CI`=0)]
        datRes1 <- datResNew[,`:=`(t0=k+1)]
    } else {
        datResNew <- melt(datRes,
                       id.vars=c("k","Outcome"),
                       measure.vars=c("NP Risk","g-form risk"),
                       variable.name="Method",
                       value.name="Cumulative risk")
        datResNew <- datResNew[,`:=`(Model=ifelse(Method=="NP Risk",
                                                  "Non-parametric estimate",
                                                  "Parametric g-formula estimate"))]
        datRes0 <- datResNew[k==0][,`:=`(`Cumulative risk`=0,t0=0)]
        datRes1 <- datResNew[,`:=`(t0=k+1)]
    }

    datRes <- rbind(datRes0,datRes1)
    datRes <- datRes[order(Outcome,t0)]
    return(datRes)
}

getTableForest <- function(datCIF) {
    datRes <- datCIF[,.(`g-form risk` = last(`g-form risk`),
                        `Risk lower 95% CI` = last(`Risk lower 95% CI`),
                        `Risk upper 95% CI` = last(`Risk upper 95% CI`),
                        `Risk difference` = last(`Risk difference`),
                        `RD lower 95% CI` = last(`RD lower 95% CI`),
                        `RD upper 95% CI` = last(`RD upper 95% CI`),
                        `Risk ratio` = last(`Risk ratio`),
                        `RR lower 95% CI` = last(`RR lower 95% CI`),
                        `RR upper 95% CI` = last(`RR upper 95% CI`),
                        `% Intervened On` = last(`% Intervened On`),
                        `Aver % Intervened On` = last(`Aver % Intervened On`),
                        RMST=29-area(age=k,ldlc=`g-form risk`)),
                     by=.(Outcome,Intervention)]
    datRes <- datRes[
        ,`:=`(Risk=paste(format95CIOneline(`g-form risk`,
                                    `Risk lower 95% CI`,
                                    `Risk upper 95% CI`)),
              RiskDiff=paste(format95CIOneline(`Risk difference`,
                                        `RD lower 95% CI`,
                                        `RD upper 95% CI`)),
              RiskRatio=paste(format95CIOneline(`Risk ratio`,
                                         `RR lower 95% CI`,
                                         `RR upper 95% CI`)),
              RMET=sprintf("%.2f",RMST),
              cumIntervened=sprintf("%.2f",`% Intervened On`),
              avgIntervened=sprintf("%.2f",`Aver % Intervened On`),
              NNT=ceiling(first(RMST)/(RMST-first(RMST))))
        ,by=.(Outcome)]
    datRes$Outcome <- with(
        datRes,ifelse(Outcome=="Cardiovascular disease","CVD",
                      ifelse(Outcome=="All-cause mortality","All-cause mortality",
                             "ASCVD")))
    datResTable <- dcast(
        datRes,
        Intervention ~ Outcome,
        value.var = c("Risk","RiskDiff",
                      "Risk difference","RD lower 95% CI","RD upper 95% CI",
                      "RiskRatio","RMET",
                      "cumIntervened","avgIntervened","NNT"))
    datResTable$idx <- with(
        datResTable,ifelse(grepl("Natural",Intervention),1,
                           ifelse(grepl("Dynamic",Intervention) & grepl("CSC",Intervention), 2,
                                  ifelse(grepl("Feasible",Intervention) & grepl("CSC",Intervention),3,
                                         ifelse(grepl("Dynamic",Intervention) & grepl("AHA",Intervention),4,
                                                ifelse(grepl("Feasible",Intervention) & grepl("AHA",Intervention),5,6))))))
    datResTable <- datResTable[order(idx)]
    names <- names(datResTable)
    vnames <- c("Intervention",
                names[grepl("_CVD",names)],
                names[grepl("All-cause",names)],
                names[grepl("ASCVD",names)])
    return(datResTable[,vnames,with=FALSE])
}

getTable <- function(datCIF) {
    datRes <- datCIF[,.(`g-form risk` = last(`g-form risk`),
                        `Risk lower 95% CI` = last(`Risk lower 95% CI`),
                        `Risk upper 95% CI` = last(`Risk upper 95% CI`),
                        `Risk difference` = last(`Risk difference`),
                        `RD lower 95% CI` = last(`RD lower 95% CI`),
                        `RD upper 95% CI` = last(`RD upper 95% CI`),
                        `Risk ratio` = last(`Risk ratio`),
                        `RR lower 95% CI` = last(`RR lower 95% CI`),
                        `RR upper 95% CI` = last(`RR upper 95% CI`),
                        `% Intervened On` = last(`% Intervened On`),
                        `Aver % Intervened On` = last(`Aver % Intervened On`),
                        RMST=29-area(age=k,ldlc=`g-form risk`)),
                     by=.(Outcome,Intervention)]
    datRes <- datRes[
        ,`:=`(Risk=paste(format95CI(`g-form risk`,
                                    `Risk lower 95% CI`,
                                    `Risk upper 95% CI`)),
              RiskDiff=paste(format95CI(`Risk difference`,
                                        `RD lower 95% CI`,
                                        `RD upper 95% CI`)),
              RiskRatio=paste(format95CI(`Risk ratio`,
                                         `RR lower 95% CI`,
                                         `RR upper 95% CI`)),
              RMET=sprintf("%.2f",RMST),
              cumIntervened=sprintf("%.2f",`% Intervened On`),
              avgIntervened=sprintf("%.2f",`Aver % Intervened On`),
              NNT=ceiling(first(RMST)/(RMST-first(RMST))))
        ,by=.(Outcome)]
    datRes$Outcome <- with(
        datRes,ifelse(Outcome=="Cardiovascular disease","CVD",
                      ifelse(Outcome=="All-cause mortality","All-cause mortality",
                             "ASCVD")))
    datResTable <- dcast(
        datRes,
        Intervention ~ Outcome,
        value.var = c("Risk","RiskDiff","RiskRatio","RMET","cumIntervened","avgIntervened","NNT"))
    datResTable$idx <- with(
        datResTable,ifelse(grepl("Natural",Intervention),1,
                    ifelse(grepl("Dynamic",Intervention) & grepl("CSC",Intervention), 2,
                    ifelse(grepl("Feasible",Intervention) & grepl("CSC",Intervention),3,
                    ifelse(grepl("Dynamic",Intervention) & grepl("AHA",Intervention),4,
                    ifelse(grepl("Feasible",Intervention) & grepl("AHA",Intervention),5,6))))))
    datResTable <- datResTable[order(idx)]
    names <- names(datResTable)
    vnames <- c("Intervention",
                names[grepl("_CVD",names)],
                names[grepl("All-cause",names)],
                names[grepl("ASCVD",names)])
    return(datResTable[,vnames,with=FALSE])
}

# datCIF$Intervention <- with(datCIF, gsub(" (2020 CSC)","",as.character(Intervention)))
pltCIF <- function(datCIF) {
    outList <- unique(datCIF$Outcome)
    if('Risk lower 95% CI' %in% names(datCIF)) {
        pCIF <- setNames(lapply(outList, function(i) {
            ggplot(
                data=datCIF[Outcome %in% i],
                aes(x=t0,y=`Cumulative risk`,color=Intervention,fill=Intervention)) +
                geom_smooth(se=FALSE,linewidth=0.5) +
                geom_ribbon(aes(ymin=`Lower 95% CI`,ymax=`Upper 95% CI`,
                                color=Intervention,fill=Intervention),
                            alpha=0.3, linetype=0) + 
                theme_classic() +
                theme(axis.line.x=element_line(colour="black",linewidth=0.75),
                      axis.line.y=element_line(colour="black",linewidth=0.75),
                      panel.grid.major.y = element_line(colour="gray90"),
                      legend.position=c(0.15,0.865),
                      legend.title=element_blank()) +
                scale_color_jama() + scale_fill_jama() +
                scale_x_continuous(name="Follow-up time (years)",
                                   breaks=seq(0,30,5),
                                   limits=c(0,30)) +
                scale_y_continuous(name="Cumulative risk",limits=c(0,0.31)) +
                labs(title=i)
        }),outList)
    } else {
        pCIF <- setNames(lapply(outList, function(i) {
            ggplot(
                data=datCIF[Outcome %in% i],
                aes(x=t0,y=`Cumulative risk`,color=Model)) +
                geom_smooth(se=FALSE,linewidth=0.5) +
                theme_classic() +
                theme(axis.line.x=element_line(colour="black",linewidth=0.75),
                      axis.line.y=element_line(colour="black",linewidth=0.75),
                      panel.grid.major.y = element_line(colour="gray90"),
                      legend.position=c(0.15,0.865),
                      legend.title=element_blank()) +
                scale_color_jama() +
                scale_x_continuous(name="Follow-up time (years)",
                                   breaks=seq(0,30,5),
                                   limits=c(0,30)) +
                scale_y_continuous(name="Cumulative risk",limits=c(0,0.30)) +
                labs(title=i)
        }),outList)        
    }

    return(pCIF)
}



## save results from excel
saveExcel <- function(x,name,path) {
    xlsx::write.xlsx(x,file=paste(path,"/",name,".xlsx",sep=""))
}
saveCSV <- function(x,name,path) {
    write.csv(x,file=paste(path,"/",name,".csv",sep=""))
}

## plot cumulative risk

## save figure to pdf file 
pltSave <- function(list.plot,name,path,width=10,height=8) {
    ggsave(filename=paste(path,"/",name,".pdf",sep=""),
           plot=list.plot,width=width,height=height,units="in")
}

## Point pgform-estimated risk
gformNNT <- function(gformula_survival.fit,int_descript=int_descript) {
    datRes <- as.data.table(gformula_survival.fit$result)
    datRes$Intervention <- with(datRes, 
                                factor(Interv.,
                                       levels=sort(unique(Interv.)),
                                       labels=c("Natural course",int_descript)))
    datRes$Intervention <- relevel(datRes$Intervention,ref=1)
    datResNNT <- datRes[,.(risk=last(`g-form risk`),
                           riskRatio=last(`Risk ratio`),
                           riskDiff=last(`Risk difference`),
                           riskNNT=1/last(`Risk difference`),
                           RMST=29-area(age=k,ldlc=`g-form risk`),
                           NPrisk=last(`NP Risk`),
                           t0=last(k)+1)
                        ,by=.(Intervention)]
    datResNNT[,`:=`(RMSTRatio=RMST/(first(RMST)),
                    RMSTDiff=RMST-first(RMST),
                    RMSTNNT=first(RMST)/(RMST-first(RMST)))]
    return(list(datRes=datRes,datResNNT=datResNNT))
}

## format p-value
formatP <- function(x)  ifelse(x < 0.001, sprintf("%.3f",x),"<.001") 

## construct 95% confidence interval
format95CI <- function(x,lci,uci) {
    paste(sprintf("%.2f",x), "\n (",
          sprintf("%.2f",lci), " to ",
          sprintf("%.2f",uci), ")", sep="")
}

format95CIOneline <- function(x,lci,uci) {
    paste(sprintf("%.2f",x), " (",
          sprintf("%.2f",lci), " to ",
          sprintf("%.2f",uci), ")", sep="")
}

format95CINNT <- function(x,lci,uci) {
    paste(ceiling(x), " \n (",
          ceiling(lci), " to ",
          ceiling(uci),")", sep="")
}

## rename the intervention label
splitCom <- function(x) {
    a <- unlist(strsplit(as.character(x)," "))
    return(paste(a[2:length(a)],collapse=" ",sep=""))
}
