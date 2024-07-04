rm(list=ls())

a <- Sys.time()
cat(paste0("== ",Sys.time()[1],"...\n"))
set.seed(123456)

## intervention
pkgs <- c("data.table","gfoRmula","discSurv","ggplot2","parallel")
used_pkgs <- suppressPackageStartupMessages(inspkgs(pkgs))

# path.res
load(paste(path.cmcs,"/datPP-20240102.RData",sep=""))
datPP$riskGrpStatus1 <- (datPP$riskGrpStatus==1)+0
datPP$riskGrpStatus2 <- (datPP$riskGrpStatus==2)+0

## prepare baseline data set
datBaseline <- datPP[t0==0,c(id,W),with=FALSE]
datBaseline <- data.table::setDT(datBaseline)
datBaseline[,`:=`(lex.id=get(id))]
## prepare outcome data set
datOutcome <- data.table::setDT(datPP[t0==0,c(id,doBirth,doExam,outY,dooutY,outYLabel),with=FALSE])
## correct intermediate outcome using the end of study period
id.alive <- datOutcome$ID[datOutcome[[outYLabel]]=="Alive" & datOutcome[[dooutY]]<as.Date("2020-12-31")]
datIntermediateOutcome <- datOutcome[,c(id,dooutY),with=FALSE]
datIntermediateOutcome <- datIntermediateOutcome[,c(id,dooutY),with=FALSE]
datIntermediateOutcome[,`:=`(doIntermediateY=get(dooutY),intermediateY=1)]
datIntermediateOutcome[,(dooutY):=(NULL)]
datOutcome[[dooutY]][datOutcome$ID %in% id.alive] <- as.Date("2020-12-31")
## prepare time-varying variables 
datTimevary <- datPP[!is.na(datPP$waveStatusLabel),c(id,visit,A,L,doExam),with=FALSE]

cat(paste0("== ",Sys.time()[1],"...\n"))
suppressWarnings(
    datWideNew <- dat4LMTP(
        datBaseline=datBaseline,datOutcome=datOutcome,
        datTimevary=datTimevary,datIntermediateOutcome=NULL,
        id=id,W=W,L=L,A=A,ATimevary=TRUE,
        doBirth=doBirth,doExam=doExam,
        dooutY=dooutY,outY=outY,outYLabel=outYLabel,
        eoiYValue=eoiYValue,crsYValue=crsYValue,
        doIntermediateY=NULL) 
)
datWide <- datWideNew
## setting parameters
AWide <- grep(A,names(datWide),value=TRUE); # AWide
tau <- length(AWide); # tau
WWide <- W; # WWide
YWide <- grep("eoiY",names(datWide),value=TRUE); # YWide
CWide <- grep("censorY",names(datWide),value=TRUE); # CWide
LWide <- lapply(fix4time(0:(tau-1)),function(i) {
    Lx <- setdiff(grep(i,names(datWide),value=TRUE),c(AWide,CWide,YWide))
    Lx <- Lx[!grepl("riskGrpStatus",Lx)]
    Lx
}); # LWide


trim <- 0.995
folds <- 2
SL_folds <- 2
k <- 2
lrnrs <- c("SL.glm","SL.step","SL.mean","SL.glmnet")
# lrnrs <- c("SL.glm","SL.mean")

cat(paste0("== ",Sys.time()[1],"...\n"))
cat(paste0("== Create the shifted data...\n"))
# Step 1. manipulate the data set with hypothetical interventions
datWideint1      <- intLDL4DatNew(datWideInt=datWide,type="1mmol")
datWideint2      <- intLDL4DatNew(datWideInt=datWide,type="Statins")
datWideint20.8   <- intLDL4DatNew(datWideInt=datWide,type="Statins",adherenceRate=0.8)
datWideIntCSC    <- intCSC4Dat(datWideInt=datWide,adherenceRate=1.0)
datWideIntCSC0.8 <- intCSC4Dat(datWideInt=datWide,adherenceRate=0.8)
# step 2. change the censor status
datWideint1[,CWide] <- apply(datWideint1[,CWide],2,function(i) 1)
datWideint2[,CWide] <- apply(datWideint2[,CWide],2,function(i) 1)
datWideint20.8[,CWide] <- apply(datWideint20.8[,CWide],2,function(i) 1)
datWideIntCSC[,CWide] <- apply(datWideIntCSC[,CWide],2,function(i) 1)
datWideIntCSC0.8[,CWide] <- apply(datWideIntCSC0.8[,CWide],2,function(i) 1)



cat(paste0(Sys.time(),"\n"))
cat(paste0("== Estimate the subsequently doubly-robust estimators:\n",
           "   Step 6. Intervention Strategy 1. Natural course...\n"))
cl <- makeCluster(detectCores())
invisible(clusterEvalQ(cl, { library(lmtp) }))
clusterExport(cl,list("datWide","WWide","AWide","LWide","YWide","CWide",
                      "datWideint1","datWideint2","datWideint20.8",
                      "datWideIntCSC","datWideIntCSC0.8",
                      "trim","folds","SL_folds","k","lrnrs"),
              envir = environment())
resNC <- setNames(parLapply(
    cl,1:(tau-1),
    fun=function(i) {
        set.seed(i)
        lmtp::lmtp_sdr(
            data=datWide,
            baseline=WWide,
            trt=AWide[seq_len(i)],
            time_vary=LWide[seq_len(i)],
            outcome=YWide[seq_len(i)],
            outcome_type=ifelse(i==1,"binomial","survival"),
            id="lex.id",
            cens=CWide[seq_len(i)],
            shift=NULL,
            learners_trt=lrnrs,
            learners_outcome=lrnrs,
            folds=folds)
    }),1:(tau-1))
stopCluster(cl)
resNC[[28]]


cat(paste0("== ",Sys.time()[1],"...\n"))
cat(paste0("   Step 7. Intervention Strategy 2. Per 1 mmol/L LDL-C reduction...\n"))
cl <- makeCluster(detectCores())
invisible(clusterEvalQ(cl, { library(lmtp) }))
clusterExport(cl,list("datWide","WWide","AWide","LWide","YWide","CWide",
                      "datWideint1","datWideint2","datWideint20.8",
                      "datWideIntCSC","datWideIntCSC0.8",
                      "trim","folds","SL_folds","k","lrnrs"),
              envir = environment())
res1mmol <- setNames(parLapply(
    cl,1:(tau-1),
    fun=function(i) {
        set.seed(i)
        lmtp::lmtp_sdr(
            data=datWide,
            baseline=WWide,
            trt=AWide[seq_len(i)],
            time_vary=LWide[seq_len(i)],
            outcome=YWide[seq_len(i)],
            outcome_type=ifelse(i==1,"binomial","survival"),
            id="lex.id",
            cens=CWide[seq_len(i)],
            shifted=datWideint1,mtp=TRUE,
            learners_trt=lrnrs,
            learners_outcome=lrnrs,
            folds=folds)
    }),1:(tau-1))
stopCluster(cl)
res1mmol[[28]]


cat(paste0("== ",Sys.time()[1],"...\n"))
cat(paste0("   Step 8. Intervention Strategy 3. Statins therapeis...\n"))
cl <- makeCluster(detectCores())
invisible(clusterEvalQ(cl, { library(lmtp) }))
clusterExport(cl,list("datWide","WWide","AWide","LWide","YWide","CWide",
                      "datWideint1","datWideint2","datWideint20.8",
                      "datWideIntCSC","datWideIntCSC0.8",
                      "trim","folds","SL_folds","k","lrnrs"),
              envir = environment())
resStatin <- setNames(parLapply(
    cl,1:(tau-1),
    fun=function(i) {
        set.seed(i)
        lmtp::lmtp_sdr(
            data=datWide,
            baseline=WWide,
            trt=AWide[seq_len(i)],
            time_vary=LWide[seq_len(i)],
            outcome=YWide[seq_len(i)],
            outcome_type=ifelse(i==1,"binomial","survival"),
            id="lex.id",
            cens=CWide[seq_len(i)],
            shifted=datWideint2,mtp=TRUE,
            learners_trt=lrnrs,
            learners_outcome=lrnrs,
            folds=folds)
    }),1:(tau-1))
stopCluster(cl)
resStatin[[28]]


cat(paste0("== ",Sys.time()[1],"...\n"))
cat(paste0("   Step 9. Intervention Strategy 4. Statin therapies (80%)...\n"))
cl <- makeCluster(detectCores())
invisible(clusterEvalQ(cl, { library(lmtp) }))
clusterExport(cl,list("datWide","WWide","AWide","LWide","YWide","CWide",
                      "datWideint1","datWideint2","datWideint20.8",
                      "datWideIntCSC","datWideIntCSC0.8",
                      "trim","folds","SL_folds","k","lrnrs"),
              envir = environment())
resStatin0.8 <- setNames(parLapply(
    cl,1:(tau-1),
    fun=function(i) {
        set.seed(i)
        lmtp::lmtp_sdr(
            data=datWide,
            baseline=WWide,
            trt=AWide[seq_len(i)],
            time_vary=LWide[seq_len(i)],
            outcome=YWide[seq_len(i)],
            outcome_type=ifelse(i==1,"binomial","survival"),
            id="lex.id",
            cens=CWide[seq_len(i)],
            shifted=datWideint20.8,mtp=TRUE,
            learners_trt=lrnrs,
            learners_outcome=lrnrs,
            folds=folds)
    }),1:(tau-1))
stopCluster(cl)
resStatin0.8[[28]]


cat(paste0("== ",Sys.time()[1],"...\n"))
cat(paste0("   Step 10. Intervention Strategy 5. risk-based intervention...\n"))
cl <- makeCluster(detectCores())
invisible(clusterEvalQ(cl, { library(lmtp) }))
clusterExport(cl,list("datWide","WWide","AWide","LWide","YWide","CWide",
                      "datWideint1","datWideint2","datWideint20.8",
                      "datWideIntCSC","datWideIntCSC0.8",
                      "trim","folds","SL_folds","k","lrnrs"),
              envir = environment())
resCSC <- setNames(parLapply(
    cl,1:(tau-1),
    fun=function(i) {
        set.seed(i)
        lmtp::lmtp_sdr(
            data=datWide,
            baseline=WWide,
            trt=AWide[seq_len(i)],
            time_vary=LWide[seq_len(i)],
            outcome=YWide[seq_len(i)],
            outcome_type=ifelse(i==1,"binomial","survival"),
            id="lex.id",
            cens=CWide[seq_len(i)],
            shifted=datWideIntCSC,mtp=TRUE,
            learners_trt=lrnrs,
            learners_outcome=lrnrs,
            folds=folds)
    }),1:(tau-1))
stopCluster(cl)
resCSC[[28]]


cat(paste0("== ",Sys.time()[1],"...\n"))
cat(paste0("   Step 10. Intervention Strategy 5. risk-based intervention...\n"))
cl <- makeCluster(detectCores())
invisible(clusterEvalQ(cl, { library(lmtp) }))
clusterExport(cl,list("datWide","WWide","AWide","LWide","YWide","CWide",
                      "datWideint1","datWideint2","datWideint20.8",
                      "datWideIntCSC","datWideIntCSC0.8",
                      "trim","folds","SL_folds","k","lrnrs"),
              envir = environment())
resCSC0.8 <- setNames(parLapply(
    cl,1:(tau-1),
    fun=function(i) {
        set.seed(i)
        lmtp::lmtp_sdr(
            data=datWide,
            baseline=WWide,
            trt=AWide[seq_len(i)],
            time_vary=LWide[seq_len(i)],
            outcome=YWide[seq_len(i)],
            outcome_type=ifelse(i==1,"binomial","survival"),
            id="lex.id",
            cens=CWide[seq_len(i)],
            shifted=datWideIntCSC0.8,mtp=TRUE,
            learners_trt=lrnrs,
            learners_outcome=lrnrs,
            folds=folds)
    }),1:(tau-1))
stopCluster(cl)
resCSC0.8[[28]]
