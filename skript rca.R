library(rms)
library(ggplot2)
library(lattice)
library(grDevices)
library(DescTools)

# read data for real data example in Humberg et al 2020 Psychological Methods - 
# Cubic Response Surface Analysis,
# this file is partly based on the script illustration_Edwards1994_data.r
# accompanying this article (hereafter cubic script)

ed <- read.csv2("mba_csv.csv")
summary(ed) # 40 variables, no missing

# original scores as in cubic script line 95 and 97, but 
# in different variables
ed$MRACTo <- ed$MRACT + 20
ed$MRPREo <- ed$MRPRE + 20

# grand mean center the predictors (adapted from cubic script line 126-132)
gm <- mean(c(ed$MRACTo,ed$MRPREo), na.rm=TRUE)
ed$MRACT.c <- ed$MRACTo - gm
ed$MRPRE.c <- ed$MRPREo - gm
# outcome seems to be MRSAT (cubic script line 117 and 148)
summary(ed[,c("MRACTo","MRPREo","MRACT.c","MRPRE.c","MRSAT")])

# find convex hull for centered variables
chull.MRACT.c.MRPRE.c <- chull(ed[,c("MRACT.c","MRPRE.c")])
# extend to closed polygon
cfull.MRACT.c.MRPRE.c <- c(chull.MRACT.c.MRPRE.c, chull.MRACT.c.MRPRE.c[1])
plot(ed[,c("MRACT.c","MRPRE.c")], pch=19, cex=.5)
lines(ed[cfull.MRACT.c.MRPRE.c,c("MRACT.c","MRPRE.c")])
# check 
PtInPoly(pnts=ed[,c("MRACT.c","MRPRE.c")], 
         poly.pnts=ed[cfull.MRACT.c.MRPRE.c,c("MRACT.c","MRPRE.c")]) 
# all within, ok
PtInPoly(pnts=rbind(ed[,c("MRACT.c","MRPRE.c")],
                    data.frame(MRACT.c=c(-10,0,10),MRPRE.c=c(-20,-10,22))), 
         poly.pnts=ed[cfull.MRACT.c.MRPRE.c,c("MRACT.c","MRPRE.c")]) 
# seems ok

# trying models with restricted cubic splines and interaction
ddm <- datadist(ed[,c("MRACT.c","MRPRE.c")],adjto.cat='first')
options(datadist="ddm")
lm1 <- ols(MRSAT~rcs(MRACT.c,4) + rcs(MRPRE.c,4) + 
             rcs(MRACT.c,4) %ia% rcs(MRPRE.c,4), 
           data=ed, x=TRUE, y=TRUE)
lm2 <- ols(MRSAT~rcs(MRACT.c,4)*rcs(MRPRE.c,4), 
           data=ed, x=TRUE, y=TRUE)
lm1
lm2
anova(lm1)
anova(lm2)
summary(lm1)
summary(lm2)
ggplot(Predict(lm1))
ggplot(Predict(lm2))
AIC(lm1)
AIC(lm2)
# lm1 is better based on AIC
BIC(lm1)
BIC(lm2)
# lm1 is also better based on BIC
bplot(Predict(lm1), lfun=wireframe)
# empty diagram, somethingdon't work

# trying to make predictions in the congruence line MRACT.c=MRPRE.c
Predict.lm1 <- Predict(lm1)
summary(Predict.lm1)
summary(ed[,c("MRACT.c","MRPRE.c")])
# trying within the common range within MRACT.c=MRPRE.c
ddm
ddmc <- ddm
ddmc$limits[4:5,"MRACT.c"] # predikjonsgrenser, endrer
ddmc$limits[4:5,"MRACT.c"] <- c(-11,13)
ddmc$limits[4:5,"MRPRE.c"] <- c(-11,13)
cbind(ddm$limits, ddmc$limits) # ok

options(datadist="ddmc")
lm1c <- ols(MRSAT~rcs(MRACT.c,4) + rcs(MRPRE.c,4) + 
              rcs(MRACT.c,4) %ia% rcs(MRPRE.c,4), 
            data=ed, x=TRUE, y=TRUE)
AIC(lm1) - AIC(lm1c) # helt likt
BIC(lm1) - BIC(lm1c) # helt likt
lm1
lm1c # looks like the same
summary(lm1)
summary(lm1c)
Predict.lm1c <- Predict(lm1c)
summary(Predict.lm1c[,1:5] - Predict.lm1[,1:5]) # some differences, ok
summary(Predict.lm1c[,6] == Predict.lm1[,6]) # the same
ggplot(Predict.lm1)
ggplot(Predict.lm1c) # some differences
summary(Predict.lm1c$MRACT.c==Predict.lm1c$MRPRE.c) # never the same
summary(unique(Predict.lm1c$MRACT.c) %in% unique(Predict.lm1c$MRPRE.c))
summary(unique(Predict.lm1c$MRPRE.c) %in% unique(Predict.lm1c$MRACT.c))
# the same with one exception for both
# Predict does not "want" to be used for specific values of the 
# independent variables. Trying the simpler function predict
# first compare Predict and predict
checkPp <- Predict.lm1
checkPp$yhatp <- predict(lm1, newdata=checkPp[,c("MRACT.c","MRPRE.c")])
summary(checkPp)
summary(checkPp$yhat - checkPp$yhatp) # the same, ok
# also check the confidence bounds
checkPp <- Predict.lm1
checkPp$yhatp <- predict(lm1, newdata=checkPp[,c("MRACT.c","MRPRE.c")])
checkPp.int <- predict(lm1, se.fit=TRUE, interval="confidence",
                          newdata=checkPp[,c("MRACT.c","MRPRE.c")])
# strange message about intervals, trying via se
checkPp.se <- predict(lm1, se.fit=TRUE, 
                       newdata=checkPp[,c("MRACT.c","MRPRE.c")])
summary((checkPp$upper-checkPp$yhat)/checkPp.se$se.fit) # 1.975 always
dim(ed) # 172 observations, no missing
lm1 # 11 df, thus n-p-1=172-11-1=160
qt((1+.95)/2, df=160) # 1.9749, closer look
summary((checkPp$upper-checkPp$yhat)/checkPp.se$se.fit) # 1.975 always
summary((checkPp$upper-checkPp$yhat)/checkPp.se$se.fit - qt((1+.95)/2, df=160)) 
# differs in 15. decimal, knows what this is, 
# now start with the congruence line
summary(ed[,c("MRACT.c","MRPRE.c","PCSAT")])
# find common min and max for the independent variables
grandmin <- max(min(ed$MRACT.c), min(ed$MRPRE.c))
grandmax <- min(max(ed$MRACT.c), max(ed$MRPRE.c))
npoints <-200 # number of intermediate points
newdata.lm1.congr <- 
  data.frame(MRACT.c=seq(from=grandmin, to=grandmax,length.out=npoints),
            MRPRE.c=seq(from=grandmin, to=grandmax,length.out=npoints))
summary(newdata.lm1.congr)
summary(newdata.lm1.congr$MRACT.c - newdata.lm1.congr$MRPRE.c) # likt, ok
predict.lm1.congr <- predict(lm1, newdata=newdata.lm1.congr, se.fit=TRUE)
conf.use <- .95
qt.lm1 <- qt((1+conf.use)/2, df=lm1$df.residual)
preddata.lm1.congr <- 
  cbind(newdata.lm1.congr, pred=predict.lm1.congr$lin,
        lower=predict.lm1.congr$lin - qt.lm1*predict.lm1.congr$se.fit,
        upper=predict.lm1.congr$lin + qt.lm1*predict.lm1.congr$se.fit)
summary(preddata.lm1.congr)
summary(ed$PCSAT) # 0 to 37
plot(x=preddata.lm1.congr$MRACT.c, y=preddata.lm1.congr$pred,
     ylim=range(ed$PCSAT), las=1, type="l",
     xlab="centered MRCT=centered MRPRE, congruence line",
     ylab="predictd MRSAT (outcome)") 
points(x=preddata.lm1.congr$MRACT.c, y=preddata.lm1.congr$lower,
       lty="dotted",type="l")
points(x=preddata.lm1.congr$MRACT.c, y=preddata.lm1.congr$upper,
       lty="dotted", type="l") # reasonably flat, except to the
# right, where there is a bit up. Thus, not large differences along the
# congruence line, in the parts with much data. fairly high


# function for finding the combination with max predicted value at each
# line of incongruence, for the selected number of means for the independent
# variables. The function max.incongr.f assumes data with two independent 
# variables and one outcome, in that order, and no missing values.
# the model used is based on restricted cubic splines with 4 knots, by the 
# two independent variables (grand mean centered), without higher order
# interactions, using the R package rms
max.incongr.f <- function(data, model, npoints=200, conf=.95,
         xlegend="first independent variable", 
         ylegend="second independent variable", zlegend="outcome") {
  names(data) <- c("X","Y","Z")
  grandmean <- mean(c(data$X,data$Y))
  data$Xc <- data$X - grandmean
  data$Yc <- data$Y - grandmean
  grandmin <- max(min(data$Xc), min(data$Yc))
  grandmax <- min(max(data$Xc), max(data$Yc))
  
  mod <- ols(Z~rcs(Xc,4) + rcs(Yc,4) + rcs(Xc,4) %ia% rcs(Yc,4),
                data=data, x=TRUE, y=TRUE)
  qt.mod <- qt((1+conf)/2, df=mod$df.residual)
  # predicted along the congruence line
  newdata.mod.congr <-
    data.frame(Xc=seq(from=grandmin, to=grandmax, length.out=npoints),
               Yc=seq(from=grandmin, to=grandmax, length.out=npoints))
  predict.mod.congr <- predict(mod, newdata=newdata.mod.congr, se.fit=TRUE)
  preddata.mod.congr <- 
    cbind(newdata.mod.congr, pred=predict.mod.congr$lin,
          lower=predict.mod.congr$lin - qt.mod*predict.mod.congr$se.fit,
          upper=predict.mod.congr$lin + qt.mod*predict.mod.congr$se.fit)
  
  # predicted along the main incongruence line
  newdata.mod.main.incongr <-
    data.frame(Xc=seq(from=grandmin, to=grandmax, length.out=npoints),
               Yc=-seq(from=grandmin, to=grandmax, length.out=npoints))
  # restrict main incongruence line to Xc and Yc within grandmin and grandmax
  newdata.mod.main.incongr <- newdata.mod.main.incongr[
    newdata.mod.main.incongr$Xc>=grandmin&newdata.mod.main.incongr$Xc<=grandmax&
      newdata.mod.main.incongr$Yc>=grandmin&newdata.mod.main.incongr$Yc<=grandmax,]
  predict.mod.main.incongr <- predict(mod, newdata=newdata.mod.main.incongr, se.fit=TRUE)
  preddata.mod.main.incongr <-
    cbind(newdata.mod.main.incongr, pred=predict.mod.main.incongr$lin,
          lower=predict.mod.main.incongr$lin - qt.mod*predict.mod.main.incongr$se.fit,
          upper=predict.mod.main.incongr$lin + qt.mod*predict.mod.main.incongr$se.fit)
          
  return(list(data=data, mod=mod, 
              xlegend=xlegend, ylegend=ylegend, zlegend=zlegend,
              pred.congr=preddata.mod.congr, pred.main.incongr=preddata.mod.main.incongr))
} # end function max.incongr.f

# plotting functions for objects returned from max.incongr.f
# plot along congruence line
plot.congr <- function(max.incongr.obj) {
  data <- max.incongr.obj$data
  xlegend=max.incongr.obj$xlegend
  ylegend=max.incongr.obj$ylegend
  zlegend=max.incongr.obj$zlegend
  preddata <- max.incongr.obj$pred.congr
  plot(x=preddata$Xc, y=preddata$pred,
             ylim=range(data$Z), las=1, type="l",
             xlab=xlegend, ylab=zlegend)
  points(x=preddata$Xc, y=preddata$lower,
               lty="dotted", type="l")
  points(x=preddata$Xc, y=preddata$upper,
               lty="dotted", type="l")
  } # end plot.congr
# plot along the main incongruence line
plot.congr <- function(max.incongr.obj) {
  data <- max.incongr.obj$data
  xlegend=max.incongr.obj$xlegend
  ylegend=max.incongr.obj$ylegend
  zlegend=max.incongr.obj$zlegend
  preddata <- max.incongr.obj$pred.congr
  plot(x=preddata$Xc, y=preddata$pred,
       ylim=range(data$Z), las=1, type="l",
       xlab=xlegend, ylab=zlegend)
  points(x=preddata$Xc, y=preddata$lower,
         lty="dotted", type="l")
  points(x=preddata$Xc, y=preddata$upper,
         lty="dotted", type="l")
} # end plot.congr

# plot along main incongruence line
plot.main.incongr <- function(max.incongr.obj) {
  data <- max.incongr.obj$data
  xlegend=max.incongr.obj$xlegend
  ylegend=max.incongr.obj$ylegend
  zlegend=max.incongr.obj$zlegend
  preddata <- max.incongr.obj$pred.main.incongr
  plot(x=preddata$Xc, y=preddata$pred,
       ylim=range(data$Z), las=1, type="l",
       xlab=xlegend, ylab=zlegend)
  box()
  points(x=preddata$Xc, y=preddata$lower,
         lty="dotted", type="l")
  points(x=preddata$Xc, y=preddata$upper,
         lty="dotted", type="l")
} # end plot.main.incongr

max.incongr1 <- max.incongr.f((data=ed[,c("MRACT","MRPRE","MRSAT")])) 
summary(max.incongr1$data)
summary(max.incongr1$pred.congr)
summary(max.incongr1$pred.main.incongr)
plot.congr(max.incongr.obj=max.incongr1)
plot.main.incongr(max.incongr.obj=max.incongr1)
