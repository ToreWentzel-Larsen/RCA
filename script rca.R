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
# all within, ok, looks at slightly displaced data
PtInPoly(pnts=ed[,c("MRACT.c","MRPRE.c")]+ matrix(runif(2*dim(ed)[1],-.5,.5), ncol=2), 
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
# empty diagram, something doesn't work

# trying to make predictions in the congruence line MRACT.c=MRPRE.c
Predict.lm1 <- Predict(lm1)
summary(Predict.lm1)
summary(ed[,c("MRACT.c","MRPRE.c")])
# trying within the common range within MRACT.c=MRPRE.c
ddm
ddmc <- ddm
ddmc$limits[4:5,"MRACT.c"] # prediction bounds, change
ddmc$limits[4:5,"MRACT.c"] <- c(-11,13)
ddmc$limits[4:5,"MRPRE.c"] <- c(-11,13)
cbind(ddm$limits, ddmc$limits) # ok
# model
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
max.incongr.f <- function(data, model="rcsia", 
                          npoints=200, incongr.show.nr=100,
                          conf=.95, xlegend="first independent variable", 
                          ylegend="second independent variable", 
                          zlegend="outcome") {
  names(data) <- c("X","Y","Z")
  grandmean <- mean(c(data$X,data$Y))
  data$Xc <- data$X - grandmean
  data$Yc <- data$Y - grandmean
  grandmin <- max(min(data$Xc), min(data$Yc))
  grandmax <- min(max(data$Xc), max(data$Yc))
  # convex hull of Xc and Yc
  chull.c <- chull(data[,c("Xc","Yc")])
  # extend to closed polygon
  cfull.c <- c(chull.c, chull.c[1])
  # model
  mod <- ols(Z~rcs(Xc,4) + rcs(Yc,4) + rcs(Xc,4) %ia% rcs(Yc,4),
                data=data, x=TRUE, y=TRUE)
  qt.mod <- qt((1+conf)/2, df=mod$df.residual)
  # predicted along the congruence line
  newdata.mod.congr <-
    data.frame(Xc=seq(from=grandmin, to=grandmax, length.out=npoints),
               Yc=seq(from=grandmin, to=grandmax, length.out=npoints))
  inhull.mod.congr <- PtInPoly(pnts=newdata.mod.congr, 
           poly.pnts=newdata.mod.congr[cfull.c,])[,3]
  predict.mod.congr <- predict(mod, newdata=newdata.mod.congr, se.fit=TRUE)
  preddata.mod.congr <- 
    cbind(newdata.mod.congr, pred=predict.mod.congr$lin,
          lower=predict.mod.congr$lin - qt.mod*predict.mod.congr$se.fit,
          upper=predict.mod.congr$lin + qt.mod*predict.mod.congr$se.fit,
          inhull=inhull.mod.congr)
  # predicted along the main incongruence line
  newdata.mod.main.incongr <-
    data.frame(Xc=seq(from=grandmin, to=grandmax, length.out=npoints),
               Yc=-seq(from=grandmin, to=grandmax, length.out=npoints))
  # restrict main incongruence line to Xc and Yc within grandmin and grandmax
  newdata.mod.main.incongr <- newdata.mod.main.incongr[
    newdata.mod.main.incongr$Xc>=grandmin&newdata.mod.main.incongr$Xc<=grandmax&
      newdata.mod.main.incongr$Yc>=grandmin&newdata.mod.main.incongr$Yc<=grandmax,]
  inhull.mod.main.incongr <- PtInPoly(pnts=newdata.mod.main.incongr, 
                               poly.pnts=newdata.mod.main.incongr[cfull.MRACT.c.MRPRE.c,])[,3]
  predict.mod.main.incongr <- predict(mod, newdata=newdata.mod.main.incongr, se.fit=TRUE)
  preddata.mod.main.incongr <-
    cbind(newdata.mod.main.incongr, pred=predict.mod.main.incongr$lin,
          lower=predict.mod.main.incongr$lin - qt.mod*predict.mod.main.incongr$se.fit,
          upper=predict.mod.main.incongr$lin + qt.mod*predict.mod.main.incongr$se.fit,
          inhull=inhull.mod.main.incongr)
  # find maximum along each incongruence line
  max.incongr <- rep(NA,npoints)
  Xc.use <- seq(from=grandmin, to=grandmax, length.out=npoints)
  for (pointnr in 1:npoints) {
    Xc.point <- Xc.use[pointnr]
    Yc.point <- Xc.use[pointnr]
    Xc.point.incongr <- Xc.use
    Yc.point.incongr <- 2*Xc.point - Xc.point.incongr
    frame.point.incongr <- data.frame(Xc=Xc.point.incongr, Yc=Yc.point.incongr)
    inhull.point.incongr <- 
      PtInPoly(pnts=frame.point.incongr, 
               poly.pnts=frame.point.incongr[cfull.c,])[,3]
    if (sum(inhull.point.incongr)>1) {
      frame.point.incongr.inhull <- frame.point.incongr[inhull.point.incongr==1,]
      Xc.inhull.min <- min(frame.point.incongr.inhull$Xc)
      Xc.inhull.max <- max(frame.point.incongr.inhull$Xc)
      Xc.point.inhull.incongr <- seq(from=Xc.inhull.min, to=Xc.inhull.max, length.out=npoints)
      Yc.point.inhull.incongr <- 2*Xc.point - Xc.point.inhull.incongr
      newdata.point.inhull.incongr <- 
        data.frame(Xc=Xc.point.inhull.incongr, Yc=Yc.point.inhull.incongr)
      predict.mod.point.inhull.incongr <- 
        predict(mod, newdata=newdata.point.inhull.incongr)
      max.incongr[pointnr] <- max(predict.mod.point.inhull.incongr)
          if (pointnr==incongr.show.nr) {
            predict.mod.point.inhull.incongr.se <- 
              predict(mod, newdata=newdata.point.inhull.incongr, se.fit=TRUE)
            preddata.mod.shownr.inhull.incongr <-
              cbind(newdata.point.inhull.incongr, 
                    pred=predict.mod.point.inhull.incongr.se$lin,
                    lower=predict.mod.point.inhull.incongr.se$lin - 
                      qt.mod*predict.mod.point.inhull.incongr.se$se.fit,
                    upper=predict.mod.point.inhull.incongr.se$lin + 
                      qt.mod*predict.mod.point.inhull.incongr.se$se.fit)
          } # end if the incongruence line to show
    } # end if at least two points on the actual incongruence line within hull
  } # end points
  frame.points.incongr.max <- data.frame(Xc=Xc.use, max=max.incongr)
  frame.points.incongr.max <- 
    frame.points.incongr.max[is.na(frame.points.incongr.max$max)==0,]
  
  return(list(
    data=data, mod=mod, 
    xlegend=xlegend, ylegend=ylegend, zlegend=zlegend,
    pred.congr=preddata.mod.congr, pred.main.incongr=preddata.mod.main.incongr,
    frame.points.incongr.max=frame.points.incongr.max,
    preddata.incongr.show=preddata.mod.shownr.inhull.incongr,
    Xcshow=Xc.use[incongr.show.nr]))
} # end function max.incongr.f
max.incongr1 <- max.incongr.f((data=ed[,c("MRACT","MRPRE","MRSAT")])) # for check


# plotting functions for objects returned from max.incongr.f
# plot along congruence line, and along main inclongruence line
plot.congr.maini <- function(max.incongr.obj, lwdhull=3, col.congr="black", 
                             col.maini="red", col.maxi="blue", col.showi="orange") {
  data <- max.incongr.obj$data
  xlegend=max.incongr.obj$xlegend
  ylegend=max.incongr.obj$ylegend
  zlegend=max.incongr.obj$zlegend
  preddata.congr <- max.incongr.obj$pred.congr
  preddata.maini <- max.incongr.obj$pred.main.incongr
  preddata.showi <- max.incongr.obj$preddata.incongr.show
  maxframe <- max.incongr.obj$frame.points.incongr.max
  Xcshow=max.incongr.obj$Xcshow
  # congruence line
  plot(x=preddata.congr$Xc, y=preddata.congr$pred, las=1, type="l", 
       xlim=range(data$Xc), ylim=range(data$Z), 
       xlab=xlegend, ylab=zlegend, col=col.congr)
  points(x=preddata.congr$Xc[preddata.congr$inhull==1], 
         y=preddata.congr$pred[preddata.congr$inhull==1],
         lty="dotted", type="l", col=col.congr, lwd=lwdhull)
  points(x=preddata.congr$Xc, y=preddata.congr$lower,
         lty="dotted", type="l", col=col.congr)
  points(x=preddata.congr$Xc, y=preddata.congr$upper,
         lty="dotted", type="l", col=col.congr)
  # main incongruence line
  points(x=preddata.maini$Xc, 
         y=preddata.maini$pred,
         type="l",xlab=xlegend, ylab=zlegend, col=col.maini)
  points(x=preddata.maini$Xc[preddata.maini$inhull==1], 
         y=preddata.maini$pred[preddata.maini$inhull==1],
         type="l",xlab=xlegend, ylab=zlegend, col=col.maini, lwd=lwdhull)
  points(x=preddata.maini$Xc, y=preddata.maini$lower,
         lty="dotted", type="l", col=col.maini)
  points(x=preddata.maini$Xc, y=preddata.maini$upper,
         lty="dotted", type="l", col=col.maini)
  # max at incongruence lines
  points(x=maxframe$Xc, y=maxframe$max,
         type="l",xlab=xlegend, ylab=zlegend, col=col.maxi)
  # chosen incongruence line
  points(x=preddata.showi$Xc, 
         y=preddata.showi$pred,
         type="l",xlab=xlegend, ylab=zlegend, col=col.showi)
  points(x=preddata.showi$Xc, y=preddata.showi$lower,
         lty="dotted", type="l", col=col.showi)
  points(x=preddata.showi$Xc, y=preddata.showi$upper,
         lty="dotted", type="l", col=col.showi)
  # mark center
  abline(v=0, lty="dotted")
  abline(v=Xcshow, lty="dotted", col=col.showi)
} # end plot.congr.maini

max.incongr1 <- max.incongr.f(data=ed[,c("MRACT","MRPRE","MRSAT")])
summary(max.incongr1$data)
summary(max.incongr1$pred.congr)
summary(max.incongr1$pred.main.incongr)
dim(max.incongr1$frame.points.incongr.max) # all 200 included
summary(max.incongr1$frame.points.incongr.max)
summary(max.incongr1$preddata.incongr.show)
max.incongr1$Xcshow
plot.congr.maini(max.incongr.obj=max.incongr1)
# with other incongruence line to show
max.incongr50 <- max.incongr.f(data=ed[,c("MRACT","MRPRE","MRSAT")],
                               incongr.show.nr=50)
max.incongr50$Xcshow
plot.congr.maini(max.incongr.obj=max.incongr50)
max.incongr150 <- max.incongr.f(data=ed[,c("MRACT","MRPRE","MRSAT")],
                               incongr.show.nr=150)
plot.congr.maini(max.incongr.obj=max.incongr150)
max.incongr175 <- max.incongr.f(data=ed[,c("MRACT","MRPRE","MRSAT")],
                                incongr.show.nr=175)
plot.congr.maini(max.incongr.obj=max.incongr175)
max.incongr190 <- max.incongr.f(data=ed[,c("MRACT","MRPRE","MRSAT")],
                                incongr.show.nr=190)
plot.congr.maini(max.incongr.obj=max.incongr190)
