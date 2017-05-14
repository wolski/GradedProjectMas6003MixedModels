## ----loadpackage, setup, include=FALSE-----------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(lme4)
library(reshape2)


## ----include=FALSE-------------------------------------------------------

rm(list=ls())
source("multiplot.R")

load('OllyBolly.rda', envir = EnvOlly <- new.env())
OllyBolly <- EnvOlly$OllyBolly
rm(EnvOlly)


## ----firstplot-----------------------------------------------------------
g1 <- ggplot(OllyBolly, aes(x = age , y=height)) + geom_line(aes(colour=Seed)) + 
  theme(legend.position="none") + labs(title="A")


## ----secondplot----------------------------------------------------------
g2 <- ggplot(OllyBolly, aes(x = age , y=height)) + geom_line(aes(colour=Seed))+ 
  scale_x_log10()+ theme(legend.position="none") + scale_y_sqrt() + labs(title="B")


## ------------------------------------------------------------------------
OllyBolly$logAge <- log(OllyBolly$age)
OllyBolly$sqrtAge <- sqrt(OllyBolly$age)
OllyBolly$sqrtHeight <- sqrt(OllyBolly$height)

## ------------------------------------------------------------------------
OllyBolly <- subset(OllyBolly, age > 4)


## ---- fig.width=8, fig.height=4 ,fig.cap="Tree height versus age (A), square root of height versus log10(age)."----

multiplot(g1,g2,cols=2)


## ------------------------------------------------------------------------
lme2 <- lmer(sqrtHeight ~ 1 + logAge + (1|Seed) ,
             data=OllyBolly, REML = F)

## ------------------------------------------------------------------------
lme3 <- lmer(sqrtHeight ~ 1 + logAge + sqrtAge +  (1|Seed) ,
             data=OllyBolly, REML = F)


## ------------------------------------------------------------------------
aovres <- anova(lme2, lme3)
obs.test.stat<- - 2*(logLik(lme2)-logLik(lme3)) 
1 - pchisq(obs.test.stat,1)


## ----results="hide"------------------------------------------------------
lme2R <- lmer(sqrtHeight ~ 1 + logAge + sqrtAge + (1|Seed)  ,
             data=OllyBolly, REML = F)
lme3R <- lmer(sqrtHeight ~ 1 + logAge + sqrtAge + (1|Seed) + (logAge - 1|Seed)  ,
             data=OllyBolly, REML = F)
lme4R <- lmer(sqrtHeight ~ 1 + logAge + sqrtAge + (1|Seed) + (logAge - 1|Seed) + 
                (sqrtAge - 1|Seed) , data=OllyBolly, REML = F)
anovas <- anova(lme2R, lme3R, lme4R)

## ------------------------------------------------------------------------
lme5R <- lmer(sqrtHeight ~ 1 + logAge + sqrtAge +  (1 + logAge|Seed) ,
             data=OllyBolly, REML = F)
anovasres <- anova(lme3R, lme5R)
obs.test.stat<- - 2*(logLik(lme3R)-logLik(lme5R)) 
1 - pchisq(obs.test.stat,1)

## ------------------------------------------------------------------------

lmeChosen <- update(lme5R,  REML = T)
summary(lmeChosen)


## ----fig.cap = "Diagnostic model plots.", fig.width=8--------------------

raneffVar <- VarCorr(lmeChosen)
stdevs <-attributes(((raneffVar))$Seed)$stddev
par(mfrow=c(2,3))

randomEffects <- ranef(lmeChosen)[[1]]
qqnorm(randomEffects[,1],main=names(randomEffects)[1],main="A")
abline(c(0,stdevs[1]),col=2)

qqnorm(randomEffects[,2],main=names(randomEffects)[2],main="B")
abline(c(0,stdevs[2]),col=2)

plot(fitted(lmeChosen), resid(lmeChosen), main="C")# + geom_jitter()

qqnorm(resid(lmeChosen), main="D")
abline(c(0,(summary(lmeChosen))$sigma),col=2)

fitted.level0<-lmeChosen@pp$X %*% fixef(lmeChosen)
resid.level0<-OllyBolly$sqrtHeight-fitted.level0
plot(jitter(fitted.level0),resid.level0, main="E")

plot(fitted(lmeChosen),OllyBolly$sqrtHeight, main="F")
abline(0,1,col=2)


## ---- fig.cap="Panel A - qqnorm plot of model residues. Panel B - Level 0 resdiues vs fitted. Panel C - Measured height vs fitted values", fig.width=8----
par(mfrow=c(1,3))


## ---- fig.cap="Histogram of the simulated heights at 25 weeks. The red dots represent the actually measured heights of the trees."----
new.height<-simulate(lmeChosen, nsim=1000)^2
new.height.25 <- new.height[OllyBolly$age==25,]
hist(unlist(new.height.25), main="")

Olb25<-subset(OllyBolly, age==25)
points(Olb25$height, rep(1000, nrow(Olb25)), col=2, type="p")


## ------------------------------------------------------------------------
(pred25 <- fixef(lmeChosen)[1] + fixef(lmeChosen)[2] *log(25) + fixef(lmeChosen)[3] * 25)^2


## ------------------------------------------------------------------------
new.height<-simulate(lmeChosen, nsim=10000)^2

## ------------------------------------------------------------------------
new.height.10 <- new.height[OllyBolly$age==10,]
height.10_2.6 <- which(new.height.10 > (2.59) & new.height.10 < (2.615),arr.ind = T)

## ----Selecting25given10 , fig.cap="Left panel, histogram of heights at age 25 (black) and histogram of heights at age 25 for those seeds which had a height of 2.6 at the age of 10 (red). Right panel, Density plot of the same data."----
new.height.25 <- new.height[OllyBolly$age==25,]
height.25_height2.6At10 <- new.height.25[height.10_2.6]
par(mfrow=c(1,2))
hist(unlist(new.height.25), main="",ylab="height")
hist(height.25_height2.6At10, add=T , col=2,main="")

plot(density(unlist(new.height.25)), ylim=c(0,4), main="",ylab="height")
lines(density(height.25_height2.6At10),col=2)

## ------------------------------------------------------------------------
mean(new.height.25[height.10_2.6])


## ------------------------------------------------------------------------
quantile(new.height.25[height.10_2.6], c(0.025,0.975))


## ---- warning=FALSE------------------------------------------------------
N<-250
boot.test.stats<-rep(0,N)
for(i in 1:N){
  if(i %% 100 == 0){
    print(i)
    }
  new.height<-unlist(simulate(lme2))
  fm.reduced.new <-update(lme2, new.height ~ . )
  fm.full.new <- update(lme3, new.height ~ .)
  boot.test.stats[i]<- -2*(logLik(fm.reduced.new)-logLik(fm.full.new))
}


## ----fig.cap="Hisogram of simulated test statistics. The green curve represents a chisq distribution with 1 df. The red vertical line are the observed test statistic."----

hist(boot.test.stats, prob=T, breaks=40, xlim = c(0, obs.test.stat+3))
curve(dchisq(x,df=1),from=0,to=40,add=T,col=3)
points(obs.test.stat,0,pch=4,col="red")
abline(v=obs.test.stat, col=2, lwd=2)
mean(boot.test.stats > obs.test.stat)


