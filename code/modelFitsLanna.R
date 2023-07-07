# Required packages
library(SoilR)
library(FME)

# Working directory should be set to the current location of this file
# Read and organize data in path relative to current location of this file
dt=read.csv2("../data/Time1001.csv")
yd=read.csv("../data/Yield1001.csv")

additionalData<-c(-21.6, -17.7, -19.6, -31.1, -22.0, -31.0, -18.8, -26.5, -13.9, -18.0, -29.2, -24.1) # Data from AMS batch Sierra_37_LAN
dt$X14C[31:42]<-additionalData

BD=1.3 # g/cm3
depth=20 #  cm

Cobs=data.frame(Year=dt[,1],  Ct=dt[,9]*BD*depth*10) # Multiply by 10 to change units to g/m2
C14obs=data.frame(Year=dt[,1],C14t=dt[,7])

avgYld=data.frame(Year=yd[,1],
                  Yield=apply(yd[,-1],MARGIN = 1, FUN = mean, na.rm=TRUE)*0.5*0.1) # Multiply by 0.5 to change from DM to C, and by 0.1 from kg/ha to g/m2

# Constant inputs from Bolinder et al. 2012, doi: 10.4141/cjss2012-036 Table 4
BolinderIn<-c(0.277,0.268,0.304,0.271,0.269,0.29, # Rotation A
              0.259, 0.291, 0.296, 0.288, 0.265, 0.256) # Rotation B

C0<-mean(Cobs[1:2,2])
F0<-mean(C14obs[1:2,2])

# Model
yr=seq(1968,2021,by=1/12)
Atm14C=Hua2021$NHZone1[,1:2]
fAtm14C=read.csv("../data/NHZ2forecast.csv")
Atm14C<-rbind(Atm14C, data.frame(Year=fAtm14C$time, mean.Delta14C=Delta14C_from_AbsoluteFractionModern(fAtm14C$F14C)))

mf=function(pars){
  md=TwopSeriesModel14(t=yr,ks=pars[1:2],C0=C0*c(pars[5], 1-pars[5]),
                       F0_Delta14C = c(F0 * pars[6], -200), 
                       In=data.frame(avgYld[,1],avgYld[,2]*pars[4]), 
                       a21=pars[1]*pars[3], inputFc = Atm14C)
  Ct=getC(md)
  C14t=getF14C(md)
  
  return(data.frame(Year=yr,Ct=rowSums(Ct),C14t=C14t))
}

mc=function(pars){
  out=mf(pars)
  Cost1=modCost(model=out,obs=C14obs,x="Year")
  return(modCost(model=out,obs=Cobs, x="Year",cost=Cost1))
}

inipars=c(0.5,0.001,0.01, 0.5, 0.1, 0.9)
# Model fit results are already saved. Uncomment the following code to run again or simply load previous results
# mFit=modFit(f=mc,p=inipars,method="Nelder-Mead",upper=c(1,0.5,1,2,1,1),lower=c(0,0,0,0,0,0))
# bestpars=mFit$par
# save(bestpars, file="bestpars.RData")
load("../data/bestpars.RData")
bestModel=mf(bestpars)

mod1=TwopSeriesModel14(t=yr,ks=bestpars[1:2],C0=C0*c(bestpars[5], 1-bestpars[5]),
                       F0_Delta14C = c(F0 * bestpars[6], -200), 
                       In=data.frame(avgYld[,1],avgYld[,2]*bestpars[4]), 
                       a21=bestpars[1]*bestpars[3], inputFc = Atm14C)
F14C1=getF14(mod1)
C1=getC(mod1)

# var0 <- mFit$var_ms_unweighted
# cov0 <- summary(mFit)$cov.scaled # The covariance matrix can be used for the jump, but wasn't used in this example.
# MCMC <- modMCMC(f=mc, p = bestpars, niter = 2500, jump = NULL, var0 = var0, wvar0 = 1)
# save(MCMC, file="MCMC.RData")
load("../data/MCMC.RData")

parsMCMC<-summary(MCMC)
# sR=sensRange(func=mf, parInput=MCMC$par)
# save(sR, file="sR.RData")
load("../data/sR.RData")
summarySR<-summary(sR)
nx<-attributes(summarySR)$nx

CtR<-as.data.frame(summarySR[1:nx,])
C14tR<-as.data.frame(summarySR[(nx+1):(2*nx),])

plot(CtR[,1:2], type="l", xlab="Year", ylab="C stock", ylim=c(0,6000)) 
polygon(c(CtR$x,rev(CtR$x)), c(CtR$Min, rev(CtR$Max)) ,col=gray(0.8), border="NA")
polygon(c(CtR$x,rev(CtR$x)), c(CtR$Mean+CtR$Sd, rev(CtR$Mean-CtR$Sd)) ,col=gray(0.5), border="NA")
lines(CtR[,1:2])
points(Cobs,pch=20)

plot(C14tR[,1:2], type="l", xlim=c(1950,2020), ylim=c(-120,400), ylab=expression(paste(Delta^14,"C ","(\u2030)")),xlab="Years")
polygon(c(C14tR$x,rev(C14tR$x)), c(C14tR$Min, rev(C14tR$Max)) ,col=gray(0.8), border="NA")
polygon(c(C14tR$x,rev(C14tR$x)), c(C14tR$Mean+C14tR$Sd, rev(C14tR$Mean-C14tR$Sd)) ,col=gray(0.5), border="NA")
lines(C14tR[,1:2])
points(C14obs, pch=20)
lines(Atm14C,col=4)

pairs(MCMC,nsample=500)

meanpars<-as.numeric(parsMCMC[1,1:5])
meanModel<-TwopSeriesModel14(t=yr,ks=meanpars[1:2],C0=C0*c(meanpars[5], 1-meanpars[5]),
                             F0_Delta14C =c(F0* meanpars[6], -200), 
                             In=data.frame(avgYld[,1],avgYld[,2]*meanpars[4]), 
                             a21=meanpars[1]*meanpars[3], inputFc = Atm14C)
meanF14C=getF14(meanModel)
meanC=getC(meanModel)

# Constant inputs
mf2=function(pars){
  md=TwopSeriesModel14(t=yr,ks=pars[1:2],C0=C0*c(pars[4], 1-pars[4]),
                       F0_Delta14C = c(F0 * pars[5], -200), 
                       #In=mean(BolinderIn)*1000, 
                       In=mean(avgYld[,2], na.rm = TRUE)*1.5,
                       a21=pars[1]*pars[3], inputFc = Atm14C)
  Ct=getC(md)
  C14t=getF14C(md)
  
  return(data.frame(Year=yr,Ct=rowSums(Ct),C14t=C14t))
}

mc2=function(pars){
  out=mf2(pars)
  Cost1=modCost(model=out,obs=C14obs,x="Year")
  return(modCost(model=out,obs=Cobs, x="Year",cost=Cost1))
}

# inipars2=c(0.5,0.001,0.01, 0.1, 0.9)
# mFit2=modFit(f=mc2,p=inipars2,method="Nelder-Mead",upper=c(1,0.5,1,1,1),lower=c(0,0,0,0,0))
# 
# bestpars2=mFit2$par
# save(bestpars2, file="bestpars2.RData")
load("../data/bestpars2.RData")
bestModel2=mf2(bestpars2)

mod2=TwopSeriesModel14(t=yr,ks=bestpars2[1:2],C0=C0*c(bestpars2[4], 1-bestpars2[4]),
                       F0_Delta14C = c(F0 * bestpars2[5], -200), 
                       #In=mean(BolinderIn)*1000, 
                       In=mean(avgYld[,2], na.rm = TRUE)*1.5,
                       a21=bestpars2[1]*bestpars2[3], inputFc = Atm14C)
F14C2=getF14(mod2)
C2=getC(mod2)

# var02 <- mFit2$var_ms_unweighted
# cov02 <- summary(mFit2)$cov.scaled # The covariance matrix can be used for the jump, but wasn't used in this example.
# MCMC2 <- modMCMC(f=mc2, p = bestpars2, niter = 2500, jump = NULL, var0 = var02, wvar0 = 1)
# save(MCMC2, file="MCMC2.RData")
load("../data/MCMC2.RData")

parsMCMC2<-summary(MCMC2)
# sR2=sensRange(func=mf2, parInput=MCMC2$par)
# save(sR2, file="sR2.RData")
load("../data/sR2.RData")
summarySR2<-summary(sR2)
nx2<-attributes(summarySR2)$nx

CtR2<-as.data.frame(summarySR2[1:nx2,])
C14tR2<-as.data.frame(summarySR2[(nx2+1):(2*nx2),])

plot(CtR2[,1:2], type="l", xlab="Year", ylab="C stock", ylim=c(0,6000)) 
polygon(c(CtR2$x,rev(CtR2$x)), c(CtR2$Min, rev(CtR2$Max)) ,col=gray(0.8), border="NA")
polygon(c(CtR2$x,rev(CtR2$x)), c(CtR2$Mean+CtR2$Sd, rev(CtR2$Mean-CtR2$Sd)) ,col=gray(0.5), border="NA")
lines(CtR2[,1:2])
points(Cobs,pch=20)

plot(C14tR2[,1:2], type="l", xlim=c(1950,2020), ylim=c(-120,400), ylab=expression(paste(Delta^14,"C ","(\u2030)")),xlab="Years")
polygon(c(C14tR2$x,rev(C14tR2$x)), c(C14tR2$Min, rev(C14tR2$Max)) ,col=gray(0.8), border="NA")
polygon(c(C14tR2$x,rev(C14tR2$x)), c(C14tR2$Mean+C14tR2$Sd, rev(C14tR2$Mean-C14tR2$Sd)) ,col=gray(0.5), border="NA")
lines(C14tR2[,1:2])
points(C14obs, pch=20)
lines(Atm14C,col=4)

meanpars2<-as.numeric(parsMCMC2[1,1:5])
meanModel2<-TwopSeriesModel14(t=yr,ks=meanpars2[1:2],C0=C0*c(meanpars[4], 1-meanpars[4]),
                             F0_Delta14C = c(F0*meanpars[5], -200), 
                             #In=mean(BolinderIn)*1000, 
                             In=mean(avgYld[,2], na.rm = TRUE)*1.5,
                             a21=meanpars[1]*meanpars[3], inputFc = Atm14C)
meanF14C2=getF14(meanModel2)
meanC2=getC(meanModel2)



# Age and transit time
tau=seq(0,500)
A1=meanModel@mat@map(1970)
A2=meanModel2@mat@map(1970)

A1min<- -1*diag(parsMCMC[3,1:2])
A1min[2,1]<- abs(parsMCMC[3,3])*parsMCMC[3,2] # alpha21 is virtually zero
A1max<- -1*diag(parsMCMC[4,1:2])
A1max[2,1]<- parsMCMC[4,3]*parsMCMC[4,2] 
A2min<- -1*diag(parsMCMC2[3,1:2])
A2min[2,1]<- parsMCMC2[3,3]*parsMCMC2[3,2] 
A2max<- -1*diag(parsMCMC2[4,1:2])
A2max[2,1]<- parsMCMC2[4,3]*parsMCMC2[4,2] 


SA1=systemAge(A=A1,u=c(1,0),a=tau)
TT1=transitTime(A=A1,u=c(1,0), a=tau, q=c(0.1,0.5,0.9))
TT1min<-transitTime(A=A1min,u=c(1,0),q=c(0.1,0.5,0.9))
TT1max<-transitTime(A=A1max,u=c(1,0),q=c(0.1,0.5,0.9))

SA2=systemAge(A=A2,u=c(1,0),a=tau)
TT2=transitTime(A=A2,u=c(1,0), a=tau, q=c(0.1,0.5,0.9))
TT2min<-transitTime(A=A2min,u=c(1,0),q=c(0.1,0.5,0.9))
TT2max<-transitTime(A=A2max,u=c(1,0),q=c(0.1,0.5,0.9))

In=data.frame(avgYld[,1],avgYld[,2]*bestpars[4])

pdf("../Figures/Fig2.pdf", encoding = 'WinAnsi.enc', width=7.5*sqrt(2), height=7.5)
par(mfrow=c(2,2), mar=c(4,4.5,1,1), cex.lab=1.2)
#plot(avgYld$Year, avgYld$Yield*2,type="o", xlab="Calendar year", ylab=expression(paste("Mean yields [g DW ", m^-2, " y", r^-1, "]" )), xlim=c(1937,2022))
plot(In, type="b", ylim=c(0,700), xlab="Calendar year", ylab=expression(paste("Carbon inputs [g C ", m^-2, " y", r^-1,"]")), col=rgb(0,0,1), pch=20)
#abline(h=mean(BolinderIn)*1000, col=2)
abline(h=mean(avgYld[,2], na.rm = TRUE)*1.5, col=rgb(0,1,0))
legend("bottomright", c("Variable inputs predicted by the model", "Constant inputs from Bolinder et al. (2012)"), lty=1, 
       col=c(rgb(0,0,1), rgb(0,1,0)), bty="n")
legend("topleft", "a", cex=1.2, text.font=2, bty="n")

par(mar=c(4,4.5,1,1))
plot(C14tR[,1:2],type="l",ylim=c(-200,400), xlab="Calendar year", ylab=expression(paste(Delta^14,"C [\U2030]")), xlim=c(1950,2022))
polygon(c(C14tR$x,rev(C14tR$x)), c(C14tR$Min, rev(C14tR$Max)) ,col=rgb(0,0,1, alpha=0.3), border="NA")
polygon(c(C14tR$x,rev(C14tR$x)), c(C14tR$Mean+C14tR$Sd, rev(C14tR$Mean-C14tR$Sd)) ,col=rgb(0,0,1,alpha=0.5), border="NA")
lines(C14tR[,1:2], col=rgb(0,0,1))
polygon(c(C14tR2$x,rev(C14tR2$x)), c(C14tR2$Min, rev(C14tR2$Max)) ,col=rgb(0,1,0, alpha=0.3), border="NA")
polygon(c(C14tR2$x,rev(C14tR2$x)), c(C14tR2$Mean+C14tR2$Sd, rev(C14tR2$Mean-C14tR2$Sd)) ,col=rgb(0,1,0, alpha=0.5), border="NA")
lines(C14tR2[,1:2], col=rgb(0,1,0))
points(C14obs, pch=20)
lines(Atm14C)
legend("topright", c("Measurements", "Model with constant inputs", "Model with variable inputs", "Atmosphere"), 
       lty=c(NA,1,1,1), pch=c(19, NA, NA, NA), col=c(1,rgb(0,1,0), rgb(0,0,1),1), bty="n")
legend("topleft", "b", cex=1.2, text.font=2, bty="n")

matplot(yr,C1,type="l",col=rgb(0,0,1),lty=2:3, lwd=2,ylim=c(0,10000), xlim=c(1966,2022), xlab="Calendar year",
        ylab=expression(paste("Carbon stock [g C ", m^-2, "]")))
polygon(c(CtR$x,rev(CtR$x)), c(CtR$Min, rev(CtR$Max)) ,col=rgb(0,0,1, alpha=0.3), border="NA")
polygon(c(CtR$x,rev(CtR$x)), c(CtR$Mean+CtR$Sd, rev(CtR$Mean-CtR$Sd)) ,col=rgb(0,0,1, alpha=0.5), border="NA")
lines(CtR[,1:2], col=rgb(0,0,1))
matlines(yr, C2, col=rgb(0,1,0), lty=2:3, lwd=2)
polygon(c(CtR2$x,rev(CtR2$x)), c(CtR2$Min, rev(CtR2$Max)) ,col=rgb(0,1,0, alpha=0.3), border="NA")
polygon(c(CtR2$x,rev(CtR2$x)), c(CtR2$Mean+CtR2$Sd, rev(CtR2$Mean-CtR2$Sd)) ,col=rgb(0,1,0, alpha=0.5), border="NA")
lines(CtR2[,1:2], col=rgb(0,1,0))
points(Cobs, pch=20)
legend(x=1975, y=10500, c("Measurements", "TOC, constant inputs", "Fast pool", "Slow pool"), lty=c(NA,1,2,3), 
       pch=c(19, NA, NA, NA), col=c(1,rep(rgb(0,1,0),3)), bty="n")
legend("topright", c("TOC, variable inputs", "Fast pool", "Slow pool"), lty=c(1,2,3), 
       col=c(rep(rgb(0,0,1),3)), bty="n")
legend("topleft", "c", cex=1.2, text.font=2, bty="n")

plot(tau,log(TT1$transitTimeDensity), type="l", xlim=c(0,200), xlab="Transit time (yr)", ylab="Log probability density", col=rgb(0,0,1))
abline(v=TT1$meanTransitTime, lty=2, col=rgb(0,0,1))
#abline(v=TT1$quantiles[2],lty=3, col=4)
lines(tau, log(TT2$transitTimeDensity), col=rgb(0,1,0))
abline(v=TT2$meanTransitTime, lty=2, col=rgb(0,1,0))
#abline(v=TT2$quantiles[2],lty=3, col=2)
legend(x=60,y=-4,c("Transit time distribution, variable inputs", paste("Mean transit time = ", round(TT1$meanTransitTime,1), " yr"))
       ,lty=1:3, col=rgb(0,0,1), bty="n")
legend(x=60,y=-6,c("Transit time distribution, constant inputs", paste("Mean transit time = ", round(TT2$meanTransitTime,1), " yr"))
       ,lty=1:3, col=rgb(0,1,0), bty="n")
legend(x=5,y=-1, "d", cex=1.2, text.font=2, bty="n")
par(mfrow=c(1,1))
dev.off()


data.frame(Parameter=c("kf", "ks", "alpha_sf", "gamma", "beta"), 
           ModelConst=paste(c(round(bestpars2[1:3], 3), NA, round(bestpars2[4], 3)), "+-", c(round(parsMCMC2[2,1:3],3), NA, round(parsMCMC2[2,4], 3))),
           ModelVar=paste(round(bestpars[1:5], 3), "+-", round(parsMCMC[2,1:5], 3))
           )

paste("Mean TT (uncertainty): ",round(TT1$meanTransitTime, 2), "(", round(TT1max$meanTransitTime, 2), "--", round(TT1min$meanTransitTime, 2),")")
paste("Mean TT (uncertainty): ",round(TT2$meanTransitTime, 2), "(", round(TT2max$meanTransitTime, 2), "--", round(TT2min$meanTransitTime, 2),")")

paste("Median TT (uncertainty): ",round(TT1$quantiles[2], 2), "(", round(TT1max$quantiles[2], 2), "--", round(TT1min$quantiles[2], 2),")")
paste("Median TT (uncertainty): ",round(TT2$quantiles[2], 2), "(", round(TT2max$quantiles[2], 2), "--", round(TT2min$quantiles[2], 2),")")

round(TT1$quantiles, 2)
round(TT1min$meanTransitTime, 2)
round(TT1min$quantiles, 2)
round(TT1max$meanTransitTime, 2)
round(TT1max$quantiles, 2)

round(TT2$meanTransitTime, 2)
round(TT2$quantiles, 2)

################################################################################
# decrease in the size of the slow C pool by about x% from 1968 to 2021
start=meanC[yr==1968,]
end=meanC[yr==2021,]
diff=end-start
100*diff/start

################################################################################
# Simulation without using radiocarbon as constraint. Only to answer reviewers' comments

mc3=function(pars){
  out=mf2(pars)
  return(modCost(model=out,obs=Cobs, x="Year"))
}

MCMC3 <- modMCMC(f=mc3, p = bestpars2, niter = 2500)

parsMCMC3<-summary(MCMC3)
sR3=sensRange(func=mf2, parInput=MCMC3$par)
#save(../data/sR3, file="sR2.RData")
#load("../data/sR3.RData")
summarySR3<-summary(sR3)
nx3<-attributes(summarySR3)$nx

CtR3<-as.data.frame(summarySR3[1:nx3,])
C14tR3<-as.data.frame(summarySR3[(nx3+1):(2*nx3),])

pdf("../Figures/No14C.pdf", encoding = 'WinAnsi.enc')
par(mfrow=c(2,1))
plot(CtR3[,1:2], type="l", xlab="Year", ylab="C stock", ylim=c(0,6000)) 
polygon(c(CtR3$x,rev(CtR3$x)), c(CtR3$Min, rev(CtR3$Max)) ,col=gray(0.8), border="NA")
polygon(c(CtR3$x,rev(CtR3$x)), c(CtR3$Mean+CtR3$Sd, rev(CtR3$Mean-CtR3$Sd)) ,col=gray(0.5), border="NA")
lines(CtR3[,1:2])
points(Cobs,pch=20)

plot(C14tR3[,1:2], type="l", xlim=c(1950,2020), ylim=c(-120,400), ylab=expression(paste(Delta^14,"C ","(\u2030)")),xlab="Years")
polygon(c(C14tR3$x,rev(C14tR3$x)), c(C14tR3$Min, rev(C14tR3$Max)) ,col=gray(0.8), border="NA")
polygon(c(C14tR3$x,rev(C14tR3$x)), c(C14tR3$Mean+C14tR3$Sd, rev(C14tR3$Mean-C14tR3$Sd)) ,col=gray(0.5), border="NA")
lines(C14tR3[,1:2])
points(C14obs, pch=20)
lines(Atm14C,col=4)
par(mfrow=c(1,1))
dev.off()
