Lanna=read.csv2("../data/R1001.csv")
head(Lanna)
LTime=read.csv2("../data/Time1001.csv")
head(LTime)

HighControl <- subset(LTime, Treatment=="Control" & Field=="High")
HighP <- subset(LTime, Treatment=="P" & Field=="High")
HighP6 <- subset(LTime, Treatment=="P6" & Field=="High")

LowControl <- subset(LTime, Treatment=="Control" & Field=="Low")
LowP <- subset(LTime, Treatment=="P" & Field=="Low")
LowP6 <- subset(LTime, Treatment=="P6" & Field=="Low")



##Fig 1 left
pdf(file="../Figures/Fig1LannaLeft.pdf", encoding='WinAnsi.enc', width=4, height=7)
par(mfrow=c(3,1), mar=c(3,4.6,2,1))
#TOC
plot(LTime$TOC~LTime$Year,col="black",ylim=c(10,26),xlab="",ylab=expression(paste("TOC [g  ", kg^-1,"]")))
modTOC=lm(LTime$TOC~LTime$Year)
abline(modTOC)
summary(modTOC)
sqrt(mean(modTOC$residuals^2))
AIC(modTOC)

modC=nls(TOC~b*exp(a*(Year-1968)), data=LTime, start=c(a=-0.1, b=21))
x=0-85
y1=coef(modC)[2]*exp(coef(modC)[1]*x)
summary(modC)
AIC(modC)

points(HighControl$Year,HighControl$TOC,pch=1, col="black")
points(HighP$Year,HighP$TOC,pch=4, col="blue")
points(HighP6$Year,HighP6$TOC,pch=8, col="red")
points(LowControl$Year,LowControl$TOC,pch=1, col="black")
points(LowP$Year,LowP$TOC,pch=4, col="blue")
points(LowP6$Year,LowP6$TOC,pch=8, col="red")
legend("bottomleft",expression(paste("y = 125.00-0.05308x, R"^2," = 0.25, P < 0.001")), bty="n")
legend("topleft", "a", bty="n", cex=1.2,text.font=2)

#TN
plot(LTime$TN~LTime$Year,col="white", ylim=c(1.2,2.2),pch=1,xlab="",ylab=expression(paste("TN [g  ", kg^-1,"]")))
modTN=lm(LTime$TN~LTime$Year)
abline(modTN)
summary(modTN)
sqrt(mean(modTN$residuals^2))
AIC(modTN)

modN=nls(TN~b*exp(a*(Year-1968)), data=LTime, start=c(a=-0.1, b=21))
x=0-85
y1=coef(modN)[2]*exp(coef(modN)[1]*x)
summary(modN)
AIC(modN)

points(HighControl$Year,HighControl$TN,pch=1, col="black")
points(HighP$Year,HighP$TN,pch=4, col="blue")
points(HighP6$Year,HighP6$TN,pch=8, col="red")
points(LowControl$Year,LowControl$TN,pch=1, col="black")
points(LowP$Year,LowP$TN,pch=4, col="blue")
points(LowP6$Year,LowP6$TN,pch=8, col="red")
legend("bottomleft",expression(paste("y = 11.06-0.004701x, R"^2," = 0.31, P < 0.001")), bty="n")
legend(x=2005, y=2.16, legend=c("Control","P","P6"), col=c("black","blue","red"), pch=c(1,4,8))
legend("topleft", "c", bty="n", cex=1.2,text.font=2)

##TOP
plot(LTime$TOP~LTime$Year,col="white", ylim=c(200,420),pch=1,xlab="",ylab=expression(paste("TOP [mg  ", kg^-1,"]")))
modTOP=lm(LTime$TOP~LTime$Year)
abline(modTOP)
summary(modTOP)
sqrt(mean(modTOP$residuals^2))
AIC(modTOP)

modOP=nls(TOP~b*exp(a*(Year-1968)), data=LTime, start=c(a=-0.1, b=21))
x=0-85
y1=coef(modOP)[2]*exp(coef(modOP)[1]*x)
summary(modOP)
AIC(modOP)

points(HighControl$Year,HighControl$TOP,pch=1, col="black")
points(HighP$Year,HighP$TOP,pch=4, col="blue")
points(HighP6$Year,HighP6$TOP,pch=8, col="red")
points(LowControl$Year,LowControl$TOP,pch=1, col="black")
points(LowP$Year,LowP$TOP,pch=4, col="blue")
points(LowP6$Year,LowP6$TOP,pch=8, col="red")
legend("bottomleft",expression(paste("y = 1742.16-0.719x, R"^2," = 0 .11, P < 0.05")), bty="n")
legend("topleft", "e", bty="n", cex=1.2,text.font=2)
dev.off()



##Fig 4 left
pdf(file="../Figures/Fig4LannaLeft.pdf", encoding='WinAnsi.enc', width=4, height=7)
par(mfrow=c(3,1), mar=c(3,4.6,2,1))
##TP
plot(LTime$TP~LTime$Year,col="white", ylim=c(400, 700),pch=1,xlab="",ylab=expression(paste("TP [mg  ", kg^-1,"]")))
points(HighControl$Year,HighControl$TP,pch=1, col="black")
points(HighP$Year,HighP$TP,pch=4, col="blue")
points(HighP6$Year,HighP6$TP,pch=8, col="red")
points(LowControl$Year,LowControl$TP,pch=1, col="black")
points(LowP$Year,LowP$TP,pch=4, col="blue")
points(LowP6$Year,LowP6$TP,pch=8, col="red")
modTP=lm(LTime$TP~LTime$Year)
abline(modTP)
summary(modTP)
AIC(modTP)
modP=nls(TP~b*exp(a*(Year-1968)), data=LTime, start=c(a=0.1, b=21))
x=0-85
y1=coef(modP)[2]*exp(coef(modP)[1]*x)
summary(modP)
AIC(modP)

sqrt(mean(modTP$residuals^2))
legend("bottomleft",expression(paste("y = 3044.00-1.2415x, R"^2," = 0.13, P < 0.05")), bty="n")
legend("topleft", "a", bty="n", cex=1.1,text.font=2)

##Phytate
plot(HighControl$Phytate~HighControl$Year,pch=1, col="black", ylim=c(45,150),xlab="",ylab=expression(paste("Phytate-P [mg  ", kg^-1,"]")))
points(HighP$Year,HighP$Phytate,pch=4, col="blue")
points(HighP6$Year,HighP6$Phytate,pch=8, col="red")
points(LowControl$Year,LowControl$Phytate,pch=1, col="black")
points(LowP$Year,LowP$Phytate,pch=4, col="blue")
points(LowP6$Year,LowP6$Phytate,pch=8, col="red")
modPhy=lm(LTime$Phytate~LTime$Year)
abline(modPhy)
summary(modPhy)
AIC(modPhy)

modPhytate=nls(Phytate~b*exp(a*(Year-1968)), data=LTime, start=c(a=0.1, b=21))
x=0-85
y1=coef(modPhytate)[2]*exp(coef(modPhytate)[1]*x)
summary(modPhytate)
AIC(modPhytate)

sqrt(mean(modPhy$residuals^2))
legend("bottomleft",expression(paste("y = 629.37-0.2676x, R"^2," = 0.07, P = 0.052")), bty="n")
legend("topleft", "c", bty="n", cex=1.1,text.font=2)


##PAL
plot(LTime$PAL~LTime$Year,col="white", ylim=c(5, 75),xlab="",ylab=expression(paste("P-AL [mg  ", kg^-1,"]")))
points(HighControl$Year,HighControl$PAL,pch=1, col="black")
points(HighP$Year,HighP$PAL,pch=4, col="blue")
points(HighP6$Year,HighP6$PAL,pch=8, col="red")
points(LowControl$Year,LowControl$PAL,pch=1, col="black")
points(LowP$Year,LowP$PAL,pch=4, col="blue")
points(LowP6$Year,LowP6$PAL,pch=8, col="red")
modPAL=lm(LTime$PAL~LTime$Year)
#abline(modPAL)
summary(modPAL)
legend(x=1990, y=75, legend=c("Control","P","P6"), col=c("black","blue","red"), pch=c(1,4,8))
legend("topleft", "e", bty="n", cex=1.1,text.font=2)
par(mfrow=c(1,1), mar=c(4,4.4,2,1))
dev.off()




#Fig. 1 right
pdf(file="../Figures/Fig1LannaRight.pdf", encoding='WinAnsi.enc', width=2, height=7)
par(mfrow=c(3,1), mar=c(3,4.6,2,1))
##TOC
MeanTOCTr=aggregate(Lanna$TOC*260, by=list(Lanna$Treatment), mean, na.rm=TRUE)
SDTOCTr=aggregate(Lanna$TOC*260, by=list(Lanna$Treatment), sd, na.rm=TRUE)

B<-barplot(MeanTOCTr[,2], axes = FALSE,cex.lab=1.1, ylim=c(0,7000), names.arg=c("Ctrl","P","P6"), col=c("black", "blue", "red"),cex.names=1.1,ylab=expression(paste("TOC [g  ", m^-2,"]")))
box()
axis(2)
arrows(B, MeanTOCTr[,2]-SDTOCTr[,2], B, MeanTOCTr[,2]+SDTOCTr[,2],angle=90,code=3,length=0.05)
legend("topleft", "b", bty="n", cex=1.1,text.font=2)
legend("topright", "P > 0.05", bty="n", cex=1)
text(y=2000, x=c(0.7,1.9,3.1), round(MeanTOCTr[,2],digits=0), col="white", pch=15)

##TN
MeanTNTr=aggregate(Lanna$TN*260, by=list(Lanna$Treatment), mean, na.rm=TRUE)
SDTNTr=aggregate(Lanna$TN*260, by=list(Lanna$Treatment), sd, na.rm=TRUE)

B<-barplot(MeanTNTr[,2], axes = FALSE,cex.lab=1.1, ylim=c(0,600), names.arg=c("Ctrl","P","P6"), col=c("black", "blue", "red"),cex.names=1.1,ylab=expression(paste("TN [g  ", m^-2,"]")))
box()
axis(2)
arrows(B, MeanTNTr[,2]-SDTNTr[,2], B, MeanTNTr[,2]+SDTNTr[,2],angle=90,code=3,length=0.05)
legend("topleft", "d", bty="n", cex=1.1,text.font=2)
legend("topright", "P > 0.05", bty="n", cex=1)
text(y=50, x=c(0.7,1.9,3.1), round(MeanTNTr[,2],digits=0), col="white", pch=15)

##TOP
MeanTOPTr=aggregate(Lanna$TOP*260/1000, by=list(Lanna$Treatment), mean, na.rm=TRUE)
SDTOPTr=aggregate(Lanna$TOP*260/1000, by=list(Lanna$Treatment), sd, na.rm=TRUE)

B<-barplot(MeanTOPTr[,2], axes = FALSE,cex.lab=1.1, ylim=c(0,110), names.arg=c("Ctrl","P","P6"), col=c("black", "blue", "red"),cex.names=1.1,ylab=expression(paste("TOP [g  ", m^-2,"]")))
box()
axis(2)
arrows(B, MeanTOPTr[,2]-SDTOPTr[,2], B, MeanTOPTr[,2]+SDTOPTr[,2],angle=90,code=3,length=0.05)
legend("topleft", "f", bty="n", cex=1.1,text.font=2)
legend("topright", "P > 0.05", bty="n", cex=1)
text(y=20, x=c(0.7,1.9,3.1), round(MeanTOPTr[,2],digits=2), col="white", pch=15)
dev.off()



##Fig. 4 right
pdf(file="../Figures/Fig4LannaRight.pdf", encoding='WinAnsi.enc', width=2, height=7)
par(mfrow=c(3,1), mar=c(3,4.6,2,1))
##TP
MeanTPTr=aggregate(Lanna$TP*260/1000, by=list(Lanna$Treatment), mean, na.rm=TRUE)
SDTPTr=aggregate(Lanna$TP*260/1000, by=list(Lanna$Treatment), sd, na.rm=TRUE)

B<-barplot(MeanTPTr[,2], axes = FALSE,cex.lab=1.1, ylim=c(0,250), names.arg=c("Ctrl","P","P6"), col=c("black", "blue", "red"),cex.names=1.1,ylab=expression(paste("TP [g  ", m^-2,"]")))
box()
axis(2)
arrows(B, MeanTPTr[,2]-SDTPTr[,2], B, MeanTPTr[,2]+SDTPTr[,2],angle=90,code=3,length=0.05)
legend("topleft", "b", bty="n", cex=1.1,text.font=2)
legend("topright", "P > 0.05", bty="n", cex=1)
text(y=20, x=c(0.7,1.9,3.1), round(MeanTPTr[,2],digits=1), col="white", pch=15)

group.factor2=as.factor(Lanna$Treatment)
aov.C3 = aov((log(Lanna$TP))~group.factor2)
TukeyC3 <- TukeyHSD(aov.C3, conf.level=.95) 
summary(aov.C3)

##Phytate
MeanPhytateTr=aggregate(Lanna$Phytate*260/1000, by=list(Lanna$Treatment), mean, na.rm=TRUE)
SDPhytateTr=aggregate(Lanna$Phytate*260/1000, by=list(Lanna$Treatment), sd, na.rm=TRUE)

B<-barplot(MeanPhytateTr[,2], axes = FALSE,cex.lab=1.1, ylim=c(0,40), names.arg=c("Ctrl","P","P6"), col=c("black", "blue", "red"),cex.names=1.1,ylab=expression(paste("Phytate-P [g  ", m^-2,"]")))
box()
axis(2)
arrows(B, MeanPhytateTr[,2]-SDPhytateTr[,2], B, MeanPhytateTr[,2]+SDPhytateTr[,2],angle=90,code=3,length=0.05)
legend("topleft", "d", bty="n", cex=1.1,text.font=2)
legend("topright", "P > 0.05", bty="n", cex=1)
text(y=5, x=c(0.7,1.9,3.1), round(MeanPhytateTr[,2],digits=2), col="white", pch=15)

group.factor2=as.factor(Lanna$Treatment)
aov.C2 = aov((log(Lanna$Phytate))~group.factor2)
TukeyC2 <- TukeyHSD(aov.C2, conf.level=.95) 
summary(aov.C2)

##PAL
MeanPALTr=aggregate(Lanna$PAL*260/1000, by=list(Lanna$Treatment), mean, na.rm=TRUE)
SDPALTr=aggregate(Lanna$PAL*260/1000, by=list(Lanna$Treatment), sd, na.rm=TRUE)

B<-barplot(MeanPALTr[,2], axes = FALSE,cex.lab=1.1, ylim=c(0,14), names.arg=c("Ctrl","P","P6"), col=c("black", "blue", "red"),cex.names=1.1,ylab=expression(paste("P-AL [g  ", m^-2,"]")))
box()
box()
axis(2)
arrows(B, MeanPALTr[,2]-SDPALTr[,2], B, MeanPALTr[,2]+SDPALTr[,2],angle=90,code=3,length=0.05)
legend("topleft", "f", bty="n", cex=1.1,text.font=2)
legend("topright", "P < 0.01", bty="n", cex=1)
text(x=c(0.7,1.9,3.1), y=11, c("B","A","A"), col=1, pch=15)
text(y=2, x=c(0.7,1.9,3.1), round(MeanPALTr[,2],digits=2), col="white", pch=15)

group.factor2=as.factor(Lanna$Treatment)
aov.C1 = aov((log(Lanna$PAL))~group.factor2)
TukeyC1 <- TukeyHSD(aov.C1, conf.level=.95) 
summary(aov.C1)
dev.off()


##Fig. S1 Pgrain 
Pgrain1=read.csv2("../data/Pgrain1001.csv")
head(Pgrain1)

HighControl <- subset(Pgrain1, Treatment=="Control" & Field=="High")
HighP <- subset(Pgrain1, Treatment=="P" & Field=="High")
HighP6 <- subset(Pgrain1, Treatment=="P6" & Field=="High")

LowControl <- subset(Pgrain1, Treatment=="Control" & Field=="Low")
LowP <- subset(Pgrain1, Treatment=="P" & Field=="Low")
LowP6 <- subset(Pgrain1, Treatment=="P6" & Field=="Low")

plot(HighControl$Pgrain~HighControl$Year,pch=1, col="black",ylim=c(0,7.5),xlab="Calendar year",ylab=expression(paste("P in grains [mg  ", kg^-1,"]")))
points(HighP$Year,HighP$Pgrain,pch=4, col="blue")
points(HighP6$Year,HighP6$Pgrain,pch=8, col="red")
points(LowControl$Year,LowControl$Pgrain,pch=1, col="black")
points(LowP$Year,LowP$Pgrain,pch=4, col="blue")
points(LowP6$Year,LowP6$Pgrain,pch=8, col="red")
legend(x=2000, y=1.9, legend=c("Control","P","P6"), col=c("black","blue","red"), pch=c(1,4,8))


##Fig. S2 Soil pH
plot(LTime$pH~LTime$Year,col="white",ylim=c(2,12),xlab="Calendar year",ylab="Soil pH")

points(HighControl$Year,HighControl$pH,pch=1, col="black")
points(HighP$Year,HighP$pH,pch=4, col="blue")
points(HighP6$Year,HighP6$pH,pch=8, col="red")
points(LowControl$Year,LowControl$pH,pch=1, col="black")
points(LowP$Year,LowP$pH,pch=4, col="blue")
points(LowP6$Year,LowP6$pH,pch=8, col="red")

legend(x=2005, y=4, legend=c("Control","P","P6"), col=c("black","blue","red"), pch=c(1,4,8))

