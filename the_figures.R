load("LocustAnalysis.RData")

load("PKref.RData")
dN2 <- paste("PK_",1:4,".fl.gz",sep="")
nb <- 58*15000
pkD <- sapply(dN2,
              function(n) {
                  mC <- gzfile(n,open="rb")
                  x <- readBin(mC,what="double",size=4,n=nb)
                  close(mC);x})
colnames(pkD) <- paste("site",1:4)

png(file="Pouzat_Fig_1.png",width=1000,height=600)
plot(c(0,15000),c(-105,15),type="n",axes = FALSE,
     xlab="",ylab="")
lines(lD[1:15000,1],lwd=0.5)
text(250,10,"Site 1",cex=2)
text(2712,9,"S",cex=1.5)
text(3620,16,"L",cex=1.5)
text(4615,16,"L",cex=1.5)
text(8642,9,"S",cex=1.5)
text(8980,12,"M",cex=1.5)
lines(lD[1:15000,2]-30,lwd=0.5)
text(250,10-30,"Site 2",cex=2)
lines(lD[1:15000,3]-60,lwd=0.5)
text(250,10-60,"Site 3",cex=2)
lines(lD[1:15000,4]-90,lwd=0.5)
text(250,10-90,"Site 4",cex=2)
segments(250,-100,250+1500,-100,lwd=2)
text(1000,-105,"100 ms",cex=2)
dev.off()

png(file="Pouzat_Fig_3.png",width=1000,height=600)
debut <- 15.35*15000
fin <- 15.55*15000
plot(c(debut,fin),c(-105,15),type="n",axes = FALSE,
     xlab="",ylab="")
lines(debut:fin,lD[debut:fin,1],lwd=0.5)
text(15.4*15000,12,"A",cex=1.5)
text(15.42*15000,12,"B",cex=1.5)
#text(debut+250,10,"Site 1",cex=2)
lines(debut:fin, lD[debut:fin,2]-30,lwd=0.5)
text(15.505*15000,-21,"C",cex=1.5)
text(15.52*15000,-21,"D",cex=1.5)
#text(debut+250,10-30,"Site 2",cex=2)
lines(debut:fin, lD[debut:fin,3]-60,lwd=0.5)
#text(debut+250,10-60,"Site 3",cex=2)
lines(debut:fin, lD[debut:fin,4]-90,lwd=0.5)
#text(debut+250,10-90,"Site 4",cex=2)
segments(debut+250,-100,debut+250+300,-100,lwd=2)
text(debut+250+75,-105,"20 ms",cex=2)
dev.off()

png(file="Pouzat_Fig_4.png",width=1000,height=500) 
plot(c(0,15000),c(-15,15),type="n",axes = FALSE,
     xlab="",ylab="")
lines(lD[1:15000,1],lwd=0.5)
abline(h=7,lty=2)
abline(h=11,lty=2,col="orangered")
abline(h=15,lty=2,col="blue")
dev.off()

png(file="Pouzat_Fig_5.png",width=1000,height=600)
layout(matrix(c(1,2:4,1,5:7),nr=2,byrow = TRUE))
par(mar=c(0,0,0,0),cex=2)
plot(median(evtsE_noj[[7]])[136:180],type="l",
     ylim = c(-10,10),axes = FALSE,
     xlab="",ylab="",lwd=2)
text(3,10,"A",cex=2)
lines(median(evtsE_noj[[6]])[136:180],
      col=2,lwd = 2)
segments(30,-8,37.5,-8,lwd=2)
text(34,-9,"0.5 ms",cex=1)
segments(30,-8,30,-7,lwd=2)
text(30,-7.5,expression("1 x"~sigma[noise]),pos=2,cex=1)

plot(median(evtsE_noj[[1]])[1:45],type="l",
     ylim = c(-15,15),axes = FALSE,
     xlab="",ylab="",lwd=3,col="grey70")
text(3,14,"B",cex=2)
##text(22,17,"Événement et motif 1",cex=0.5)
lines(evtsE_noj[[2]][,2][1:45],
      lwd=2)

plot(median(evtsE_noj[[2]])[1:45],type="l",
     ylim = c(-15,15),axes = FALSE,
     xlab="",ylab="",lwd=3,col="grey70")
lines(evtsE_noj[[2]][,2][1:45],
      lwd=2)

plot(median(evtsE_noj[[3]])[1:45],type="l",
     ylim = c(-15,15),axes = FALSE,
     xlab="",ylab="",lwd=3,col="grey70")
lines(evtsE_noj[[2]][,2][1:45],
      lwd=2)

plot(evtsE_noj[[2]][,2][1:45]-median(evtsE_noj[[1]])[1:45],
     type="l",
     ylim = c(-15,15),axes = FALSE,
     xlab="",ylab="",lwd=3,col="grey50")
text(30,-8,paste(round(sum((evtsE_noj[[2]][,2][1:45]-median(evtsE_noj[[1]])[1:45])^2))),
     cex=1,col=2)

plot(evtsE_noj[[2]][,2][1:45]-median(evtsE_noj[[2]])[1:45],
     type="l",
     ylim = c(-15,15),axes = FALSE,
     xlab="",ylab="",lwd=3,col="grey50")
segments(30,8,37.5,8,lwd=2)
text(34,6,"0.5 ms",cex=1)
segments(30,8,30,11,lwd=2)
text(30,9,expression("3 x"~sigma[noise]),pos=2,cex=1)
text(30,-8,paste(round(sum((evtsE_noj[[2]][,2][1:45]-median(evtsE_noj[[2]])[1:45])^2))),
     cex=1,col=2)

plot(evtsE_noj[[2]][,2][1:45]-median(evtsE_noj[[3]])[1:45],
     type="l",
     ylim = c(-15,15),axes = FALSE,
     xlab="",ylab="",lwd=3,col="grey50")
text(30,-8,paste(round(sum((evtsE_noj[[2]][,2][1:45]-median(evtsE_noj[[3]])[1:45])^2))),
     cex=1,col=2)
dev.off()

png(file="Pouzat_Fig_6.png",width=800,height=800)
quatreSurUn <- cbind(unclass(evtsE_noj[[1]])[1:45,],
                     unclass(evtsE_noj[[2]])[1:45,],
                     unclass(evtsE_noj[[3]])[1:45,],
                     unclass(evtsE_noj[[7]])[1:45,])


layout(matrix(1:4,nc=2,byrow=TRUE))
par(mar=c(4,4,1,1))
matplot(quatreSurUn,type='l',col=1,lwd=0.5,lty=1,
        axes=FALSE,xlab="",ylab = "")
mtext("A",side=2,at=17.5,cex=1.5,las=1,line=2)
abline(v=6,lwd=2,col="grey70")
text(6,-17,expression("400"~mu*s),pos=4,cex=1.5)
abline(v=18,lwd=2,col="grey70")
text(18,-17,"1.2 ms",pos=4,cex=1.5)

plot(quatreSurUn[6,],quatreSurUn[18,],pch=19,
     xlab=expression("Amplitude at 400"~mu*s),
     ylab="Amplitude at 1.2 ms",cex.lab=1.5)
mtext("B",side=2,at=15.5,las=1,cex=1.5,line=2)
un <- (1:dim(quatreSurUn)[2])[quatreSurUn[6,] <= -11]
deux <- (1:dim(quatreSurUn)[2])[quatreSurUn[6,] > -11 & quatreSurUn[18,] > 8]
trois <- (1:dim(quatreSurUn)[2])[quatreSurUn[6,] > -11 & quatreSurUn[18,] <=  8 & quatreSurUn[18,] > 5+quatreSurUn[6,]]
quatre <- (1:dim(quatreSurUn)[2])[quatreSurUn[18,] <= 5+quatreSurUn[6,]]
couleurs <- rep("black",dim(quatreSurUn)[2])
couleurs[deux] <- "brown4"
couleurs[trois] <- "royalblue4"
couleurs[quatre] <- "orangered"
plot(quatreSurUn[6,],quatreSurUn[18,],pch=19,col=couleurs,
     xlab=expression("Amplitude at 400"~mu*s),
     ylab="Amplitude at 1.2 ms",cex.lab=1.5)
mtext("C",side=2,at=15.5,las=1,cex=1.5,line=2)
abline(v=-11)
segments(-11,8,5,8)
abline(5,1)

matplot(quatreSurUn,type='l',col=couleurs,lwd=0.5,lty=1,
        axes=FALSE,xlab="",ylab = "")
mtext("D",side=2,at=17.5,cex=1.5,las=1,line=2)
segments(25,-15,32.5,-15,lwd=2)
text(29,-17,"0.5 ms",cex=2)
segments(25,-15,25,-10,lwd=2)
text(25,-12.5,expression("5 x"~sigma[noise]),pos=2,cex=2)
dev.off()

png(file="Pouzat_Fig_7.png",width=800,height=600)
quatreSurUn.pc <- prcomp(t(quatreSurUn))
layout(matrix(1:2,nc=2,byrow=TRUE))
par(mar=c(5,5,1,1))
moyenne <- quatreSurUn.pc$center
moyennePC1 <- moyenne - 10*quatreSurUn.pc$rotation[,1]
moyennePC2 <- moyenne + 10*quatreSurUn.pc$rotation[,2]
plot(moyenne,type="l",ylim=range(quatreSurUn),
     axes=FALSE,xlab="",ylab = "",lwd=3)
lines(moyennePC1,col=2,lwd=3)
lines(moyennePC2,col=4,lwd=3)
mtext("A",side=2,at=17.5,cex=1.5,las=1,line=2)
legend(1,17.5,c("Mean"),
       col=c(1),lwd=3,bty = "n",cex=2)
legend(1,-10.5,c("Mean + 10 x PC1", "Mean + 10 x PC2"),
       col=c(2,4),lwd=3,bty = "n",cex=2)


plot(quatreSurUn.pc$x[,2],-quatreSurUn.pc$x[,1],pch=19,cex.lab=1.5,
     xlab="Projection onto PC2",ylab = "Projection onto PC1")
mtext("B",side=2,at=35,las=1,cex=1.3,line=2)
abline(h=0)
segments(0,0,0,40)
abline(a=-17,b=1.1)
dev.off()

png(file="Pouzat_Fig_8.png",width=700,height=400)
predAll <- ts(pred0+pred1+pred2+pred3+pred4,start=0,freq=15e3)
superList <- round0[sapply(round0, function(l) 209100 < l[2] & l[2] < 209190)]
cA <- centers[[superList[[1]][[1]]]]$center[1:130]+
    superList[[1]][[3]]*centers[[superList[[1]][[1]]]]$centerD[1:130]+
    superList[[1]][[3]]^2/2*centers[[superList[[1]][[1]]]]$centerDD[1:130]
tA <- centers[[superList[[1]][[1]]]]$center_idx+superList[[1]][[2]]-209099
gA <- 0 < tA & tA <= 90
cB <- centers[[superList[[2]][[1]]]]$center[1:130]+
    superList[[2]][[3]]*centers[[superList[[2]][[1]]]]$centerD[1:130]+
    superList[[2]][[3]]^2/2*centers[[superList[[2]][[1]]]]$centerDD[1:130]
tB <- centers[[superList[[2]][[1]]]]$center_idx+superList[[2]][[2]]-209099
gB <- 0 < tB & tB <= 90
layout(matrix(1:2,nc=2))
par(mar=c(1,1,1,1))
plot(lD[209100:209190,1],type="l",lwd=3,col="black",
     axes = FALSE, xlab = "",ylab = "",ylim = range(lD[209100:209190,1]))
lines(cA[gA]+cB[gB],col="grey50",lwd=2)
text(1,15.5,"A",cex=2)
legend(40,-9,c("Data","Neurons 1 \n& 2"),
       col=c("black","grey50"),pch=15,bty="n",cex=1.5)

plot(tA[gA],cA[gA],col=2,lwd=3,type="l",
     axes = FALSE, xlab = "",ylab = "",ylim = range(lD[209100:209190,1]))
lines(tB[gB],cB[gB],col=4,lwd=3)
text(1,15.5,"B",cex=2)
segments(50,-15,65,-15,lwd=2)
text(57.5,-17,"1 ms",cex=2)
segments(50,-15,50,-10,lwd=2)
text(50,-12.5,expression("5 x"~sigma[noise]),pos=4,cex=2)
legend(50,15.5,c("Neuron 1","Neuron 2"),
       col=c(2,4),pch=15,bty="n",cex=1.5)
dev.off()

png(file="Pouzat_Fig_9.png",width=1000,height=700)
layout(matrix(c(1,2,3,3),nc=2))
par(mar=c(1,1,2,1),cex=2)
domaine <- 7501:15000
plot(PKref[1,domaine],type="l",axes=FALSE,
     xlab = "",ylab="",lwd=1,
     main="Cell attached")
segments(1,-1.2,1+50*15,-1.2,lwd=2)
text(25*15,-1.5,"50 ms")
text(1,1.5,"A")
plot(pkD[domaine,3],type="l",axes=FALSE,
     xlab = "",ylab="",lwd=1,
     main="Extracellular")
text(1,0.71,"B")

PKrefRec <- PKref[1,]
PKrefRec[PKrefRec<1] <- 0
SPref <- peaks(PKrefRec,25)
debut_bouffée <- diff(SPref[-489]) > 5000 & diff(SPref[-1]) < 500
debut_bouffée <- c(TRUE,debut_bouffée,FALSE)
SPref2 <- SPref[debut_bouffée]
PKevtsBouffées <- mkEvents(SPref2,PKref[1,],25,500)
plot(PKevtsBouffées,evts.lwd=0.5,medAndMad=FALSE)
title(main="489 bursts aligned \n on their first AP")
text(1,1.65,"C")
segments(60,-1.2,60+5*15,-1.2,lwd=2)
text(60+2.5*15,-1.45,"5 ms")
dev.off()

png(file="Pouzat_Fig_10.png",width=1000,height=600)
domaine_super <- c(1:130,391:520)
F7 <- centers[[7]]$center[domaine_super]-sum(centers[[7]]$center[domaine_super])/260
F7 <- F7/sum(centers[[7]]$center[domaine_super]*F7)
F2 <- centers[[2]]$center[domaine_super]-sum(centers[[2]]$center[domaine_super])/260
F2 <- F2/sum(centers[[2]]$center[domaine_super]*F2)
dist2_7 <- outer(spike_trains[[2]],spike_trains[[7]],function(x,y) abs(x-y))
which(5 <= dist2_7 & dist2_7 <= 10)
which(5 <= dist2_7 & dist2_7 <= 10) %/% length(spike_trains[[2]])
which(5 <= dist2_7 & dist2_7 <= 10) %% length(spike_trains[[2]])
spike_trains[[2]][which(5 <= dist2_7 & dist2_7 <= 10) %% length(spike_trains[[2]])]
spike_trains[[7]][which(5 <= dist2_7 & dist2_7 <= 10) %/% length(spike_trains[[2]]) + 1]
predAll <- ts(pred0+pred1+pred2+pred3+pred4,start=0,freq=15e3)
le_centre <- round(spike_trains[[2]][which(5 <= dist2_7 & dist2_7 <= 10) %% length(spike_trains[[2]])])
s1 <- c(numeric(130),lD[(le_centre-130):(le_centre+130),1],numeric(130))
s4 <- c(numeric(130),lD[(le_centre-130):(le_centre+130),4],numeric(130))
s1F2 <- filter(s1,rev(F2[1:130]),method="convolution",circular = TRUE)
s4F2 <- filter(s4,rev(F2[-(1:130)]),method="convolution",circular = TRUE)
s1_4F2 <- s1F2+s4F2
s1F7 <- filter(s1,rev(F7[1:130]),method="convolution",circular = TRUE)
s4F7 <- filter(s4,rev(F7[-(1:130)]),method="convolution",circular = TRUE)
s1_4F7 <- s1F7+s4F7
noiseEL <- mkNoise(sp0,lD,49,80,safetyFactor=2,2000)
F2limites <- quantile(apply(unclass(noiseEL)[domaine_super,],2,function(x) sum((x+centers[[2]]$center[domaine_super])*F2)),c(0.005,0.995))
F7limites <- quantile(apply(unclass(noiseEL)[domaine_super,],2,function(x) sum((x+centers[[7]]$center[domaine_super])*F7)),c(0.005,0.995))
layout(matrix(1:6,nc=3,nr=2))
par(mar=c(1,1,1,1),cex=2)
plot(s1[131:(length(s1)-130)],type="l",lwd=2,
     axes=FALSE,xlab="",ylab="",
     ylim=range(c(s1,s4)))
text(0,10,"Site 1",pos=4)
text(0,15,"A")
segments(10,-10,41,-10,lwd=3)
text(26,-12,"2 ms")
plot(s4[131:(length(s4)-130)],type="l",lwd=2,
     axes=FALSE,xlab="",ylab="",
     ylim=range(c(s1,s4)))
text(0,10,"Site 4",pos=4)
plot(s1_4F2[131:(length(s1_4F2)-130)],
     axes=FALSE,xlab="",ylab="",type="n",
     ylim=range(c(s1_4F2,s1_4F7,F2limites,F7limites)))
text(0,0.4,"filter 2\n output",pos=4)
text(-2,1.1,"B")
rect(0,F2limites[1],260,F2limites[2],col="grey70",border=NA,density=50)
lines(s1_4F2[131:(length(s1_4F2)-130)],lwd=2,col=2)
plot(s1_4F7[131:(length(s1_4F7)-130)],
     axes=FALSE,xlab="",ylab="",type="n",
     ylim=range(c(s1_4F2,s1_4F7,F2limites,F7limites)))
text(0,0.4,"filter 7\n output",pos=4)
rect(0,F7limites[1],260,F7limites[2],col="grey70",border=NA,density=50)
lines(s1_4F7[131:(length(s1_4F7)-130)],lwd=2,col=4)

p2 <- round(spike_trains[[2]][which(5 <= dist2_7 & dist2_7 <= 10) %% length(spike_trains[[2]])])
p7_1 <- round(spike_trains[[7]][which(5 <= dist2_7 & dist2_7 <= 10) %/% length(spike_trains[[2]]) + 1])-le_centre
#p7_2 <- round(spike_trains[[7]][which(5 <= dist2_7 & dist2_7 <= 10) %/% length(spike_trains[[2]]) + 2])-le_centre
p7_2 <- 72
plot(s1[131:(length(s4)-130)],type="l",lwd=2,
     axes=FALSE,xlab="",ylab="",
     ylim=range(c(s1,s4)))
text(-20,7,"Resolution\n site 1",pos=4)
text(0,15,"C")
lines(centers[[2]]$center_idx+130,centers[[2]]$center[1:130],col=2,lwd=2)
lines(centers[[7]]$center_idx+130+p7_1,centers[[7]]$center[1:130],col=4,lwd=2)
lines(centers[[7]]$center_idx+130+p7_2,centers[[7]]$center[1:130],col=4,lwd=2)
plot(s4[131:(length(s4)-130)],type="l",lwd=2,
     axes=FALSE,xlab="",ylab="",
     ylim=range(c(s1,s4)))
text(-20,7,"Resolution\n site 4",pos=4)
lines(centers[[2]]$center_idx+130,centers[[2]]$center[391:520],col=2,lwd=2)
lines(centers[[7]]$center_idx+130+p7_1,centers[[7]]$center[391:520],col=4,lwd=2)
lines(centers[[7]]$center_idx+130+p7_2,centers[[7]]$center[391:520],col=4,lwd=2)
dev.off()

pkD <- ts(pkD,start=0,freq=15e3)
##pkDd <- apply(pkD,2,function(x) c(0,diff(x,2),0))
##pkDd <- ts(pkDd,start=0,freq=15e3)
##pkDd.mad <- apply(pkDd,2,mad)
##pkDd <- t(t(pkDd)/pkDd.mad)
##pkDd <- t(t(pkDd)-apply(pkDd,2,median))
##pkDd <- ts(pkDd,start=0,freq=15e3)
##stereo2 <- pkDd[,3:4]
stereo2 <- pkD[,3:4]
stereo2.mad <- apply(stereo2,2,mad)
stereo2 <- t(t(stereo2)/stereo2.mad)
##rm(pkD,pkDd)
stereo2r <- stereo2
stereo2r <- filter(stereo2r,rep(1,5)/5)
stereo2r[is.na(stereo2r)] <- 0
stereo2r.mad <- apply(stereo2r,2,mad)
stereo2r <- t(t(stereo2r)/stereo2r.mad)
##stereo2r[stereo2r < 3.5] <- 0
stereo2r[stereo2r < 7.5] <- 0
##stereo2r <- ts(stereo2r,start=0,freq=15e3)
sp2 <- peaks(apply(stereo2r,1,sum),30)
rm(stereo2r)
evtsPK <- mkEvents(sp2,stereo2,14,25)
##goodEvtsPK <- goodEvtsFct(evtsPK,8)
##les_angles <- apply(evtsPK[,goodEvtsPK],2,function(x) atan(coef(lm(x[41:80] ~ x[1:40]-1))[1]))*90/pi
##les_angles <- atan(unclass(evtsPK)[55,goodEvtsPK]/unclass(evtsPK)[15,goodEvtsPK])*90/pi
##les_angles <- apply(evtsPK[,goodEvtsPK],2,function(x) atan(coef(lm(x[50:60] ~ x[10:20]-1))[1]))*90/pi
les_angles <- apply(evtsPK,2,function(x) atan(coef(lm(x[10:20] ~ x[50:60]-1))[1]))*180/pi
les_classes <- integer(length(les_angles))
#frontieres_des_classes <- c(0,8,17,23.36,29.95)
frontieres_des_classes <- c(0,28.6,42.7,51.5,73.3)
les_classes[frontieres_des_classes[1] > les_angles] <- 6
les_classes[frontieres_des_classes[1] <= les_angles & les_angles < frontieres_des_classes[2]] <- 2
les_classes[frontieres_des_classes[2] <= les_angles & les_angles < frontieres_des_classes[3]] <- 3
les_classes[frontieres_des_classes[3] <= les_angles & les_angles < frontieres_des_classes[4]] <- 1
les_classes[frontieres_des_classes[4] <= les_angles & les_angles < frontieres_des_classes[5]] <- 4
les_classes[frontieres_des_classes[5] <= les_angles] <- 5
mes_couleurs <- c("black","royalblue","brown4","orangered","sienna3","white")

png(file="Pouzat_Fig_11.png",width=800,height=400)
par(mar=c(1,1,1,1))
domaine_stereo <- (2.3*15e3):(2.5*15e3)
plot(domaine_stereo,stereo2[domaine_stereo,1],type="n",
     axes=FALSE,xlab="",ylab="",ylim = c(-90,30))
pos_bouffee1 <- c(2.379636,2.386250)
pos_bouffee2 <- c(2.346563, 2.398377, 2.436410, 2.469115)
pos_fond <- c(2.333334,2.370633,2.405175, 2.454600)
abline(v=pos_bouffee1*15e3,lwd=2,col="grey50")
abline(v=pos_bouffee2*15e3,lwd=2,col="grey50",lty=2)
lines(domaine_stereo,stereo2[domaine_stereo,1],lwd=3)
lines(domaine_stereo,stereo2[domaine_stereo,2]-60,lwd=3)
segments(2.31*15e3,-80,2.33*15e3,-80,lwd=3)
text(2.32*15e3,-85,"20 ms",cex=2)
text(2.32*15e3,15,"Site 1",cex=2)
text(2.32*15e3,-45,"Site 2",cex=2)
dev.off()

png(file="Pouzat_Fig_12.png",width=800,height=400)
layout(matrix(c(1,2),nc=2))
par(mar=c(5,5,1,1))
plot(stereo2[sp2,2],stereo2[sp2,1],cex.lab=2,xlim=c(0,25),ylim=c(0,45),
     type="n",xlab="Amplitude on site 2",ylab="Amplitude on site 1")
sapply(frontieres_des_classes[-1],
       function(θ) segments(0,0,25,25*tan(θ*pi/180),lwd=2,col="grey50"))
segments(0,0,0,45,lwd=2,col="grey50")
points(stereo2[sp2,2],stereo2[sp2,1],pch=".")
hist(les_angles,breaks=200,xlim=c(0,90),xlab=expression(theta),
     main="",prob=TRUE,ylab="Estimated density",cex.lab=2)
abline(v=frontieres_des_classes[-1],lwd=2,col="grey50")
dev.off()

png(file="Pouzat_Fig_13.png",width=800,height=400) 
layout(matrix(c(1,2,2,2),nc=4))
par(mar=c(5,5,1,1))
assez_grand <- stereo2[sp2,1] > 10 | stereo2[sp2,2] > 10
plot(stereo2[sp2[assez_grand],2],stereo2[sp2[assez_grand],1],cex.lab=2,xlim=c(0,25),ylim=c(0,45),
     type="n",xlab="Amplitude on site 2",ylab="Amplitude on site 1")
points(stereo2[sp2[assez_grand],2],stereo2[sp2[assez_grand],1],cex=2,
     col = mes_couleurs[les_classes[assez_grand]],pch=".")
mtext("A",2,at=43,las=1,cex=2,line=2)
domaine_stereo <- (2.3*15e3):(2.5*15e3)
les_PAs <- unclass(sp2)
les_PAs_idx <- (1:length(les_PAs))[(domaine_stereo[1] <= les_PAs &
                                    les_PAs <= domaine_stereo[length(domaine_stereo)]) &
                                   (stereo2[les_PAs,1] > 10 | stereo2[les_PAs,2] > 10)] 
plot(domaine_stereo,stereo2[domaine_stereo,1],type="n",
     axes=FALSE,xlab="",ylab="",ylim = c(-90,30))
mtext("B",2,at=25,las=1,cex=2)
lines(domaine_stereo,stereo2[domaine_stereo,1],lwd=3,col="grey50")
lines(domaine_stereo,stereo2[domaine_stereo,2]-60,lwd=3,col="grey50")
sapply(les_PAs_idx,
       function(idx) {
           origine <- les_classes[idx]
           position <- les_PAs[idx]
           lines((position-14):(position+30),
                 stereo2[(position-14):(position+30),1],
                 col=mes_couleurs[origine],
                 lwd = 3)
           lines((position-14):(position+30),
                 stereo2[(position-14):(position+30),2]-60,
                 col=mes_couleurs[origine],
                 lwd = 3)})
segments(2.31*15e3,-80,2.33*15e3,-80,lwd=3)
text(2.32*15e3,-85,"20 ms",cex=2)
text(2.32*15e3,15,"Site 1",cex=2)
text(2.32*15e3,-45,"Site 2",cex=2)
dev.off()

png(file="Pouzat_Fig_14.png",width=800,height=400)
ideal2 <- splinefun(centers[[2]][["center_idx"]],centers[[2]][["center"]][1:130])
dense <- (-150:300)/10
xx1 <- seq(-14,30,len=31)
xx2 <- xx1-0.75

layout(matrix(c(1,1,2,3),nc=2))

## Graphe 1
par(mar=c(1,1,1,1))
plot(dense,ideal2(dense),type="n",
     axes=FALSE,xlab="",ylab="",
     ylim=c(-45,15),xlim = c(-9,25))
abline(v=xx2,col="grey70")
lines(dense,ideal2(dense),lwd=2)
points(xx1,ideal2(xx1),pch=16,col=1,cex=2)
lines(dense,ideal2(dense)-30,lwd=2)
points(xx2,ideal2(xx2)-30,pch=16,col=1,cex=2)
text(-9,14,"A",cex=3)
segments(17.5,-15,25,-15,lwd=3)
text(21,-17,"0.5 ms",cex=2)
segments(17.5,-15,17.5,-10,lwd=3)
text(17.5,-12.5,expression("5 x"~sigma[noise]),pos=2,cex=2)

set.seed(20110928)
centre <- ideal2(xx2)
le_filtre <- centre - mean(centre)
le_filtre <- le_filtre/sum(le_filtre*centre)
filtre_longueur <- length(le_filtre)
nrep <- 1000
simu1 <- sapply(1:nrep,
                function(rep_idx) {
                    gigue <- runif(1,-0.75,0.75)
                    bruit <- rnorm(filtre_longueur)
                    evtA <- ideal2(xx2+gigue)+bruit
                    evtB <- centre+bruit
                    c(sum(bruit^2),
                      sum((evtA-centre)^2),
                      sum(evtB*le_filtre),
                      sum(evtA*le_filtre))})
## Graphe 2 motifs
par(mar=c(1,5,1,1))
plot(xx1,ideal2(xx1)-ideal2(xx2),type="n",
     axes = FALSE, xlab="",ylab = "",
     ylim = c(-10,10),xlim=c(-9,25))
abline(h=0,col="grey70",lwd=2)
lines(xx1,ideal2(xx1)-ideal2(xx2),lwd=3)
points(xx1,ideal2(xx1)-ideal2(xx2),pch=16,cex=2)
text(-9,9,"B1",cex=2)
segments(0,-8,7.5,-8,lwd=3)
text(3.5,-9.5,"0.5 ms",cex=2)
segments(0,-8,0,-3,lwd=3)
text(0,-5,expression("5 x"~sigma[noise]),pos=2,cex=2)

par(mar=c(5,5,1,1),cex=1.2)
plot(sort(simu1[2,]),(1:nrep)/nrep,
     type="l",xlab="Residual sum of squares",
     ylab="Probability",lwd = 3,col=2)
lines(sort(simu1[1,]),(1:nrep)/nrep,lwd = 3,col=1)
legend(100,0.6,c("Without jitter","With jitter"),col=c(1,2),lwd=3,bty="n")
mtext("B2",side=2,at=1,cex=1.5,las=1,line=2)
dev.off()

png(file="Pouzat_Fig_15.png",width=800,height=800)
layout(matrix(cbind(1:3,c(0,4,5),c(0,0,6)),nc=3))
par(mar=c(2,2,1,1))
plot(evtsE.pc$x[,1],evtsE.pc$x[,2],pch=".",axes = FALSE)
box()
mtext("CP1",side=1,line=1,cex=1.5)
mtext("CP2",side=2,line=0.4,cex=1.5)
plot(evtsE.pc$x[,1],evtsE.pc$x[,3],pch=".",axes = FALSE)
box()
mtext("CP1",side=1,line=1,cex=1.5)
mtext("CP3",side=2,line=0.4,cex=1.5)
plot(evtsE.pc$x[,1],evtsE.pc$x[,4],pch=".",axes = FALSE)
box()
mtext("CP1",side=1,line=1,cex=1.5)
mtext("CP4",side=2,line=0.4,cex=1.5)
plot(evtsE.pc$x[,2],evtsE.pc$x[,3],pch=".",axes = FALSE)
box()
mtext("CP2",side=1,line=1,cex=1.5)
mtext("CP3",side=2,line=0.4,cex=1.5)
plot(evtsE.pc$x[,2],evtsE.pc$x[,4],pch=".",axes = FALSE)
box()
mtext("CP2",side=1,line=1,cex=1.5)
mtext("CP4",side=2,line=0.4,cex=1.5)
plot(evtsE.pc$x[,3],evtsE.pc$x[,4],pch=".",axes = FALSE,
     xlab="CP3",ylab="CP4")
box()
mtext("CP3",side=1,line=1,cex=1.5)
mtext("CP4",side=2,line=0.4,cex=1.5)
dev.off()

png(file="Pouzat_Fig_16.png",width=800,height=800)
garde7 <- c10b != 5 & c10b != 7 & c10b != 10
how_many <- 7
X <- evtsE.pc$x[,1][garde7]
Y <- evtsE.pc$x[,3][garde7]
origine <- c10b[garde7]
echant_idx <- 1:length(X)
library(colorspace)
my_pallette <- rainbow_hcl(7)

layout(matrix(1:4,nc=2,byrow=TRUE))
par(cex=2,mar=c(1,1,1,1))

plot(X,Y,pch=3,
     col="black",
     xlab="",
     ylab="",axes=FALSE,
     main="Initialization",
     xlim=range(X),ylim=range(Y))
##text(-23,18,"A",cex=0.75)
segments(-15,-20,-5,-20,lwd=2)
text(-10,-23,"PC1",cex=0.5)
segments(-15,-20,-15,-15,lwd=2)
text(-22,-17.5,"PC3",cex=0.5)
idx0 <- c(111,188,228,46,308,133,35)+1
points(X[idx0],Y[idx0],pch=16,col=my_pallette,cex=1)

classe0 <- sapply(echant_idx,
                  function(idx) which.min(sqrt((X[idx0]-X[idx])^2+(Y[idx0]-Y[idx])^2)))

plot(X,Y,pch=3,
     xlab="",
     ylab="",axes=FALSE,
     main="Assignment",
     col=my_pallette[classe0],
     xlim=range(X),ylim=range(Y))
##text(-23,18,"B",cex=0.75)
segments(-15,-20,-5,-20,lwd=2)
text(-10,-23,"PC1",cex=0.5)
segments(-15,-20,-15,-15,lwd=2)
text(-22,-17.5,"PC3",cex=0.5)

plot(X,Y,pch=3,
     xlab="",
     ylab="",axes=FALSE,
     main="New centroids",
     col=my_pallette[classe0],
     xlim=range(X),ylim=range(Y))
##text(-23,18,"C",cex=0.75)
segments(-15,-20,-5,-20,lwd=2)
text(-10,-23,"PC1",cex=0.5)
segments(-15,-20,-15,-15,lwd=2)
text(-22,-17.5,"PC3",cex=0.5)
points(X[idx0],Y[idx0],pch=1,col=1)
X_1 <- sapply(1:how_many, function(idx) mean(X[classe0==idx]))
Y_1 <- sapply(1:how_many, function(idx) mean(Y[classe0==idx]))
points(X_1,Y_1,pch=16,col=1,cex=1)

mk_kmeans_step <- function(current_pos,Data=cbind(X,Y)) {
    ne <- nrow(Data)
    nclass <- nrow(current_pos)
    new_c <- apply(Data,1,
                   function(x) which.min(apply((matrix(x,nr=nclass,nc=2,byrow=TRUE)-current_pos)^2,1,sum)))
    new_pos <- sapply(1:nclass,
                      function(c_idx) apply(Data[new_c == c_idx,],2,mean))
    t(new_pos) }

trajectory <- array(0,c(how_many,2,20))
trajectory[,,1] <- cbind(X[idx0],Y[idx0])
for (i in 2:20) trajectory[,,i] <- mk_kmeans_step(trajectory[,,i-1])
final_class <- apply(cbind(X,Y),1,
                     function(x)
                     which.min(apply((matrix(x,nr=how_many,nc=2,byrow=TRUE)-trajectory[,,20])^2,1,sum)))

plot(X,Y,pch=3,
     xlab="",
     ylab="",axes=FALSE,
     main="Results",
     col=my_pallette[final_class],
     xlim=range(X),ylim=range(Y))
##text(-23,18,"D",cex=0.75)
segments(-15,-20,-5,-20,lwd=2)
text(-10,-23,"PC1",cex=0.5)
segments(-15,-20,-15,-15,lwd=2)
text(-22,-17.5,"PC3",cex=0.5)
points(X[idx0],Y[idx0],pch=1,col=1)
apply(trajectory,1,
      function(y) sapply(1:(ncol(y)-1),
                         function(idx) segments(y[1,idx],y[2,idx],y[1,idx+1],y[2,idx+1],lwd=2)))
dev.off()
