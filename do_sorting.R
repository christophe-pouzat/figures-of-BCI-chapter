reposName <-"http://xtof.disque.math.cnrs.fr/data/"
dN <- paste("Locust_",1:4,".dat.gz",sep="")
sapply(1:4, function(i)
       download.file(paste(reposName,dN[i],sep=""),
                     dN[i],mode="wb")
       )

dN <- paste("Locust_",1:4,".dat.gz",sep="")
nb <- 20*15000
lD <- sapply(dN,
             function(n) {
                 mC <- gzfile(n,open="rb")
                 x <- readBin(mC,what="double",n=nb)
                 close(mC);x
             }
             )
colnames(lD) <- paste("site",1:4)
dim(lD)

source("http://xtof.perso.math.cnrs.fr/code/sorting_with_r.R")

lD.mad <- apply(lD,2,mad)
lD <- t((t(lD)-apply(lD,2,median))/lD.mad)
lD <- ts(lD,start=0,freq=15e3)

lDf <- filter(lD,rep(1,5)/5)
lDf[is.na(lDf)] <- 0
lDf.mad <- apply(lDf,2,mad)
lDf <- t(t(lDf)/lDf.mad)
thrs <- c(4,4,4,4)
bellow.thrs <- t(t(lDf) < thrs)
lDfr <- lDf
lDfr[bellow.thrs] <- 0
remove(lDf)

(sp0 <- peaks(apply(lDfr,1,sum),15))

remove(lDfr)

(sp0E <- as.eventsPos(sp0[sp0 <= dim(lD)[1]/2]))
(sp0L <- as.eventsPos(sp0[sp0 > dim(lD)[1]/2]))

evtsE <- mkEvents(sp0E,lD,14,30)

summary(evtsE)

noiseE <- mkNoise(sp0E,lD,14,30,safetyFactor=2.5,2000)

summary(noiseE)

goodEvtsFct <- function(samp,thr=3) {
    samp.med <- apply(samp,1,median)
    samp.mad <- apply(samp,1,mad)
    above <- samp.med > 0
    samp.r <- apply(samp,2,function(x) {x[above] <- 0;x})
    apply(samp.r,2,function(x) all(abs(x-samp.med) < thr*samp.mad))
}

goodEvts <- goodEvtsFct(evtsE,8)

evtsE.pc <- prcomp(t(evtsE[,goodEvts]))

sapply(1:15, function(n) sum(diag(cov(t(noiseE))))+sum(evtsE.pc$sdev[1:n]^2)-sum(evtsE.pc$sdev^2))

set.seed(20061001,kind="Mersenne-Twister")
km10 <- kmeans(evtsE.pc$x[,1:3],centers=10,iter.max=100,nstart=100)
c10 <- km10$cluster

cluster.med <- sapply(1:10, function(cIdx) median(evtsE[,goodEvts][,c10==cIdx]))
sizeC <- sapply(1:10,function(cIdx) sum(abs(cluster.med[,cIdx])))
newOrder <- sort.int(sizeC,decreasing=TRUE,index.return=TRUE)$ix
cluster.mad <- sapply(1:10, function(cIdx) {ce <- t(evtsE)[goodEvts,];ce <- ce[c10==cIdx,];apply(ce,2,mad)})
cluster.med <- cluster.med[,newOrder]
cluster.mad <- cluster.mad[,newOrder]
c10b <- sapply(1:10, function(idx) (1:10)[newOrder==idx])[c10]

evtsE_noj = lapply(1:10, function(i) mk_aligned_events(sp0E[goodEvts][c10b==i],lD))

t(sapply(evtsE_noj, function(l) summary(attr(l,"delta"),digits=3)))

centers = lapply(1:10, function(i) mk_center_list(attr(evtsE_noj[[i]],"positions"),lD))
names(centers) = paste("Cluster",1:10)

round0 = lapply(as.vector(sp0),classify_and_align_evt, data=lD,centers=centers)
sum(sapply(round0, function(l) l[[1]] == '?'))

pred0 = predict_data(round0,centers)
lD1 = lD - pred0

lDf <- filter(lD1,rep(1,3)/3)
lDf[is.na(lDf)] <- 0
lDf.mad <- apply(lDf,2,mad)
lDf <- t(t(lDf)/lDf.mad)
bellow.thrs <- t(t(lDf) < thrs)
lDfr <- lDf
lDfr[bellow.thrs] <- 0
remove(lDf)
(sp1 <- peaks(lDfr[,1],15))

round1 = lapply(as.vector(sp1),classify_and_align_evt, data=lD1,centers=centers)
sum(sapply(round1, function(l) l[[1]] == '?'))

pred1 = predict_data(round1,centers)
lD2 = lD1 - pred1
lDf <- filter(lD2,rep(1,3)/3)
lDf[is.na(lDf)] <- 0
lDf.mad <- apply(lDf,2,mad)
lDf <- t(t(lDf)/lDf.mad)
bellow.thrs <- t(t(lDf) < thrs)
lDfr <- lDf
lDfr[bellow.thrs] <- 0
remove(lDf)
(sp2 <- peaks(lDfr[,2],15))

round2 = lapply(as.vector(sp2),classify_and_align_evt, data=lD2,centers=centers)
sum(sapply(round2, function(l) l[[1]] == '?'))

pred2 = predict_data(round2,centers)
lD3 = lD2 - pred2
lDf <- filter(lD3,rep(1,3)/3)
lDf[is.na(lDf)] <- 0
lDf.mad <- apply(lDf,2,mad)
lDf <- t(t(lDf)/lDf.mad)
bellow.thrs <- t(t(lDf) < thrs)
lDfr <- lDf
lDfr[bellow.thrs] <- 0
remove(lDf)
(sp3 <- peaks(lDfr[,3],15))

round3 = lapply(as.vector(sp3),classify_and_align_evt, data=lD3,centers=centers)
sum(sapply(round3, function(l) l[[1]] == '?'))

pred3 = predict_data(round3,centers)
lD4 = lD3 - pred3
lDf <- filter(lD4,rep(1,3)/3)
lDf[is.na(lDf)] <- 0
lDf.mad <- apply(lDf,2,mad)
lDf <- t(t(lDf)/lDf.mad)
bellow.thrs <- t(t(lDf) < thrs)
lDfr <- lDf
lDfr[bellow.thrs] <- 0
remove(lDf)
(sp4 <- peaks(lDfr[,4],15))

round4 = lapply(as.vector(sp4),classify_and_align_evt, data=lD4,centers=centers)
sum(sapply(round4, function(l) l[[1]] == '?'))

pred4 = predict_data(round4,centers)
lD5 = lD4 - pred4
lDf <- filter(lD5,rep(1,3)/3)
lDf[is.na(lDf)] <- 0
lDf.mad <- apply(lDf,2,mad)
lDf <- t(t(lDf)/lDf.mad)
bellow.thrs <- t(t(lDf) < thrs)
lDfr <- lDf
lDfr[bellow.thrs] <- 0
remove(lDf)
(sp5 <- peaks(apply(lDfr,1,sum),15))

round5 = lapply(as.vector(sp5),classify_and_align_evt, data=lD5,centers=centers)
sum(sapply(round5, function(l) l[[1]] == '?'))

round_all = c(round0,round1,round2,round3,round4,round5)
spike_trains = lapply(paste("Cluster",1:10), function(cn) sapply(round_all[sapply(round_all, function(l) l[[1]]==cn)], function(l) l[[2]]+l[[3]]))
names(spike_trains) = paste("Cluster",1:10)

save.image(file="LocustAnalysis.RData")
