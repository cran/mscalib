#testing massvector
rm(list=ls())
source("all.R")
#Load a set of spectra.
mvl <- readBruker(massvectorlist(),"D:/DATA/020529r_dt_Pirellula/")
#Hist of massvector lengths.
hist(mvl)
#Hist of mass Frequencies.
#get peptide masses discarding chemical noise
print("1")
mvlf <- WsFilter(mvl,peptides = TRUE)
print("2")
#Recalibration.
res <- getrecalib(mvlf) # get recalibration model for filtered data
print("3")
plot(res)
print("4")
# Refining the calibration model
# Select Calibration objects with a PQM > 10.
res2<-subset(res,PQM>10)
print("5")
#Applying the calibration.

mvlfr<-applyrecalib(mvl, res2)
print("6")
plot(mvlf,xlim=c(1712,1714))
print("7")
plot(mvlfr,xlim=c(1712,1714),add=TRUE,col=3)


#Splitting the collection
#A part for which the calibration worked well.
print("8")
tmvlfr<-mvlfr[names(subset(res,PQM>10))]
#A part for which it worked worse.
print("9")
tmvlfr2<-mvlfr[names(subset(res,PQM<=10))]
print("10")
plot(tmvlfr2,xlim=c(1712,1714))
print("11")
#stop()
##Try for this part of the collection an affine calibration.
tmvlfr2 <- correctinternal(tmvlfr2,calib)
print("12")
plot(tmvlfr2,add=T,col=3)
print("13")
mvlfA<-c(tmvlfr,tmvlfr2)

#remove nod needed massvecotrlists
rm(tmvlfr,tmvlfr2,mvlf,mvlfr)

print("14")
plot(mvl,xlim=c(1712.5,1713.2))
print("15")
plot(mvlfA[names(mvl)],add=T,col=3)


#External Calibration.
#read calibration samples
source("all.R")
print("16")
ppg<-readBruker(massvectorlist(),"D:/DATA/020529r_dt_Pirellula_ppg/")
print("17")
image(ppg,what="lengthmv")
print("18")
tmp <- getextcalib(ppg,error=200)
plot(tmp)
print("19")
mvlfB <- applyextcalib(mvlfA,tmp)
print("20")
plot(mvl,xlim=c(1712.7,1713.1))
print("21")
plot(mvlfA,add=T,col=3)
print("22")
plot(mvlfB,add=T,col=4)


#get the calibration objects.

# The length of the new calibration vector of abundant masses
# can be limited by the parameter labund.
print("23")
cres<- getglobalcalib(mvlfB,calib,error=150,labund=20)
#the calibration list is even longer.
print("24")
length(calib)
# For comparison we perform a normal internal calibrarion.
print("25")
ires<- getintcalib(mvlfB,calib,error=250)
# We compare the number of peaks matched in every peaklist
# be the global and internal calibration.
print("26")
summary(as.data.frame(cres)$nrmatch)
summary(as.data.frame(ires)$nrmatch)

print("27")
cres<- getglobalcalib(mvlfB,calib,error=150,labund=12)
image(subset(cres,nrmatch>4),what="nrmatch")
# for calibration we use only peaklists where at least
# 4 peaks matched the calibration list.
print("28")
cres <- subset(cres,nrmatch>4)
cres
#the error model of this "calibintstat" objects are
#applied to correct the error of the massvectors.
mvlfC <- applyintcalib(mvlfB,cres)
print("29")
plot(mvl,xlim=c(1712.7,1713.1))
plot(mvlfA,add=T,col="green")
plot(mvlfB,add=T,col="blue")
plot(mvlfC,add=T,col="darkorange")
