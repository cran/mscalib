#This is a sample session demonstrating the usage of the package mscalib.

library(mscalib)
data(mvl)
#mvl <- readBruker( massvectorlist(), "D:/DATA/020529r_dt_Pirellula/" )
#mvl <- readBruker(massvectorlist(), "~/DATA/020529r_dt_Pirellula/" )

                                        #Hist of massvector lengths.
hist(mvl)
plot(mvl)


#Hist of mass Frequencies.
#get peptide masses discarding chemical noise
mvlf <- wsFilter( mvl , peptides = TRUE )
plot(mvlf)
#get nonpeptide masses
res2 <- wsFilter(mvl,peptides = FALSE)
plot(res2,col=2,add=TRUE)
#show mass frequencies
hist(res2)

#Take a look on two peaklists as an clustering example
par(mfrow=c(2,2))
plot(mvl[[3]])
plot(wsFilter(mvl[[3]],peptides=FALSE),add=TRUE,col=2)
plot(hclust(wsdist(mvl[[3]]),method="single"),cex=0.6)
plot(mvl[[232]])
plot(wsFilter(mvl[[232]],peptides=FALSE),add=TRUE,col=2)
plot(hclust(wsdist(mvl[[232]]),method="single"),cex=0.5)
par(mfrow=c(1,1))


plot(mvl[[3]])
plot(wsFilter(mvl[[3]],peptides=FALSE),add=TRUE,col=2)

plot(mvl[[232]])
plot(wsFilter(mvl[[232]],peptides=FALSE),add=TRUE,col=2)



#Recalibration.
res <- getrecalib(mvlf) # get recalibration model for filtered data
plot(res)
res2 <- getrecalib(mvl) # get recalibration model for not filtered data
## Look how the models for filtered and not filtered data differ.
## First convert the caliblists to a dataframe.
dres<-as.data.frame(res)
dres2<-as.data.frame(res2)
par(mfrow=c(1,2))
par(mar=c(4.2,4,1,1))
#Compare only data where 
tmp<-dres$lengthmv!=dres2$lengthmv
plot(dres$Coef.Intercept[tmp],dres2$Coef.Intercept[tmp],xlab="b_0 (filtered)",ylab="b_0 (not filtered)",main="Intercept",pch="*")
plot(dres$Coef.Slope[tmp],dres2$Coef.Slope[tmp],xlab="(1-a)*1e6 (filtered)",ylab="(1-a)*1e6 (not filtered)", main="Slope",pch="*")

par(mfrow=c(1,1))

#look at the calibration results.
res <- getrecalib(mvlf) # get recalibration model for filtered data
res
plot(res)
hist(res)
# To take closer look at the Intercept and PQM
# convert to an R Structure.
dres <- as.data.frame(res)
plot(dres$Coef.Intercept,dres$PQM,xlab="Coef.Intercept",ylab="PQM",pch="*")
# Refining the calibration model
# Select Calibration objects with a PQM > 10.
res2<-subset(res,PQM>10)
res2
image(res2,what="Coef.Slope")
image(res2,what="Coef.Intercept")
plot(res2)


#Applying the calibration.

mvlfr<-applyrecalib(mvl, res2)
tm <- gamasses(mvl,accur=0.3,abund=60)
plot(mvlf,xlim=c(1712,1714))
plot(mvlfr,xlim=c(1712,1714),add=TRUE,col=3)


#Splitting the collection
#A part for which the calibration worked well.
tmvlfr<-mvlfr[names(subset(res,PQM>10))]


#A part for which it worked worse.
tmvlfr2<-mvlfr[names(subset(res,PQM<=10))]
plot(tmvlfr2,xlim=c(1712,1714))

##Try for this part of the collection an affine calibration.
data(cal)
tmvlfr3 <- correctinternal(tmvlfr2,cal,ppm=TRUE)
plot(tmvlfr3,add=T,col=3)

mvlfA<-c(tmvlfr,tmvlfr3)

#remove nod needed massvecotrlists
#rm(tmvlfr,tmvlfr2,mvlf,mvlfr)
plot(mvl,xlim=c(1712.5,1713.2))
plot(mvlfA[names(mvl)],add=T,col=3)


#External Calibration.
#read calibration samples


#ppg<-readBruker(massvectorlist(),"D:/DATA/020529r_dt_Pirellula_ppg/")
#ppg<-readBruker(massvectorlist(),"~/DATA/020529r_dt_Pirellula_ppg/")
data(ppg)
image(ppg,what="lengthmv")
tmp <- getextcalib(ppg,error=200)
plot(tmp)

mvlfB <- applyextcalib(mvlfA,tmp)
plot(mvl,xlim=c(1712.7,1713.1))
plot(mvlfA,add=T,col=3)
plot(mvlfB,add=T,col=4)


# Get the calibration objects.
# The length of the new calibration vector of abundant masses
# can be limited by the parameter labund.
cres<- getglobalcalib(mvlfB,cal,error=150,labund=20)
#the calibration list is even longer.
length(cal)
# For comparison we perform a normal internal calibrarion.
ires<- getintcalib(mvlfB,cal,error=250)
# We compare the number of peaks matched in every peaklist
# be the global and internal calibration.
summary(as.data.frame(cres)$nrmatch)
summary(as.data.frame(ires)$nrmatch)


cres<- getglobalcalib(mvlfB,cal,error=250,labund=8)

image(subset(cres,nrmatch>4),what="nrmatch")
# for calibration we use only peaklists where at least
# 4 peaks matched the calibration list.
cres <- subset(cres,nrmatch>4)
cres
#the error model of this "calibintstat" objects are
#applied to correct the error of the massvectors.
mvlfC <- applyintcalib(mvlfB,cres)

plot(mvl,xlim=c(1712.7,1713.1))
plot(mvlfA,add=T,col="green")
plot(mvlfB,add=T,col="blue")
plot(mvlfC,add=T,col="darkorange")


# Now when all calibration steps are performed
# its time to estimxate the error.
# Either Plot the histogram for a mass range
# where in all samples a peak occurs
hist(mvlfC,xlim=c(1712.7,1713.1),accur=0.05)
# or unlist the massvectorlist
# unlist makes one massvector out of a massvectorlist
tt<-unlist(mvlfC)
m<-tt[,1]
nt<-m[m>1712.7 & m<1713.4]
hist(nt,freq=F,ylim=c(0,35),xlim=c(1712.8,1713))
lines(density(nt),col=2)
#the error are +- 0.05 Da
#> 0.05/1700*1e6
#29.41176
# 30 ppm

#get abundant masses that can occure
#in at least than 30 massvectors
res<-gamasses(mvlfC,abund=25)
res
plot(res)
res2 <- mvFilter(mvlfC,res,abundant=TRUE)
image(res2,what="lengthmv")
image(res2,what="area.Max.")

##Now remove abundant masses!
mvlff <- mvFilter(mvlfC,res,abundant=FALSE)
par(mfrow=c(1,2))
hist(mvlfC)
hist(mvlff)
par(mfrow=c(1,1))

#Compute all intra massvector mass differences!!
tmp<-getdiff(mvlff,range=c(0,100))
gtmp<-gamasses(tmp,abund=50)
plot(gtmp)

