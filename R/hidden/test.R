#testing massvector
rm(list=ls())
source("all.R")
#constructors for massvector

massvector(NULL,"hello march")
massvector(1:10,"hello march")
massvector(cbind(1:10,10:1),"hello march")
tmp<-cbind(1:10,1:10)
colnames(tmp)<-c("mass","test")
rr<-massvector(tmp,"hello march")
rr2<-massvector(cbind(1:12,1:12),"hello bart")


# plot functions for massvector
plot(rr)
hist(rr)
summary(rr)
image(rr)
info(rr)

# setting new masses
mass(rr)
mass(rr,1:10)
peaks(rr,cbind(1:10,11:20))


#testing with real peaklists
mv1<-readBruker(massvector(),"D:/DATA/020529r_dt_Pirellula/0_A1_1SRef/pdata/1/peaklist.xml")
mv2<-readBruker(massvector(),"D:/DATA/020529r_dt_Pirellula/0_B10_1SRef/pdata/1/peaklist.xml")
# plotting with masses
plot(mv1,mv2)
image(mv1  ,mv2)
summary(mv1)
print(mv1)
hist(mv1)
plot(mv1)
image(mv1,error=199,ppm=F)

#testing with massvectorlist
mvl <- readBruker(massvectorlist(),"D:/DATA/020529r_dt_Pirellula/")
plot(mvl)
summary(mvl)
hist(mvl)
hist(mvl,length=F)
image(mvl) #fehler.

mvl2<-mvl[1:100]
plot(mvl2)
summary(mvl2)
hist(mvl2)
image(mvl2,what="lengthmv")


##testing assingments
mvl2[[11]]<-mvl2[[1]]
plot(mvl2[[11]],mvl2[[1]])
image(mvl2[[11]],mvl2[[1]])

##make one massvector out of the peaklist
tt<-unlist(mvl)
plot(tt)


##Wool and smilanski filtering
length(WsFilter(mvl[[3]]))
length(mvl[[3]])


res <- WsFilter(mvl,peptides = TRUE)
plot(res)
res2 <- WsFilter(mvl,peptides = FALSE)
plot(res2,col=2,add=TRUE)
image(res2)

res <- getrecalib(mv1)
print(res)
as.vector(res)
summary(res)
image(res)
plot(res)



#Filtering for abundant masses.
res<-gamasses(mvl,abund=50)
plot(res)
mvFilter(mvl[[1]],res)
res2<-mvFilter(mvl,res,abundant=TRUE)
image(res2,what="lengthmv")
X11()
image(mvl,what="lengthmv")
image(image(res2,what="lengthmv")/image(mvl,what="lengthmv"))

hist(mvl,length=F,accur=0.3)
hist(res2,length=F,add=TRUE,col=2,accur=0.3)






res <- getrecalib(mvlf)
print(res)
summary(res)
image(res)
#image(res,what="Coef.Intercept")
#image(res,what="Coef.Slope")
plot(res)
hist(res)
#tmp <- as.matrix(res)
#tmp2<- as.data.frame(res)
dres<-as.data.frame(res)
plot(dres$Coef.Intercept,dres$PQM,xlab="Coef.Intercept",ylab="PQM")

res2<-subset(res,PQM>10)
length(res2)
plot(res2)
test<-applyrecalib(mvl, res2)
mvl<-getcalib(test)
plot(tmp)

                                        #plot(unlist(lapply(mvl2,length)),unlist(lapply(mvl,length))

pp<-recalibrate(mvl)
stat<-getcalib(pp)
plot(stat)
hist(stat)
summary(stat)
image(stat,what="Coef.Intercept")

res<-getintcalib(mv1,calib)
plot(res)
hist(res)
print(res)
summary(res)
as.vector(res)


res<-getintcalib(mvl,calib)
plot(res)
hist(res)
print(res)
summary(res)
image(res)
tmp<-as.matrix(res)

res2<-subset(res)
res2<-subset(res,error.stdv < 80)
res3<-subset(res,error.stdv > 80)

plot(res2)
hist(res2)
print(res2)
summary(res2)

tt<- applyintcalib(mvl,res)

res<-getintcalib(mvl[names(res3)],calib,error=200)
res<- applyintcalib(mvl,res)

stat<-getcalib(res)
plot(stat)

amv<-correctinternal(mv1,calib)
amvl1 <- correctinternal(mvl,calib)
compare(amvl1[[1]],mvl[[1]])


get(amvl1,"calib")
getcalib(amvl1[[1]])


stat<-getcalib(amvl1)
plot(stat)
hist(stat)
summary(stat)
image(stat)

plot(stat[[1]])
hist(stat[[1]])
summary(stat[[1]])
image(stat[[1]])


#testing calib spline
ppg<-readBruker(massvectorlist(),"D:/DATA/020529r_dt_Pirellula_ppg/")
tmp <- getCalibSpline(ppg[[1]])
summary(tmp)
plot(tmp)
hist(tmp)
image(tmp)
print(tmp)

tmp <- getCalibSpline(ppg,error=200)
plot(tmp)
hist(tmp)
image(tmp)
print(tmp)
res <- applySpline(mv1,tmp)
compare(res,mv1,plot=T)

res <- applySpline(mvl,tmp)
compare(res,mvl)



res<-calibspline(pp[[1]],ppg)
res<-calibspline(mvl,ppg)
stat<-getcalib(res)
print(stat)


res<-cor(mvl,abund=50)
getintcalib(res,calib)
tt <- globalCalib(mvl,calib,abund=0)

res<- getglobalCalib(mvl,calib,abund=50)
plot(res)
hist(res)



res<- applyintcalib(mvl,res)
stat<-getcalib(res)
plot(stat)





tmp<-seq(1000,1003,1)*1.000423
plot(tmp,rep(1,length(tmp)),type="h",ylim=c(0,2),xlab="m/z",ylab="")
abline(h=0)
#curve(cos(x*2*pi/1.000423),add=T)
curve(cos(x*2*pi/1.000423)^2,add=T)
curve(sin(x*2*pi/1.000423)^2,add=T,col=1,lty=2)

curve(cos(x*2*pi/1.000495)^2,add=T,col=2)
curve(sin(x*2*pi/1.000495)^2,add=T,col=2,lty=2)
legend(1001,2,c("1000423","1000495"),col=c(1,2),lty=c(1,1))
sum(cos(tmp*2*pi/1.000423))
sum(cos(tmp*2*pi/1.000495))

