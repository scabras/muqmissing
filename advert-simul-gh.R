rm(list = ls())

#We use a smaller data set (100)
advert=read.table("advert-subsample.csv",sep=" ",header=TRUE,dec=".")
n=nrow(advert)
#################################################################
#Proceso de simulacion, lo saltamos
 lm.expl.radio<-lm(advert$radio~advert$sales+advert$newspaper)
 summary(lm.expl.radio)

 new.radio<-rnorm(n,mean=lm.expl.radio$coefficients[1]+lm.expl.radio$coefficients[2]*advert$sales+lm.expl.radio$coefficients[2]/6*advert$newspaper,sd=8)
 TV=scale(advert$TV)
 new.radio<-scale(new.radio)
 newspaper=scale(advert$newspaper)
 X1=as.matrix(cbind(TV,new.radio))
 alpha=3
 beta=c(2,1.5)
 sigma=3
 y=alpha+X1%*%beta+rnorm(n,0,sigma*sqrt(abs(TV)))

 #we make y positive:
 y<-y-(min(y)-1)

 new.advert=data.frame(cbind(X1,newspaper,y))
 names(new.advert)<-names(advert)

 