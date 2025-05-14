library(spatstat)

codes=unique(poly$code)
length(codes)

polygs=list()

for (i in 1:length(codes)){
    admunit=poly[poly$code==codes[i],]

    polygs[[i]]=owin(poly=list(x=admunit$x,y=admunit$y))
}

tt <- tess(tiles = polygs)

#ALR transformation as Jeff proposes
temp=data
attach(temp)
temp[,'prop0.15']=log(temp[,'prop0.15']/temp[,'prop65.more'])
temp[,'prop16.24']=log(temp[,'prop16.24']/temp[,'prop65.more'])
temp[,'prop25.64']=log(temp[,'prop25.64']/temp[,'prop65.more'])


#create list of images for all covariates
covariates=list()

relevant=c(3:5,7,c(10:17))
l=1
for (i in relevant){
    tmp=temp[,i]
    m=mean(na.omit(tmp))
    s=sd(na.omit(tmp))
    tmp=(tmp-m)/s
    covariates[[l]]=as.im.tess(tt,values=tmp,dimyx=c(256,256))
    l=l+1
}

save(covariates,file="covariates.rdata")
#load("covariates.rdata")


#create list of functions

covariatesfct=list()
l=1
for (i in relevant){
    tmp=temp[,i]
    m=mean(na.omit(tmp))
    s=sd(na.omit(tmp))
    tmp=(tmp-m)/s
    covariatesfct[[l]]=as.function(tt, values=tmp)
    l=l+1
}

names(covariatesfct)=names(data)[relevant]
save(covariatesfct,file="covariatesfct.rdata")
#load("covariatesfct.rdata")
#plot(covariatesfct[[1]])
#der er nogle "hvide" dele i funktionen.

#load("covariates.rdata")
#plot(covariates[[1]])
#det er der ogsaa i im.
#tesselationen ser OK ud.

#standardized data
stddata=data
l=1
for (i in relevant){
    tmp=temp[,i]
    m=mean(na.omit(tmp))
    s=sd(na.omit(tmp))
    tmp=(tmp-m)/s
    stddata[,i]=tmp
    l=l+1
}
save(stddata,file="stddata.rdata")
load
mydensity=stddata$pop.level/stddata$area#lig med density.


logdensity.fct=as.function(tt, values=log(data$density))
plot(logdensity.fct)

save(logdensity.fct,file="logdensity.rdata")
