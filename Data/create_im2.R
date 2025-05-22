library(spatstat)

load("Data/dataALR.RData")
load("Data/pointsAndPoly.RData")

codes=unique(data.alr$code)
length(codes)

polygs=list()

for (i in 1:length(codes)){
  cat(i, "\r")
  admunit=poly[poly$code==codes[i],]
  polygs[[i]]=owin(poly=list(x=admunit$x,y=admunit$y))
}

tt <- tess(tiles = polygs)

#create list of functions
relevant=c(1:2, 11, 13:24)

covariatesfct=list()
l=1
for (i in relevant){
    tmp=data.alr[,i]
    # m=mean(na.omit(tmp))
    # s=sd(na.omit(tmp))
    # tmp=(tmp-m)/s
    covariatesfct[[l]]=as.function(tt, values=tmp)
    l=l+1
}


names(covariatesfct)=c("code", "city", "density", "prop0.15", "prop16.24",
                       "prop25.64", "evolution", "poverty", "activity",
                       "median", "decile", "prop19", "proppublic",
                       "propindustry", "proptrade")
save(covariatesfct, file = "Data/covariatesfct_alr.rdata")
save(tt, file = "Data/regions.rdata")
#load("covariatesfct.rdata")
#plot(covariatesfct[[1]])
#der er nogle "hvide" dele i funktionen.

#load("covariates.rdata")
#plot(covariates[[1]])
#det er der ogsaa i im.
#tesselationen ser OK ud.

# par(mfrow=c(2,3))
# for (i in 1:6)
#     plot(covariatesfct[[i]])
# for (i in 7:12)
#     plot(covariatesfct[[i]])
