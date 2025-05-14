                                        #fit=glm(lucrative~log(prop0.15/prop65.more)+log(prop16.24/prop65.more)+log(prop25.64/prop65.more)+evolution+median+decile.ratio+activity+poverty+prop.public+prop.1to9+prop.industry+prop.trade,family="poisson",offset=log(pop.level),data=stddata)
fitluc=glm(lucrative~prop0.15+prop16.24+prop25.64+evolution+median+decile.ratio+activity+poverty+prop.public+prop.1to9+prop.industry+prop.trade,family="poisson",offset=log(pop.level),data=stddata)

fitcoop=glm(cooperative~prop0.15+prop16.24+prop25.64+evolution+median+decile.ratio+activity+poverty+prop.public+prop.1to9+prop.industry+prop.trade,family="poisson",offset=log(pop.level),data=stddata)

coef(fitluc)

 (Intercept)      prop0.15     prop16.24     prop25.64     evolution
 -1.126029485  -0.175648330   0.140224509  -0.170331902  -0.019646935
       median  decile.ratio      activity       poverty   prop.public
  0.089931663   0.078502114   0.002244272   0.035138214   0.216196920
    prop.1to9 prop.industry    prop.trade
0.061920506   0.327446003   0.470239715

as.numeric(-round(coef(fitluc)-coef(fitcoop),2))

 [1]  1.74 -0.47  0.00  0.34  0.04  0.05 -0.29  0.24  0.22 -0.22 -0.11 -0.33
[13] -0.58
