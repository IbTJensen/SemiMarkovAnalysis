load('../Data/pointsAndPoly.RData')

library(ggplot2)
library(glmnet)
library(reshape2)

###### a graph example ######
##############################
pl=function(v,trans='log10',size=5){
  p=ggplot(poly,aes(x=x,y=y,group=code,fill=.data[[v]]))+
    geom_polygon()+xlab('')+ylab('')+theme_minimal()+
    theme(legend.position=c(.18,.25),legend.text=element_text(size=s))+
    labs(fill='')+ggtitle(v)+
    guides(fill=guide_colourbar(barheight=4))+
    scale_fill_distiller(palette="Spectral",trans=trans)

  p
}
s=7
p1=pl('median',size=s)
p1
p2=pl('decile.ratio',size=s)
p2
p3=pl('activity',size=s)
p4=pl('poverty',size=s)
p5=pl('evolution',trans='identity',size=s)
p6=pl('prop0.15',size=s)
p7=pl('prop16.24',size=s)
p8=pl('prop25.64',size=s)
p9=pl('prop65.more',size=s)
p10=pl('prop.public',size=s)
p11=pl('prop.1to9',size=s)
p12=pl('prop.industry',size=s)
p13=pl('prop.trade',size=s)
p13
install.packages("ggpubr")
require(ggpubr)
pI=ggarrange(p1,p2,p3,p4);
h=6
pI
ggsave(pI,file='~/Dropbox/CooperativeBanks/R/covarI.pdf',height=h,width=h)
pII=ggarrange(p6,p7,p8,p5);
pII
ggsave(pII,file='~/Dropbox/CooperativeBanks/R/covarII.pdf',height=h,width=h)
pIII=ggarrange(p11,p10,p12,p13);
pIII
ggsave(pIII,file='~/Dropbox/CooperativeBanks/R/covarIII.pdf',height=h,width=h)
