#####PCoA analysis 202201

library(ecodist)
library(vegan)
library(ggplot2)

##### put soil and root samples raw sequence data together ####
ngs.all<-rbind(ngs.root,ngs.soil)
env.all<-rbind(env.root,env.soil)

ngs.log<-log2(ngs.all+1)

ngs.bray.all<-vegdist(ngs.log, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the dataset native to vegan
pcoaVS.all<-pco(ngs.bray.all,negvals="zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999

plot(pcoaVS.all$vectors[,2]~pcoaVS.all$vectors[,1],xlab = "PCoA1", ylab="PCoA2",pch=19,col=as.factor(env.all$type),
     axes = TRUE, main = "PCoA (ecodist) on Soil and Root")

pcoa.all<-cbind(pcoaVS.all$vectors[,2],pcoaVS.all$vectors[,1],env.all)


#### PCoA for root samples ####
ngs.log1<-log2(ngs.root+1)

ngs.bray.root<-vegdist(ngs.log1, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the dataset native to vegan
pcoaVS.root<-pco(ngs.bray.root,negvals="zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999

plot(pcoaVS.root$vectors[,2]~pcoaVS.root$vectors[,1],xlab = "PCoA1", ylab="PCoA2",pch=19,col=env.root$field,
     axes = TRUE, main = "PCoA (ecodist) on Roots")

pcoa.root<-cbind(pcoaVS.root$vectors[,2],pcoaVS.root$vectors[,1],env.root)

## code from Marine to add ellipse
#p1=p + stat_ellipse(geom = "polygon", type="norm",  level = 0.95,alpha=0.2, aes(fill=Pratique))

#### PCoA for soil samples ####
ngs.log2<-log2(ngs.soil+1)

ngs.bray.soil<-vegdist(ngs.log2, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the dataset native to vegan
pcoaVS.soil<-pco(ngs.bray.soil,negvals="zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999

plot(pcoaVS.soil$vectors[,2]~pcoaVS.soil$vectors[,1],xlab = "PCoA1", ylab="PCoA2",pch=19,col=env.soil$field,
     axes = TRUE, main = "PCoA (ecodist) on Soil")

pcoa.soil<-cbind(pcoaVS.soil$vectors[,2],pcoaVS.soil$vectors[,1],env.soil)

pdf("PCoA of root and soil samples in field raw data.pdf",height=5,width=7)
ggplot(data=pcoa.all,aes(x=pcoaVS.all$vectors[, 1],y=pcoaVS.all$vectors[, 2]))+
  xlim(-0.08,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(type)),size=3)+ 
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=type))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on root and soil samples")
dev.off()

pdf("PCoA of root and soil samples in field 1.pdf",height=5,width=7)
ggplot(data=pcoa.all,aes(x=pcoaVS.all$vectors[, 1],y=pcoaVS.all$vectors[, 2]))+
  xlim(-0.1,0.1)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(field)),size=3)+ 
  #stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=field))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on root and soil samples")
dev.off()


pdf("PCoA of root samples raw data.pdf",height=4,width=7)
ggplot()+
  geom_point(data=pcoa.root,aes(x=pcoaVS.root$vectors[, 1],y=pcoaVS.root$vectors[, 2],
                                 color=as.factor(field)),position=position_dodge(width=0.75),size=3)+
  xlim(-0.25,0.25)+
     #stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=weed.neighbour))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on Wheat roots")
dev.off()


pdf("PCoA of soil samples raw data.pdf",height=4,width=7)
ggplot()+
  geom_point(data=pcoa.soil,aes(x=pcoaVS.soil$vectors[, 1],y=pcoaVS.soil$vectors[, 2],
                                 color=as.factor(field)),position=position_dodge(width=0.75),size=3)+
  xlim(-0.4,0.6)+
  #stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=weed.neighbour))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on Weed roots")
dev.off()

#### statistics to check the effect of weed species on root mycobiota community composition ####
adonis.all<-adonis(ngs.log~type,env.all,by=NULL, method="bray") 
adonis.all
rda(ngs.log)

adonis.root<-adonis(ngs.log1~field,env.root,by=NULL, method="bray") 
adonis.root

adonis.soil<-adonis(ngs.log2~field,env.soil,by=NULL, method="bray") 
adonis.soil

######################################################################
##### put soil and root samples rarefied sequence data together ####
ngs.all.r<-rbind(ngs.broot,ngs.bsoil)

ngs.log.r<-log2(ngs.all.r+1)

ngs.bray.all.r<-vegdist(ngs.log.r, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the dataset native to vegan
pcoaVS.all.r<-pco(ngs.bray.all.r,negvals="zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999

pcoa.all.r<-cbind(pcoaVS.all.r$vectors[,2],pcoaVS.all.r$vectors[,1],env.all)

#### PCoA for root samples ####
ngs.log1.r<-log2(ngs.broot+1)

ngs.bray.root.r<-vegdist(ngs.log1.r, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the dataset native to vegan
pcoaVS.root.r<-pco(ngs.bray.root.r,negvals="zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999

plot(pcoaVS.root.r$vectors[,2]~pcoaVS.root.r$vectors[,1],xlab = "PCoA1", ylab="PCoA2",pch=19,col=env.root$field,
     axes = TRUE, main = "PCoA (ecodist) on Roots after rarefaction")

pcoa.root.r<-cbind(pcoaVS.root.r$vectors[,2],pcoaVS.root.r$vectors[,1],env.root)

#### PCoA for soil samples ####
ngs.log2.r<-log2(ngs.bsoil+1)

ngs.bray.soil.r<-vegdist(ngs.log2.r, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the dataset native to vegan
pcoaVS.soil.r<-pco(ngs.bray.soil.r,negvals="zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999

plot(pcoaVS.soil.r$vectors[,2]~pcoaVS.soil.r$vectors[,1],xlab = "PCoA1", ylab="PCoA2",pch=19,col=env.soil$field,
     axes = TRUE, main = "PCoA (ecodist) on Soil after rarefaction")

pcoa.soil.r<-cbind(pcoaVS.soil.r$vectors[,2],pcoaVS.soil.r$vectors[,1],env.soil)

pdf("PCoA of root and soil samples in field after rarefaction.pdf",height=5,width=7)
ggplot(data=pcoa.all.r,aes(x=pcoaVS.all.r$vectors[, 1],y=pcoaVS.all.r$vectors[, 2]))+
  xlim(-0.08,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(type)),size=3)+ 
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=type))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on root and soil samples after rarefaction")
dev.off()

pdf("PCoA of root and soil samples in field after rarefaction 1.pdf",height=5,width=7)
ggplot(data=pcoa.all,aes(x=pcoaVS.all$vectors[, 1],y=pcoaVS.all$vectors[, 2]))+
  xlim(-0.055,0.035)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(field)),size=3)+ 
  #stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=field))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on root and soil samples")
dev.off()


pdf("PCoA of root samples after rarefaction.pdf",height=4,width=7)
ggplot()+
  geom_point(data=pcoa.root.r,aes(x=pcoaVS.root.r$vectors[, 1],y=pcoaVS.root.r$vectors[, 2],
                                color=as.factor(field)),position=position_dodge(width=0.75),size=3)+
  xlim(-0.3,0.3)+
  #stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=as.factor(field)))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on Wheat roots after rarefaction")
dev.off()


pdf("PCoA of soil samples after rarefaction.pdf",height=4,width=7)
ggplot()+
  geom_point(data=pcoa.soil.r,aes(x=pcoaVS.soil.r$vectors[, 1],y=pcoaVS.soil.r$vectors[, 2],
                                color=as.factor(field)),position=position_dodge(width=0.75),size=3)+
  xlim(-0.6,0.4)+
  #stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=weed.neighbour))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on soil after rarefaction")
dev.off()

#### statistics to check the effect of weed species on root mycobiota community composition ####

adonis.root<-adonis(ngs.log1.r~otu.richness.f+otu.richness.soil,richness.all1,by=NULL, method="bray") 
adonis.root

adonis.root<-adonis(ngs.log1.r~otu.richness.f,env.root,by=NULL, method="bray") 
adonis.root

adonis.soil<-adonis(ngs.log2.r~otu.richness.f,env.soil,by=NULL, method="bray") 
adonis.soil

adonis.root<-adonis(ngs.log1.r~otu.shannon.f,env.root,by=NULL, method="bray") 
adonis.root

adonis.soil<-adonis(ngs.log2.r~otu.shannon.f,env.soil,by=NULL, method="bray") 
adonis.soil

adonis.root<-adonis(ngs.log1.r~otu.evenness.f,env.root,by=NULL, method="bray") 
adonis.root

adonis.soil<-adonis(ngs.log2.r~otu.evenness.f,env.soil,by=NULL, method="bray") 
adonis.soil

adonis.root<-adonis(ngs.log1.r~Cirsium_arvense+Galium_aparine+Lamium_purpureum+Matricaria_sp.+Papaver_rhoeas+Poa_annua+Poa_trivialis+Trifolium_repens+Veronica_persica+Vicia_sativa,all.weed,by=NULL, method="bray") 
adonis.root

adonis.soil<-adonis(ngs.log2.r~Cirsium_arvense+Galium_aparine+Lamium_purpureum+Matricaria_sp.+Papaver_rhoeas+Poa_annua+Poa_trivialis+Trifolium_repens+Veronica_persica+Vicia_sativa,all.weed,by=NULL, method="bray") 
adonis.soil

all.weed<-cbind(richness.phylum,otu.richness.asco.s,otu.shannon.asco.s,otu.evenness.asco.s,
                otu.richness.basi.s,otu.shannon.basi.s,otu.evenness.basi.s,
                otu.richness.chyt.s,otu.shannon.chyt.s,otu.evenness.chyt.s,
                otu.richness.glom.s,otu.shannon.glom.s,otu.evenness.glom.s,
                otu.richness.zygo.s,otu.shannon.zygo.s,otu.evenness.zygo.s,floristic.t)

#### For soil mycobiota
lm.s<-glmer.nb(otu.richness.soil~Cirsium_arvense+Galium_aparine+Lamium_purpureum+Matricaria_sp.+Papaver_rhoeas+Poa_annua+Poa_trivialis+Trifolium_repens+Veronica_persica+Vicia_sativa+(1|field),all.weed)
Anova(lm.s)
AIC(lm.s)
r.squaredGLMM(lm.s)

lm.s<-lmer(otu.evenness.soil~Cirsium_arvense+Galium_aparine+Lamium_purpureum+Matricaria_sp.+Papaver_rhoeas+Poa_annua+Poa_trivialis+Trifolium_repens+Veronica_persica+Vicia_sativa+(1|field),all.weed)
Anova(lm.s)
AIC(lm.s)
r.squaredGLMM(lm.s)

###
lm.s<-glmer.nb(otu.richness.asco.s~Cirsium_arvense+Galium_aparine+Lamium_purpureum+Matricaria_sp.+Papaver_rhoeas+Poa_annua+Poa_trivialis+Trifolium_repens+Veronica_persica+Vicia_sativa+(1|field),all.weed)
Anova(lm.s)
AIC(lm.s)
r.squaredGLMM(lm.s)

lm.s<-lmer(otu.evenness.asco.s~Cirsium_arvense+Galium_aparine+Lamium_purpureum+Matricaria_sp.+Papaver_rhoeas+Poa_annua+Poa_trivialis+Trifolium_repens+Veronica_persica+Vicia_sativa+(1|field),all.weed)
Anova(lm.s)
AIC(lm.s)
r.squaredGLMM(lm.s)

#### For root mycobiota
lm.r<-glmer.nb(otu.richness.root~Cirsium_arvense+Galium_aparine+Lamium_purpureum+Matricaria_sp.+Papaver_rhoeas+Poa_annua+Poa_trivialis+Trifolium_repens+Veronica_persica+Vicia_sativa+(1|field),all.weed)
Anova(lm.r)
AIC(lm.r)
r.squaredGLMM(lm.r)

lm.r<-lmer(otu.evenness.root~Cirsium_arvense+Galium_aparine+Lamium_purpureum+Matricaria_sp.+Papaver_rhoeas+Poa_annua+Poa_trivialis+Trifolium_repens+Veronica_persica+Vicia_sativa+(1|field),all.weed)
Anova(lm.r)
AIC(lm.r)
r.squaredGLMM(lm.r)


