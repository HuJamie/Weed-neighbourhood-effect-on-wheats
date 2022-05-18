
#####Statistics from mixed-linear effect models #####
#####impact of weed diversity on soil and wheat microbiota #####
#### loading the packages 
library(nlme)
library(lme4)
library(MuMIn)
library(car)
library(lmerTest)
library(lsmeans)
library(MASS)
library(rsq)

############################################################################
####put floristic data and soil/root sequence data together ####

richness.r1<-richness.root[(row.names(richness.root) %in% row.names(richness.soil)), ]
richness.s1<-richness.soil[(row.names(richness.soil) %in% row.names(richness.root)), ]

richness.f1<-richness.f[(row.names(richness.f) %in% row.names(richness.s1)), ]

row.names(richness.s1) %in% row.names(richness.r1)
row.names(richness.f1) %in% row.names(richness.r1)

richness.all1<-cbind(richness.f,richness.soil,richness.root)

richness.f2<-richness.f[(row.names(richness.f) %in% row.names(richness.root)), ]
richness.fr<-cbind(richness.f2,richness.root)
floristic.root<-floristic.t[(row.names(richness.f) %in% row.names(richness.root)), ]

richness.f3<-richness.f[(row.names(richness.f) %in% row.names(richness.soil)), ]
richness.fs<-cbind(richness.f3,richness.soil)
floristic.soil<-floristic.t[(row.names(richness.f) %in% row.names(richness.soil)), ]

write.table(richness.fr,file="richness.fr.txt",sep="\t")
write.table(richness.fs,file="richness.fs.txt",sep="\t")
write.table(richness.all,file="richness.all.txt",sep="\t")

richness.fr1<-read.table("richness.fr1.txt", header=T, row.names=1, sep="\t")  
richness.fs1<-read.table("richness.fs1.txt", header=T, row.names=1, sep="\t")  

richness.all1<-read.table("richness.all.i.txt", header=T, row.names=1, sep="\t")

##### model 1 ####
mod<-glmer.nb(otu.richness.root~otu.richness.f+(1|field), data=richness.all1)
mod<-lmer(otu.richness.root~otu.richness.f+otu.shannon.f+otu.evenness.f+(1|field/Champs), data=richness.fr1)
mod<-lmer(otu.richness.root~otu.richness.soil+(1|field/Champs), data=richness.all1)
mod<-lmer(otu.richness.soil~otu.richness.f+otu.shannon.f+otu.evenness.f+(1|field/Champs), data=richness.fs1)

mod<-lmer(otu.shannon.root~otu.richness.f+otu.shannon.f+otu.evenness.f+otu.richness.soil+otu.shannon.soil+otu.evenness.soil+(1|field/Champs), data=richness.all1)
mod<-lmer(otu.shannon.root~otu.richness.f+otu.shannon.f+otu.evenness.f+(1|field/Champs), data=richness.fr1)
mod<-lmer(otu.shannon.root~otu.shannon.soil+(1|field/Champs), data=richness.all1)
mod<-lmer(otu.shannon.soil~otu.richness.f+otu.shannon.f+otu.evenness.f+(1|field/Champs), data=richness.fs1)

mod<-lmer(otu.evenness.root~otu.richness.f+otu.shannon.f+otu.evenness.f+otu.richness.soil+otu.shannon.soil+otu.evenness.soil+(1|field/Champs), data=richness.all1)
mod<-lmer(otu.evenness.root~otu.richness.f+otu.shannon.f+otu.evenness.f+(1|field/Champs), data=richness.fr1)
mod<-lmer(otu.evenness.root~otu.evenness.soil+(1|field/Champs), data=richness.all1)
mod<-lmer(otu.evenness.soil~otu.richness.f+otu.shannon.f+otu.evenness.f+(1|field/Champs), data=richness.fs1)

summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

####group comparison after model corrections
lsmeans(mod,pairwise~otu.richness.f,by="otu.richness.root",data=richness.all1)
lsmeans(mod,pairwise~otu.richness.root,by="otu.richness.f",data=richness.all1)

#### Try to draw some figures to show trend ####
xyplot(otu.richness.root~otu.richness.f,richness.all1)
xyplot(otu.shannon.root~otu.shannon.f,richness.all1)
xyplot(otu.evenness.root~otu.evenness.f,richness.all1)

xyplot(otu.richness.root~otu.richness.soil,richness.all1)
xyplot(otu.shannon.root~otu.shannon.soil,richness.all1)
xyplot(otu.evenness.root~otu.evenness.soil,richness.all1)

xyplot(otu.richness.soil~otu.richness.f,richness.all1)
xyplot(otu.shannon.soil~otu.shannon.f,richness.all1)
xyplot(otu.evenness.soil~otu.evenness.f,richness.all1)

####### COA random analysis ####
mod<-lmer(pcoaVS.root.r$vectors[,1]~otu.shannon.f+(1|field), data=pcoa.root.r)

summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

####group comparison after model corrections
lsmeans(mod,pairwise~otu.shannon.f,by="otu.shannon.soil",data=pcoa.root.r)
lsmeans(mod,pairwise~otu.shannon.soil,by="otu.shannon.f",data=richness.all1)

pcoaVS.root.r$vectors[,2]

pcoaVS.soil.r$vectors[,2]


