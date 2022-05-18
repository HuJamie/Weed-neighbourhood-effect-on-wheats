####taxonomic summary

ngs.asco<-subset(t(ngs.root),taxo$Phylum=="Ascomycota")
ngs.basi<-subset(t(ngs.root),taxo$Phylum=="Basidiomycota")
ngs.chyt<-subset(t(ngs.root),taxo$Phylum=="Chytridiomycota")
ngs.glom<-subset(t(ngs.root),taxo$Phylum=="Glomeromycota")
ngs.zygo<-subset(t(ngs.root),taxo$Phylum=="Zygomycota")

ngs.asco.s<-subset(t(ngs.soil),taxo$Phylum=="Ascomycota")
ngs.basi.s<-subset(t(ngs.soil),taxo$Phylum=="Basidiomycota")
ngs.chyt.s<-subset(t(ngs.soil),taxo$Phylum=="Chytridiomycota")
ngs.glom.s<-subset(t(ngs.soil),taxo$Phylum=="Glomeromycota")
ngs.zygo.s<-subset(t(ngs.soil),taxo$Phylum=="Zygomycota")

#####
##### calculate diversity index inside each phylum for wheat root ####
library(vegan)
otu.shannon.asco.r<- diversity(t(ngs.asco),index="shannon") 
otu.richness.asco.r<-rowSums(t(ngs.asco)>0)
otu.evenness.asco.r<- otu.shannon.asco.r/log(otu.richness.asco.r)  ## Pielou's evenness

########Basidiomycota
otu.shannon.basi.r<- diversity(t(ngs.basi),index="shannon") 
otu.richness.basi.r<-rowSums(t(ngs.basi)>0)
otu.evenness.basi.r<- otu.shannon.basi.r/log(otu.richness.basi.r)  ## Pielou's evenness

##### Chytridiomycota
otu.shannon.chyt.r<- diversity(t(ngs.chyt),index="shannon") 
otu.richness.chyt.r<-rowSums(t(ngs.chyt)>0)
otu.evenness.chyt.r<- otu.shannon.chyt.r/log(otu.richness.chyt.r)  ## Pielou's evenness

######Glomeromycota
otu.shannon.glom.r<- diversity(t(ngs.glom),index="shannon") 
otu.richness.glom.r<-rowSums(t(ngs.glom)>0)
otu.evenness.glom.r<- otu.shannon.glom.r/log(otu.richness.glom.r)  ## Pielou's evenness

######Zygomycota
otu.shannon.zygo.r<- diversity(t(ngs.zygo),index="shannon") 
otu.richness.zygo.r<-rowSums(t(ngs.zygo)>0)
otu.evenness.zygo.r<- otu.shannon.zygo.r/log(otu.richness.zygo.r)  ## Pielou's evenness


richness.phylum<-cbind(richness.all1,otu.richness.asco.r,otu.shannon.asco.r,otu.evenness.asco.r,
                            otu.richness.basi.r,otu.shannon.basi.r,otu.evenness.basi.r,
                            otu.richness.chyt.r,otu.shannon.chyt.r,otu.evenness.chyt.r,
                            otu.richness.glom.r,otu.shannon.glom.r,otu.evenness.glom.r,
                            otu.richness.zygo.r,otu.shannon.zygo.r,otu.evenness.zygo.r)

library(ggplot2)
library(gridExtra)
pdf("Richness in wheat root phylum v.s floristic richness.pdf",height=4,width=18)
a<-ggplot(data=richness.phylum,aes(x=otu.richness.f,y=otu.richness.asco.r))+ 
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster richness of wheat root Ascomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

b<-ggplot(data=richness.phylum,aes(x=otu.richness.f,y=otu.richness.basi.r))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster richness of wheat root Basidiomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

c<-ggplot(data=richness.phylum,aes(x=otu.richness.f,y=otu.richness.chyt.r))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster richness of wheat root Chytridiomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

g<-ggplot(data=richness.phylum,aes(x=otu.richness.f,y=otu.richness.glom.r))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster richness of wheat root Glomeromycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

z<-ggplot(data=richness.phylum,aes(x=otu.richness.f,y=otu.richness.zygo.r))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster richness of wheat root Zygomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(a,b,c,g,z,ncol=5)
dev.off()

pdf("Evenness in wheat root phylum v.s floristic richness.pdf",height=4,width=18)
ae<-ggplot(data=richness.phylum,aes(x=otu.richness.f,y=otu.evenness.asco.r))+ 
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster evenness of wheat root Ascomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

be<-ggplot(data=richness.phylum,aes(x=otu.richness.f,y=otu.evenness.basi.r))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster evenness of wheat root Basidiomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

ce<-ggplot(data=richness.phylum,aes(x=otu.richness.f,y=otu.evenness.chyt.r))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster evenness of wheat root Chytridiomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

ge<-ggplot(data=richness.phylum,aes(x=otu.richness.f,y=otu.evenness.glom.r))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster evenness of wheat root Glomeromycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

ze<-ggplot(data=richness.phylum,aes(x=otu.richness.f,y=otu.evenness.zygo.r))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster evenness of wheat root Zygomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(ae,be,ce,ge,ze,ncol=5)
dev.off()

### Statistics for diversity inside each phylum ####
lm.r<-glmer.nb(otu.richness.root~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.r)
AIC(lm.r)
r.squaredGLMM(lm.r)

lm.r<-lmer(otu.evenness.root~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.r)
AIC(lm.r)
r.squaredGLMM(lm.r)

### effect for Ascomycota####
#lm.asc<-lmer(otu.richness.asco~log2(otu.richness.f)*log2(otu.richness.soil)+(1|field)+(1|quadrats),richness.phylum)
lm.asc<-glmer.nb(otu.richness.asco.r~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.evenness.asco.r~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

#### effect for Basidiomycota ####
lm.bas<-glmer.nb(otu.richness.basi.r~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
#lm.bas<-glmer.nb(otu.richness.basi~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lmer(otu.evenness.basi.r~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

#### effect for Chytridiomycota####
lm.chy<-glmer.nb(otu.richness.chyt.r~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
#lm.chy<-glmer.nb(otu.richness.chyt~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.evenness.chyt.r~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

#### effect for Glomeromycota ####
lm.glom<-glmer.nb(otu.richness.glom.r~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
#lm.glom<-glmer.nb(otu.richness.glom~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.glom)
AIC(lm.glom)
r.squaredGLMM(lm.glom)

lm.glom<-lmer(otu.evenness.glom.r~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.glom)
AIC(lm.glom)
r.squaredGLMM(lm.glom)

####effect on Zygomycota ####
lm.zygo<-glmer.nb(otu.richness.zygo.r~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.zygo)
AIC(lm.zygo)
r.squaredGLMM(lm.zygo)

lm.zygo<-lmer(otu.evenness.zygo.r~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.zygo)
AIC(lm.zygo)
r.squaredGLMM(lm.zygo)

############ calculate diversity index inside each phylum for soil mycobiota ####
library(vegan)
otu.shannon.asco.s<- diversity(t(ngs.asco.s),index="shannon") 
otu.richness.asco.s<-rowSums(t(ngs.asco.s)>0)
otu.evenness.asco.s<- otu.shannon.asco.s/log(otu.richness.asco.s)  ## Pielou's evenness

########Basidiomycota
otu.shannon.basi.s<- diversity(t(ngs.basi.s),index="shannon") 
otu.richness.basi.s<-rowSums(t(ngs.basi.s)>0)
otu.evenness.basi.s<- otu.shannon.basi.s/log(otu.richness.basi.s)  ## Pielou's evenness

##### Chytridiomycota
otu.shannon.chyt.s<- diversity(t(ngs.chyt.s),index="shannon") 
otu.richness.chyt.s<-rowSums(t(ngs.chyt.s)>0)
otu.evenness.chyt.s<- otu.shannon.chyt.s/log(otu.richness.chyt.s)  ## Pielou's evenness

######Glomeromycota
otu.shannon.glom.s<- diversity(t(ngs.glom.s),index="shannon") 
otu.richness.glom.s<-rowSums(t(ngs.glom.s)>0)
otu.evenness.glom.s<- otu.shannon.glom.s/log(otu.richness.glom.s)  ## Pielou's evenness

######Zygomycota
otu.shannon.zygo.s<- diversity(t(ngs.zygo.s),index="shannon") 
otu.richness.zygo.s<-rowSums(t(ngs.zygo.s)>0)
otu.evenness.zygo.s<- otu.shannon.zygo.s/log(otu.richness.zygo.s)  ## Pielou's evenness


richness.phylum.s<-cbind(richness.all1,otu.richness.asco.s,otu.shannon.asco.s,otu.evenness.asco.s,
                       otu.richness.basi.s,otu.shannon.basi.s,otu.evenness.basi.s,
                       otu.richness.chyt.s,otu.shannon.chyt.s,otu.evenness.chyt.s,
                       otu.richness.glom.s,otu.shannon.glom.s,otu.evenness.glom.s,
                       otu.richness.zygo.s,otu.shannon.zygo.s,otu.evenness.zygo.s)

library(ggplot2)
library(gridExtra)
pdf("Richness in soil phylum v.s floristic richness.pdf",height=4,width=18)
as<-ggplot(data=richness.phylum.s,aes(x=otu.richness.f,y=otu.richness.asco.s))+ 
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster richness of soil Ascomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

bs<-ggplot(data=richness.phylum.s,aes(x=otu.richness.f,y=otu.richness.basi.s))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster richness of soil Basidiomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

cs<-ggplot(data=richness.phylum.s,aes(x=otu.richness.f,y=otu.richness.chyt.s))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster richness of soil Chytridiomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

gs<-ggplot(data=richness.phylum.s,aes(x=otu.richness.f,y=otu.richness.glom.s))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster richness of soil Glomeromycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

zs<-ggplot(data=richness.phylum.s,aes(x=otu.richness.f,y=otu.richness.zygo.s))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster richness of soil Zygomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(as,bs,cs,gs,zs,ncol=5)
dev.off()

pdf("Evenness in phylum v.s floristic richness.pdf",height=4,width=18)
ae<-ggplot(data=richness.phylum.s,aes(x=otu.richness.f,y=otu.evenness.asco.s))+ 
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster evenness of wheat root Ascomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

be<-ggplot(data=richness.phylum.s,aes(x=otu.richness.f,y=otu.evenness.basi.s))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster evenness of wheat root Basidiomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

ce<-ggplot(data=richness.phylum.s,aes(x=otu.richness.f,y=otu.evenness.chyt.s))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster evenness of wheat root Chytridiomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

ge<-ggplot(data=richness.phylum.s,aes(x=otu.richness.f,y=otu.evenness.glom.s))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster evenness of wheat root Glomeromycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

ze<-ggplot(data=richness.phylum.s,aes(x=otu.richness.f,y=otu.evenness.zygo.s))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  labs(x='Weed species richness',y='Sequence cluster evenness of wheat root Zygomycota')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(ae,be,ce,ge,ze,ncol=5)
dev.off()

### Statistics for diversity inside each phylum for soil####
lm.s<-glmer.nb(otu.richness.soil~otu.richness.f+(1|field)+(1|quadrats),richness.phylum.s)
Anova(lm.s)
AIC(lm.s)
r.squaredGLMM(lm.s)

lm.s<-lmer(otu.evenness.soil~otu.richness.f+(1|field)+(1|quadrats),richness.phylum.s)
Anova(lm.s)
AIC(lm.s)
r.squaredGLMM(lm.s)

### effect for Ascomycota####
#lm.asc<-lmer(otu.richness.asco~log2(otu.richness.f)*log2(otu.richness.soil)+(1|field)+(1|quadrats),richness.phylum)
lm.asc.s<-glmer.nb(otu.richness.asco.s~otu.richness.f+(1|field)+(1|quadrats),richness.phylum.s)
Anova(lm.asc.s)
AIC(lm.asc.s)
r.squaredGLMM(lm.asc.s)

lm.asc.s<-lmer(otu.evenness.asco.s~otu.richness.f+(1|field)+(1|quadrats),richness.phylum.s)
Anova(lm.asc.s)
AIC(lm.asc.s)
r.squaredGLMM(lm.asc.s)

#### effect for Basidiomycota ####
lm.bas.s<-glmer.nb(otu.richness.basi.s~otu.richness.f+(1|field)+(1|quadrats),richness.phylum.s)
#lm.bas<-glmer.nb(otu.richness.basi~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.bas.s)
AIC(lm.bas.s)
r.squaredGLMM(lm.bas.s)

lm.bas.s<-lmer(otu.evenness.basi.s~otu.richness.f+(1|field)+(1|quadrats),richness.phylum.s)
Anova(lm.bas.s)
AIC(lm.bas.s)
r.squaredGLMM(lm.bas.s)

#### effect for Chytridiomycota####
lm.chy.s<-glmer.nb(otu.richness.chyt.s~otu.richness.f+(1|field)+(1|quadrats),richness.phylum.s)
#lm.chy<-glmer.nb(otu.richness.chyt~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.chy.s)
AIC(lm.chy.s)
r.squaredGLMM(lm.chy.s)

lm.chy.s<-lmer(otu.evenness.chyt.s~otu.richness.f+(1|field)+(1|quadrats),richness.phylum.s)
Anova(lm.chy.s)
AIC(lm.chy.s)
r.squaredGLMM(lm.chy.s)

#### effect for Glomeromycota ####
lm.glom.s<-glmer.nb(otu.richness.glom.s~otu.richness.f+(1|field)+(1|quadrats),richness.phylum.s)
#lm.glom<-glmer.nb(otu.richness.glom~otu.richness.f+(1|field)+(1|quadrats),richness.phylum)
Anova(lm.glom.s)
AIC(lm.glom.s)
r.squaredGLMM(lm.glom.s)

lm.glom.s<-lmer(otu.evenness.glom.s~otu.richness.f+(1|field)+(1|quadrats),richness.phylum.s)
Anova(lm.glom.s)
AIC(lm.glom.s)
r.squaredGLMM(lm.glom.s)

####effect on Zygomycota ####
lm.zygo.s<-glmer.nb(otu.richness.zygo.s~otu.richness.f+(1|field)+(1|quadrats),richness.phylum.s)
Anova(lm.zygo.s)
AIC(lm.zygo.s)
r.squaredGLMM(lm.zygo.s)

lm.zygo.s<-lmer(otu.evenness.zygo.s~otu.richness.f+(1|field)+(1|quadrats),richness.phylum.s)
Anova(lm.zygo.s)
AIC(lm.zygo.s)
r.squaredGLMM(lm.zygo.s)

