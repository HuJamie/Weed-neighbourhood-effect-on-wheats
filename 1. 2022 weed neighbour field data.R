##############################################
###### Data format, 202201, Jie Hu ####

ngs.soil<-read.table("ngs.soil.txt", header=T, row.names=1, sep="\t")  
ngs.root<-read.table("ngs.root.txt", header=T, row.names=1, sep="\t")  
floristic<-read.table("floristic.txt", header=T, row.names=1, sep="\t")

env.soil<-read.table("soil.env.txt", sep="\t", header=T,row.names=1)        ## sample properties, we need this
env.root<-read.table("root.env.txt", sep="\t", header=T,row.names=1)        ## sample properties, we need this
env.floristic<-read.table("floristic.env.txt", sep="\t", header=T,row.names=1)        ## sample properties, we need this

taxo<-read.table("taxa.txt", sep="\t",header=T,row.names=1)        ## sample properties, we need this

ngs.soil<-ngs.soil[(row.names(ngs.soil) %in% row.names(richness.all1)), ]
ngs.root<-ngs.root[(row.names(ngs.root) %in% row.names(richness.all1)), ]
floristic<-floristic[,(colnames(floristic) %in% row.names(richness.all1))]

env.soil<-env.soil[(row.names(env.soil) %in% row.names(richness.all1)), ]
env.root<-env.root[(row.names(env.root) %in% row.names(richness.all1)), ]
env.floristic<-env.floristic[(row.names(env.floristic) %in% row.names(richness.all1)), ]

write.table(floristic,file="floristic.i.txt",sep="\t")
write.table(ngs.soil,file="ngs.soil.i.txt",sep="\t")
write.table(ngs.root,file="ngs.root.i.txt",sep="\t")

write.table(env.floristic,file="env.floristic.i.txt",sep="\t")
write.table(env.soil,file="env.soil.i.txt",sep="\t")
write.table(env.root,file="env.root.i.txt",sep="\t")

######Reading the dataset ####
rm(list = ls()) #remove the data in the environment
ngs.soil<-read.table("ngs.soil.ir.txt", header=T, row.names=1, sep="\t")  
ngs.root<-read.table("ngs.root.ir.txt", header=T, row.names=1, sep="\t")  
floristic<-read.table("floristic.i.txt", header=T, row.names=1, sep="\t")

env.soil<-read.table("env.soil.i.txt", sep="\t", header=T,row.names=1)        ## sample properties, we need this
env.root<-read.table("env.root.i.txt", sep="\t", header=T,row.names=1)        ## sample properties, we need this
env.floristic<-read.table("env.floristic.i.txt", sep="\t", header=T,row.names=1)        ## sample properties, we need this

taxo<-read.table("taxa.txt", sep="\t",header=T,row.names=1)        ## sample properties, we need this

##### Microbial diversity of soil ####
#### overview, calculate value for method description ####
rowSums(ngs.soil)
mean(rowSums(ngs.soil))
min(rowSums(ngs.soil))
max(rowSums(ngs.soil))

#### normalization with minimum read
ngs.relat<-0*ngs.soil                                                     ## create a empty matrix
for(i in 1:dim(ngs.relat)[1]) ngs.relat[i,]<-ngs.soil[i,]/rowSums(ngs.soil)[i] ## loop to calculate relative abundance
rowSums(ngs.relat)                                                       ## check if calculate correctly

#### Standardization of abundance = rarefaction ####
ngs.bsoil<-round(ngs.relat*min(rowSums(ngs.soil)),0)   ## recover otu number based on relative otu abundance, rarefy
ngs.btsoil<-as.data.frame(cbind(t(ngs.bsoil),taxo))            ## combine otu and taxonomy
ngs.tsoil<-data.frame(t(ngs.btsoil))                     ## subset based on dimensionality of env.txt

mean(rowSums(ngs.bsoil))
min(rowSums(ngs.bsoil))
max(rowSums(ngs.bsoil))

library(vegan)
rarecurve(ngs.bsoil,step=100,xlim=c(0,22000),ylim=c(0,340),xlab="Sequence depth", ylab="Number of sequence clusters in soil mycobiota",col="blue",label=FALSE)

rarecurve(ngs.soil,step=100,xlim=c(0,32000),ylim=c(0,340),xlab="Sequence depth", ylab="Number of sequence clusters in soil mycobiota",col="blue",label=FALSE)

#write.table(ngs.btsoil,file="rarefied_matrix.txt",sep="\t")

otu.shannon.soil<- diversity(ngs.bsoil, index = "shannon") 
otu.richness.soil<-rowSums(ngs.bsoil>0)
otu.evenness.soil<- otu.shannon.soil/log(otu.richness.soil)        ## Pielou's evenness

richness.soil<-cbind(env.soil,otu.richness.soil,otu.shannon.soil,otu.evenness.soil)

##### microbial diversity of roots ####
#### overview, calculate value for method description ####
rowSums(ngs.root)
mean(rowSums(ngs.root))
min(rowSums(ngs.root))
max(rowSums(ngs.root))

#### normalization with minimum read
ngs.root1<-ngs.root
ngs.relat.root<-0*ngs.root1                                                     ## create a empty matrix
for(i in 1:dim(ngs.relat.root)[1]) ngs.relat.root[i,]<-ngs.root1[i,]/rowSums(ngs.root1)[i] ## loop to calculate relative abundance
rowSums(ngs.relat.root)                                                       ## check if calculate correctly

#### Standardization of abundance = rarefaction ####
ngs.broot<-round(ngs.relat.root*min(rowSums(ngs.root)),0)   ## recover otu number based on relative otu abundance, rarefy
ngs.btroot<-as.data.frame(cbind(t(ngs.broot),taxo))         ## combine otu and taxonomy
ngs.troot<-data.frame(t(ngs.btroot))                     ## subset based on dimensionality of env.txt

mean(rowSums(ngs.broot))
min(rowSums(ngs.broot))
max(rowSums(ngs.broot))

rarecurve(ngs.broot,step=100,xlim=c(0,14800),ylim=c(0,200),xlab="Sequence depth", ylab="Number of sequence clusters in wheat root mycobiota",col="blue",label=FALSE)

otu.shannon.root<- diversity(ngs.broot, index = "shannon") 
otu.richness.root<-rowSums(ngs.broot>0)
otu.evenness.root<- otu.shannon.root/log(otu.richness.root)        ## Pielou's evenness

richness.root<-cbind(env.root,otu.richness.root,otu.shannon.root,otu.evenness.root)

##### calaulation of floristic diversity ####
#### overview, calculate value for method description ####
colSums(floristic[,-1])

sum(colSums(floristic[,-1]))
test<-as.data.frame(rowSums(floristic[,-1])/sum(colSums(floristic[,-1])))

occurance<-floristic[,-1]/colSums(floristic[,-1])
write.table(occurance,file="occurance.txt",sep="\t")


mean(colSums(floristic[,-1]))
min(colSums(floristic[,-1]))
max(colSums(floristic[,-1]))
floristic.t<-as.data.frame(t(floristic[,-1]))

rarecurve(floristic.t,step=1,xlab="Sequence depth", ylab="Number of weed species in field",col="blue",label=FALSE)



otu.shannon.f<- diversity(floristic.t, index = "shannon") 
otu.richness.f<-rowSums(floristic.t>0)
otu.evenness.f<- otu.shannon.f/log(otu.richness.f)        ## Pielou's evenness

richness.f<-cbind(env.floristic,otu.richness.f,otu.shannon.f,otu.evenness.f)

richness.all1<-cbind(richness.f,richness.soil,richness.root)
#write.table(richness.all1,file="richness.all.i.txt",sep="\t")

#####To check which plant species is )####
otu.richness.ft<-as.data.frame(rowSums(floristic>0))
write.table(otu.richness.ft,file="otu.richness.ft.txt",sep="\t")

###### For the summary descriptions for sequence data summary####
ngs.psoil<-aggregate(ngs.btsoil[,1:60], list(taxo$Phylum), FUN = sum)    ## otu reads inside phylum,based on normalized matrix
rownames(ngs.psoil)<-ngs.psoil[,1]  
ngs.ptsoil<-t(ngs.psoil[,-1])

ngs.proot<-aggregate(ngs.btroot[,1:60], list(taxo$Phylum), FUN = sum)    ## otu reads inside phylum,based on normalized matrix
rownames(ngs.proot)<-ngs.proot[,1]  
ngs.ptroot<-t(ngs.proot[,-1])

sum(rowSums(ngs.psoil[,-1]))
sum(rowSums(ngs.proot[,-1]))

library(ggplot2)
library(gridExtra)

#### Soil 
ngs.psoil$seq<-rowSums(ngs.psoil[,-1])
ngs.psoil$prop<-rowSums(ngs.psoil[,-1])/sum(rowSums(ngs.psoil[,-1]))
ngs.psoil$number<-rowSums(ngs.psoil[,-1])
ngs.psoil$richness<-c(229,105,51,18,5,59)
ngs.psoil$phylum<-c("Ascomycota","Basidiomycota","Chytridiomycota","Glomeromycota","Multi-affiliation","GZygomycota")
ngs.psoil[order(ngs.psoil$phylum),]
ngs.psoil$variable<-c("Phyla","Phyla","Phyla","Phyla","Phyla","Phyla")

p.ps= ggplot(ngs.psoil, aes(x=variable,y=prop, fill = phylum )) + 
  geom_bar(stat = "identity",position="stack",width=0.8)+ 
  #scale_y_continuous(labels = scales::percent) + 
  xlab("Soil mycobiota in field")+
  ylab("Percentage of sequence clusters")+ theme_classic()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values =rev(c("Multi-affiliation"="#000000","GZygomycota"="#44A8DB",
                                  "Glomeromycota"="#8DD3C7","Chytridiomycota"="#57B78C",
                                  "Basidiomycota"="#BEBADA","Ascomycota"="#3C3A8D"))) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())
p.ps

p.rs= ggplot(ngs.psoil, aes(x=variable,y=richness, fill = phylum )) + 
  geom_bar(stat = "identity",position="stack",width=0.8)+ 
  #scale_y_continuous(labels = scales::percent) + 
  xlab("Soil mycobiota in field")+
  ylab("Richness of sequence clusters")+ theme_classic()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values =rev(c("Multi-affiliation"="#000000","GZygomycota"="#44A8DB",
                                  "Glomeromycota"="#8DD3C7","Chytridiomycota"="#57B78C",
                                  "Basidiomycota"="#BEBADA","Ascomycota"="#3C3A8D"))) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())
p.rs

#### wheat roots 
ngs.proot$seq<-rowSums(ngs.proot[,-1])
ngs.proot$prop<-rowSums(ngs.proot[,-1])/sum(rowSums(ngs.proot[,-1]))
ngs.proot$number<-rowSums(ngs.proot[,-1])
ngs.proot$phylum<-c("Ascomycota","Basidiomycota","Chytridiomycota","Glomeromycota","Multi-affiliation","GZygomycota")
ngs.proot$richness<-c(229,105,51,18,5,59)
ngs.psoil[order(ngs.proot$phylum),]
ngs.proot$variable<-c("Phyla","Phyla","Phyla","Phyla","Phyla","Phyla")

p.pr= ggplot(ngs.proot, aes(x=variable,y=prop, fill = phylum )) + 
  geom_bar(stat = "identity",position="stack",width=0.8)+ 
  #scale_y_continuous(labels = scales::percent) + 
  xlab("Wheat root mycobiota in field")+
  ylab("Percentage of sequence clusters")+ theme_classic()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values =rev(c("Multi-affiliation"="#000000","GZygomycota"="#44A8DB",
                                  "Glomeromycota"="#8DD3C7","Chytridiomycota"="#57B78C",
                                  "Basidiomycota"="#BEBADA","Ascomycota"="#3C3A8D"))) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())
p.pr

p.rr= ggplot(ngs.proot, aes(x=variable,y=richness, fill = phylum )) + 
  geom_bar(stat = "identity",position="stack",width=0.8)+ 
  #scale_y_continuous(labels = scales::percent) + 
  xlab("Wheat root mycobiota in field")+
  ylab("Richness of sequence clusters")+ theme_classic()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values =rev(c("Multi-affiliation"="#000000","GZygomycota"="#44A8DB",
                                  "Glomeromycota"="#8DD3C7","Chytridiomycota"="#57B78C",
                                  "Basidiomycota"="#BEBADA","Ascomycota"="#3C3A8D"))) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())
p.rr

pdf("Figure S3 new.pdf",height=4,width=12)
grid.arrange(p.ps,p.rs,p.pr,p.rr,ncol=4)
dev.off() 

ngs.sasco<-subset(ngs.btsoil,taxo$Phylum=="Ascomycota")
ngs.sb<-subset(ngs.btsoil,taxo$Phylum=="Basidiomycota")
ngs.sc<-subset(ngs.btsoil,taxo$Phylum=="Chytridiomycota")
ngs.sg<-subset(ngs.btsoil,taxo$Phylum=="Glomeromycota")
ngs.sz<-subset(ngs.btsoil,taxo$Phylum=="Zygomycota")
ngs.sm<-subset(ngs.btsoil,taxo$Phylum=="Multi-affiliation")

ngs.rasco<-subset(ngs.btroot,taxo$Phylum=="Ascomycota")
ngs.rb<-subset(ngs.btroot,taxo$Phylum=="Basidiomycota")
ngs.rc<-subset(ngs.btroot,taxo$Phylum=="Chytridiomycota")
ngs.rg<-subset(ngs.btroot,taxo$Phylum=="Glomeromycota")
ngs.rz<-subset(ngs.btroot,taxo$Phylum=="Zygomycota")
ngs.rm<-subset(ngs.btroot,taxo$Phylum=="Multi-affiliation")

root.cluster<-as.data.frame(rowSums(ngs.btroot[,1:60]))
write.table(root.cluster,file="root.cluster.txt",sep="\t")

mean(richness.all1$otu.richness.soil)
sd(richness.all1$otu.richness.soil)

mean(richness.all1$otu.richness.root)
sd(richness.all1$otu.richness.root)

mean(richness.all1$otu.richness.f)
sd(richness.all1$otu.richness.f)

x<-cbind(as.data.frame(richness.all1$otu.richness.soil),as.data.frame(richness.all1$otu.richness.root))
colnames(x)<-c("soil","root")
rsr<-as.data.frame(stack(x))
colnames(rsr)<-c("richness","sample")

rich.sr<-ggplot(rsr, aes(x=sample, y=richness))+ 
  geom_point(position=position_dodge(width=1))+ 
  geom_boxplot(aes(fill=as.factor(sample)))+
  ylim(0,350)+
  labs(y='Mycobiota sequence-clusters richness',x='')+
  ggprism::theme_prism()+theme(axis.text.x=element_text(angle = 0))


y<-cbind(as.data.frame(richness.all1$otu.evenness.soil),as.data.frame(richness.all1$otu.evenness.root))
colnames(y)<-c("soil","root")
esr<-as.data.frame(stack(y))
colnames(esr)<-c("evenness","sample")

even.sr<-ggplot(esr, aes(x=sample, y=evenness))+ 
  geom_point(position=position_dodge(width=1))+ 
  geom_boxplot(aes(fill=as.factor(sample)))+
  ylim(0,1)+
  labs(y='Mycobiota sequence-clusters evenness',x='')+
  ggprism::theme_prism()+theme(axis.text.x=element_text(angle = 0))

pdf("Figure S4.pdf",height=4.5,width=12)
grid.arrange(rich.sr,even.sr,ncol=2)
dev.off() 

#### Will try to get some figures ####
library(ggplot2)

pdf("Weed richness v.s wheat root mycobiota richness.pdf",height=6,width=12)
fr<-ggplot(richness.all1,aes(x=otu.richness.f,y=otu.richness.root))+#geom_boxplot()+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  #geom_hline(yintercept=0.3797,color="blue")+
  #ylim(0.4,0.7)+
  #geom_hline(yintercept=0.7008,color="red")+
  labs(x='Weed species richness',y='Mycobiota diversity of wheat roots')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))
dev.off()


pdf("Weed richness v.s soil mycobiota richness.pdf",height=6,width=12)
fs<-ggplot(richness.all1,aes(x=otu.richness.f,y=otu.richness.soil))+#geom_boxplot()+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  #geom_hline(yintercept=0.3797,color="blue")+
  #ylim(0.4,0.7)+
  #geom_hline(yintercept=0.7008,color="red")+
  labs(x='weed species richness',y='Mycobiota diversity of soil')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))
dev.off()

pdf("Richness of soil mycobiota v.s root mycobiota.pdf",height=6,width=15)
sr<-ggplot(richness.all1,aes(x=otu.richness.soil,y=otu.richness.root))+
  geom_point(position=position_dodge(width=0.75))+
  geom_smooth(method = 'lm',se = FALSE)+
  #geom_hline(yintercept=0.3797,color="blue")+
  #ylim(0.4,0.7)+
  #geom_hline(yintercept=0.7008,color="red")+
  labs(x='Mycobiota diversity of soil',y='Mycobiota diversity of wheat roots')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))
dev.off()

library(gridExtra)

pdf("Floristic richness, soil mycobiota richness and root mycobiota richness.pdf",height=4,width=12)
grid.arrange(fr,sr,fs,ncol=3)
dev.off() 


#### loading the packages #### 
library(nlme)
library(lme4)
library(MuMIn)
library(car)
library(lmerTest)
library(lsmeans)
library(MASS)
library(rsq)

richness.all1<-read.table("richness.all.i.txt", header=T, row.names=1, sep="\t")

lm1<-lm(otu.richness.root~otu.richness.f+otu.richness.soil,richness.all1)
anova(lm1)

lm1<-lm(otu.richness.root~otu.richness.f,richness.all1)
anova(lm1)

lm1<-lm(otu.richness.soil~otu.richness.f,richness.all1)
anova(lm1)


#mod<-lmer(otu.evenness.root~otu.richness.f*otu.richness.soil+(1|field)+(1|quadrats), data=richness.all1)
#mod<-glmer.nb(otu.richness.soil~otu.richness.f+(1|field)+(1|quadrats), data=richness.all1)
mod<-glmer.nb(otu.richness.root~otu.richness.f+otu.richness.soil+(1|field)+(1|quadrats), data=richness.all1)
#summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

####group comparison after model corrections
lsmeans(mod,pairwise~otu.shannon.f,by="otu.shannon.soil",data=pcoa.root.r)
lsmeans(mod,pairwise~otu.shannon.soil,by="otu.shannon.f",data=richness.all1)

mod<-lmer(otu.evenness.root~otu.richness.f+otu.richness.soil+(1|field)+(1|quadrats), data=richness.all1)
#summary(mod)
Anova(mod)
AIC(mod)
r.squaredGLMM(mod)
hist(residuals(mod))#check normality of residuals

####group comparison after model corrections
lsmeans(mod,pairwise~otu.shannon.f,by="otu.shannon.soil",data=pcoa.root.r)
lsmeans(mod,pairwise~otu.shannon.soil,by="otu.shannon.f",data=richness.all1)


