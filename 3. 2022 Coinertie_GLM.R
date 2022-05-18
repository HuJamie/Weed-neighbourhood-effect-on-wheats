#####################################################################################################################
#####################################################################################################################
######
###### Coinertie entre composition et variables paysag?res
library(ade4)
#library(omicade4)

#### Reading the data ####
ngs.soil.i<-read.table("ngs.soil.ir.txt", header=T, row.names=1, sep="\t")  
ngs.root.i<-read.table("ngs.root.ir.txt", header=T, row.names=1, sep="\t")  
floristic.i<-read.table("floristic.txt", header=T, row.names=1, sep="\t")

#ngs.soil.i<-ngs.soil[(row.names(ngs.soil) %in% row.names(richness.all1)), ]
#ngs.root.i<-ngs.root[(row.names(ngs.root) %in% row.names(richness.all1)), ]
floristic.ri<-floristic.i[,(colnames(floristic.i) %in% row.names(ngs.root.i))]
floristic.si<-floristic.i[,(colnames(floristic.i) %in% row.names(ngs.soil.i))]

## CIA 
#AFCE=dudi.coa(ngs.broot,scannf=F,nf=2)
AFCE=dudi.pca(ngs.root.i,scannf=F,nf=2)
AFCE

# ACP landscape (centered reduced by default) by weighting the lines by AFC bota
#ACPP=dudi.pca(floristic.root,row.w=AFCE$lw,scannf=F,nf=2)
ACPP=dudi.pca(t(floristic.ri),scannf=F,nf=2)
ACPP

# Do the CoInertia analysis
Coi=coinertia(ACPP,AFCE,scann=F,nf=2)
Coi

summary(Coi)
# Tester l'analyse de CoInertie (RV= ?observations)
test=randtest(Coi,nrepet = 999,fixed=2) #Performs a Monte-Carlo test on a Co-inertia analysis
#plot(test)
test

# Repr?sentation graphique de la CoInertie
plot(Coi)
plot(Coi$co)

plot(Comp1~Comp2, main = 'Floristic species v.s root', data = Coi$co)
with(Coi$co, text(Comp1~Comp2, labels=row.names(Coi$co), pos=3))

coiner.co1<-read.table("coiner.co1.txt", header=T, row.names=1, sep="\t") 
coiner.cocs<-read.table("coiner.cocs.txt", header=T, row.names=1, sep="\t")

s.arrow(coiner.co1[,1:2], clab = 0.7)
s.arrow(coiner.co1[,3:4], clab = 0.7)
plot(coiner.co1[,3:4])

s.arrow(coiner.cocs[,1:2], clab = 0.5)
plot(coiner.cocs[,1:2])

plot(CS1~CS2, main = 'Floristic species v.s root', data = Coi$c1)
with(Coi$c1, text(CS1~CS2,labels=row.names(Coi$c1), pos=3))

# Inertia on the CoInertia axes 
Coi$eig[1]/sum(Coi$eig)
Coi$eig[2]/sum(Coi$eig)
inertia.dudi(Coi)

#### coinertia for root Ascomycota####
ngs.asco<-subset(t(ngs.root.i),taxo$Phylum=="Ascomycota")

AFCE.a=dudi.pca(t(ngs.asco),scannf=F,nf=2)
AFCE.a

ACPP=dudi.pca(t(floristic.ri),scannf=F,nf=2)
ACPP

# Do the CoInertia analysis
Coi.a=coinertia(ACPP,AFCE.a,scann=F,nf=2)
Coi.a

summary(Coi.a)
# Tester l'analyse de CoInertie (RV= ?observations)
test=randtest(Coi.a,fixed=2) #Performs a Monte-Carlo test on a Co-inertia analysis
test

# Inertia on the CoInertia axes 
Coi.a$eig[1]/sum(Coi.a$eig)
Coi.a$eig[2]/sum(Coi.a$eig)
inertia.dudi(Coi.a)

#### coinertia for root Basidiomycota####
ngs.basi<-subset(t(ngs.root.i),taxo$Phylum=="Basidiomycota")

AFCE.b=dudi.pca(t(ngs.basi),scannf=F,nf=2)
AFCE.b

ACPP=dudi.pca(t(floristic.ri),scannf=F,nf=2)
ACPP

# Do the CoInertia analysis
Coi.b=coinertia(ACPP,AFCE.b,scann=F,nf=2)
Coi.b

summary(Coi.b)
# Tester l'analyse de CoInertie (RV= ?observations)
test=randtest(Coi.b,fixed=2) #Performs a Monte-Carlo test on a Co-inertia analysis
test

# Inertia on the CoInertia axes 
Coi.b$eig[1]/sum(Coi.b$eig)
Coi.b$eig[2]/sum(Coi.b$eig)
inertia.dudi(Coi.b)

#### coinertia for root Chytridiomycota ####
ngs.chyt<-subset(t(ngs.root.i),taxo$Phylum=="Chytridiomycota")

AFCE.c=dudi.pca(t(ngs.chyt),scannf=F,nf=2)
AFCE.c

ACPP=dudi.pca(t(floristic.ri),scannf=F,nf=2)
ACPP

# Do the CoInertia analysis
Coi.c=coinertia(ACPP,AFCE.c,scann=F,nf=2)
Coi.c

summary(Coi.c)
# Tester l'analyse de CoInertie (RV= ?observations)
test=randtest(Coi.c,fixed=2) #Performs a Monte-Carlo test on a Co-inertia analysis
test

# Inertia on the CoInertia axes 
Coi.c$eig[1]/sum(Coi.c$eig)
Coi.c$eig[2]/sum(Coi.c$eig)
inertia.dudi(Coi.c)

#### coinertia for root Glomeromycota ####
ngs.glom<-subset(t(ngs.root.i),taxo$Phylum=="Glomeromycota")

AFCE.g=dudi.pca(t(ngs.glom),scannf=F,nf=2)
AFCE.g

ACPP=dudi.pca(t(floristic.ri),scannf=F,nf=2)
ACPP

# Do the CoInertia analysis
Coi.g=coinertia(ACPP,AFCE.g,scann=F,nf=2)
Coi.g

summary(Coi.g)
# Tester l'analyse de CoInertie (RV= ?observations)
test=randtest(Coi.g,fixed=2) #Performs a Monte-Carlo test on a Co-inertia analysis
test

# Inertia on the CoInertia axes 
Coi.g$eig[1]/sum(Coi.g$eig)
Coi.g$eig[2]/sum(Coi.g$eig)
inertia.dudi(Coi.g)

#### coinertia for root Zygomycota ####
ngs.zygo<-subset(t(ngs.root.i),taxo$Phylum=="Zygomycota")

AFCE.z=dudi.pca(t(ngs.zygo),scannf=F,nf=2)
AFCE.z

ACPP=dudi.pca(t(floristic.ri),scannf=F,nf=2)
ACPP

# Do the CoInertia analysis
Coi.z=coinertia(ACPP,AFCE.z,scann=F,nf=2)
Coi.z

summary(Coi.z)
# Tester l'analyse de CoInertie (RV= ?observations)
test=randtest(Coi.z,fixed=2) #Performs a Monte-Carlo test on a Co-inertia analysis
test

# Inertia on the CoInertia axes 
Coi.z$eig[1]/sum(Coi.z$eig)
Coi.z$eig[2]/sum(Coi.z$eig)
inertia.dudi(Coi.z)

######## coinertia for soil #########
AFCE.s=dudi.pca(ngs.soil.i,scannf=F,nf=2)
AFCE.s

#ACPP=dudi.pca(floristic.root,row.w=AFCE$lw,scannf=F,nf=2)
ACPP.s=dudi.pca(t(floristic.si),scannf=F,nf=2)
ACPP.s

Coi.s=coinertia(ACPP.s,AFCE.s,scann=F,nf=2)
Coi.s

summary(Coi.s)
# Tester l'analyse de CoInertie (RV= ?observations)
test.s=randtest(Coi.s,fixed=2) #Performs a Monte-Carlo test on a Co-inertia analysis
test.s

Coi.s$co
Coi.s$c1

# Inertie sur les axes de la CoInertie
# Inertia on the CoInertia axes 
# CoInertia 轴上的惯性
Coi.s$eig[1]/sum(Coi.s$eig)
Coi.s$eig[2]/sum(Coi.s$eig)
inertia.dudi(Coi.s)

#### coinertia for soil Ascomycota####
ngs.asco.s<-subset(t(ngs.soil.i),taxo$Phylum=="Ascomycota")

AFCE.a=dudi.pca(t(ngs.asco.s),scannf=F,nf=2)
AFCE.a

ACPP=dudi.pca(t(floristic.si),scannf=F,nf=2)
ACPP

# Do the CoInertia analysis
Coi.a=coinertia(ACPP,AFCE.a,scann=F,nf=2)
Coi.a

summary(Coi.a)
# Tester l'analyse de CoInertie (RV= ?observations)
test=randtest(Coi.a,fixed=2) #Performs a Monte-Carlo test on a Co-inertia analysis
test

# Inertia on the CoInertia axes 
Coi.a$eig[1]/sum(Coi.a$eig)
Coi.a$eig[2]/sum(Coi.a$eig)
inertia.dudi(Coi.a)

#### coinertia for soil Basidiomycota####
ngs.basi.s<-subset(t(ngs.soil.i),taxo$Phylum=="Basidiomycota")

AFCE.b=dudi.pca(t(ngs.basi.s),scannf=F,nf=2)
AFCE.b

ACPP=dudi.pca(t(floristic.si),scannf=F,nf=2)
ACPP

# Do the CoInertia analysis
Coi.b=coinertia(ACPP,AFCE.b,scann=F,nf=2)
Coi.b

summary(Coi.b)
# Tester l'analyse de CoInertie (RV= ?observations)
test=randtest(Coi.b,fixed=2) #Performs a Monte-Carlo test on a Co-inertia analysis
test

# Inertia on the CoInertia axes 
Coi.b$eig[1]/sum(Coi.b$eig)
Coi.b$eig[2]/sum(Coi.b$eig)
inertia.dudi(Coi.b)

#### coinertia for spol Chytridiomycota ####
ngs.chyt.s<-subset(t(ngs.soil.i),taxo$Phylum=="Chytridiomycota")

AFCE.c=dudi.pca(t(ngs.chyt.s),scannf=F,nf=2)
AFCE.c

ACPP=dudi.pca(t(floristic.si),scannf=F,nf=2)
ACPP

# Do the CoInertia analysis
Coi.c=coinertia(ACPP,AFCE.c,scann=F,nf=2)
Coi.c

summary(Coi.c)
# Tester l'analyse de CoInertie (RV= ?observations)
test=randtest(Coi.c,fixed=2) #Performs a Monte-Carlo test on a Co-inertia analysis
test

# Inertia on the CoInertia axes 
Coi.c$eig[1]/sum(Coi.c$eig)
Coi.c$eig[2]/sum(Coi.c$eig)
inertia.dudi(Coi.c)

#### coinertia for soil Glomeromycota ####
ngs.glom.s<-subset(t(ngs.soil.i),taxo$Phylum=="Glomeromycota")

AFCE.g=dudi.pca(t(ngs.glom.s),scannf=F,nf=2)
AFCE.g

ACPP=dudi.pca(t(floristic.si),scannf=F,nf=2)
ACPP

# Do the CoInertia analysis
Coi.g=coinertia(ACPP,AFCE.g,scann=F,nf=2)
Coi.g

summary(Coi.g)
# Tester l'analyse de CoInertie (RV= ?observations)
test=randtest(Coi.g,fixed=2) #Performs a Monte-Carlo test on a Co-inertia analysis
test

# Inertia on the CoInertia axes 
Coi.g$eig[1]/sum(Coi.g$eig)
Coi.g$eig[2]/sum(Coi.g$eig)
inertia.dudi(Coi.g)

#### coinertia for soil Zygomycota ####
ngs.zygo.s<-subset(t(ngs.soil.i),taxo$Phylum=="Zygomycota")

AFCE.z=dudi.pca(t(ngs.zygo.s),scannf=F,nf=2)
AFCE.z

ACPP=dudi.pca(t(floristic.si),scannf=F,nf=2)
ACPP

# Do the CoInertia analysis
Coi.z=coinertia(ACPP,AFCE.z,scann=F,nf=2)
Coi.z

summary(Coi.z)
# Tester l'analyse de CoInertie (RV= ?observations)
test=randtest(Coi.z,fixed=2) #Performs a Monte-Carlo test on a Co-inertia analysis
test

# Inertia on the CoInertia axes 
Coi.z$eig[1]/sum(Coi.z$eig)
Coi.z$eig[2]/sum(Coi.z$eig)
inertia.dudi(Coi.z)

