#install.packages("ape")
#install.packages("phangorn")
#install.packages("picante")
require(ape)
require(phangorn)
require(picante)
require(phytools)
#Read in the data from the fasta file.
x<-read.dna("marsup_rag1.fasta",format="fasta")
#change the rownames to account for doubles
rownames(x)<-((1:103))

#Compare the number of substitutions seperating any pair of sequences(distance), this uses deafult K80
d<-dist.dna(x)

#Write this into a viewable table
write.table(as.matrix(d),"distances1.csv")
#
#Plotting trees from distance matrices
#
#choose method of construction(bionj in this case)
tr.bionj<-bionj(d)
#plot the tree
plot(tr.bionj)
#root the tree at the node=outgroup
tr.bionjr<-root(tr.bionj,outgroup="1", resolve.root=TRUE )
#plot the rooted tree
plot(tr.bionjr);add.scale.bar(length=0.001)
#
#Tree distortion
#
#Calculate the distortion using cophenetic analysis
dt.bionj<-cophenetic(tr.bionj)
#pull out the taxa data as a matrix
dmat<-as.matrix(d)
#pul out the rownames of dmat and assign value
nms<-rownames(dmat)
#update the cophenetic analysis with the rownames
dt.bionj<-dt.bionj[nms, nms]
#Calculate the distances
dt.bionj<-as.dist(dt.bionj)
#plot new distance - old distance = residuals
plot(dt.bionj-d,ylab="residuals", cex=0.5,main="BIONJ")
abline(h=0,lty=3)
#

#Redo for nj tree
#Plotting trees from distance matrices
#
#choose method of construction(nj in this case)
tr.nj<-nj(d)
#plot the tree
plot(tr.nj)
#root the tree at the node=outgroup
tr.njr<-root(tr.nj,outgroup="1", resolve.root=TRUE )
#plot the rooted tree
plot(tr.nj);add.scale.bar(length=0.001)
#
#
#
#Tree distortion
#
#
#Calculate the distortion using cophenetic analysis
dt.nj<-cophenetic(tr.nj)
#pull out the taxa data as a matrix
dmat<-as.matrix(d)
#pul out the rownames of dmat and assign value
nms<-rownames(dmat)
#update the cophenetic analysis with the rownames
dt.nj<-dt.nj[nms, nms]
#Calculate the distances
dt.nj<-as.dist(dt.nj)
#plot new distance - old distance = residuals
plot(dt.nj-d,ylab="residuals", cex=0.5,main="NJ")
abline(h=0,lty=3)
#
#
#redo with upgma
#choose method of construction(upgma in this case)
tr.upgma<-upgma(d)
#plot the tree
plot(tr.upgma)
#root the tree at the node=outgroup
tr.upgmar<-root(tr.upgma,outgroup="1", resolve.root=TRUE )
#plot the rooted tree
plot(tr.upgmar);add.scale.bar(length=0.001)
#
#
#
#Tree distortion
#
#
#Calculate the distortion using cophenetic analysis
dt.upgma<-cophenetic(tr.upgma)
#pull out the taxa data as a matrix
dmat<-as.matrix(d)
#pul out the rownames of dmat and assign value
nms<-rownames(dmat)
#update the cophenetic analysis with the rownames
dt.upgma<-dt.upgma[nms, nms]
#Calculate the distances
dt.upgma<-as.dist(dt.upgma)
#plot new distance - old distance = residuals
plot(dt.upgma-d,ylab="residuals", cex=0.5,main="upgma")
abline(h=0,lty=3)
#
#compare the 3
par(mfrow=c(1,3))
plot(dt.upgma-d,ylab="residuals", cex=0.5,main="upgma")
abline(h=0,lty=3)
plot(dt.nj-d,ylab="residuals", cex=0.5,main="nj")
abline(h=0,lty=3)
plot(dt.bionj-d,ylab="residuals", cex=0.5,main="bionj")
abline(h=0,lty=3)
#Bionj has lowest residuals
#
#BOOTSTRAPPING
#
#Fit the tree to the data
fit<-pml(tr.bionj,as.phyDat(x))
#optimise plot and set random seed for bootstrap process
fit=optim.pml(fit,T)
plot(fit)
set.seed(8)
#bootstrap the data and plot the results
bs<-bootstrap.pml(fit,bs=100,optNni=T)
treeBS<-plotBS(fit$tree, type="p", bs)
#
#Substitution Models: find the most appropriate model
mt<-modelTest(as.phyDat(x),G=F,I=F)
View(mt)#pick lowest AIC values
#lowest AIC value is GTR model
#This is the best model, ,recalculate distance with this method.
#
fittedtreeGTR<-pml(tr.bionj,as.phyDat(x),k=4,inv=0.2)
plot(fittedtreeGTR)
#boostrap
bs<-bootstrap.pml(fittedtreeGTR,bs=100,optNni=T)
plotBS(fittedtreeGTR$tree,type='p',bs)
#
#cophenetic analysis
dt.fittedtreeGTR<-cophenetic(fittedtreeGTR$tree)
dt.fittedtreeGTR<-dt.fittedtreeGTR[nms, nms]
#Calculate the distances
dt.fittedtreeGTR<-as.dist(dt.fittedtreeGTR)
#plot new distance - old distance = residuals
plot(dt.fittedtreeGTR-d,ylab="residuals", cex=0.5,main="GTR")
abline(h=0,lty=3)

#BAYES
#Was done on xming
#ED scores
#
#Calculate the evolutionary distinctiveness of ech gene (must be done on a rooted tree)
orig<-evol.distinct(tr.bionjr,type="fair.proportion")
orig
View(orig)
#plot histogram
 plot(orig, col=" blue", border=" blue" ,main="ED score", density=100)
