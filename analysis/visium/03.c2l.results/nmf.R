library(Seurat)
library(visutils) # devtools::install_github("iaaka/visutils")
library(NMF)
library(plyr)
library(doMC)
source('src/my/visium/c2l.utils.R')

# load data ########
nmf.frq10 = readRDS('data.nfs/my/visium/nmf.frq.r10.n100.rds')
meta = read.csv('src/my/visium/samples.csv',row.names = 1)

vsf = readRDS('data.nfs/my/visium/viss.all.filtered.rds')
#c2l_sd1 = loadC2L('data.nfs/my/visium/skin.c2l/pred/fetal_skin.norm.maternal_removed.20220202.filtered.denorm.20/predmodel/',1)
c2l_sd0 = loadC2L('data.nfs/my/visium/skin.c2l/pred/fetal_skin.norm.maternal_removed.20220202.filtered.denorm.20/predmodel/',0)

vsf = vsf[meta$Sanger_id]
#c2l_sd1 = c2l_sd1[meta$Sanger_id]
c2l_sd0 = c2l_sd0[meta$Sanger_id]

for(n in names(vsf))
  vsf[[n]] = vsf[[n]][,rownames(c2l_sd0[[n]])]

c2l = c2l_sd0
# NMF #####
c2la = do.call(rbind,lapply(names(c2l),function(n){x=c2l[[n]];rownames(x)=paste0(n,'|',rownames(x));x}))
c2lf = sweep(c2la,1,apply(c2la,1,sum),'/')
# makes  sense for sd1
f = apply(is.na(c2lf),1,sum) > 0
table(f)
c2la[f,]
#c2lf[is.na(c2lf)] = 0
# paralellization
doMC::registerDoMC(9)
# repeat N times for rank microenvironments
N = 100 # set to 100
rank = 10
set.seed(1234)
# nmf.frq10 = llply(1:N,function(i){nmf(c2lf[!f,], rank = rank)},.parallel = T)
# saveRDS(nmf.frq10,'data.nfs/my/visium/nmf.frq.r10.n100.rds')
# nmf.abs10 = llply(1:N,function(i){nmf(c2la[!f,], rank = rank)},.parallel = T)
# saveRDS(nmf.frq10,'data.nfs/my/visium/nmf.abs.r10.n100.rds')
# 
# 
# c2lal = do.call(rbind,lapply(meta$Sanger_id[meta$site=='limb'],function(n){x=c2l[[n]];rownames(x)=paste0(n,'|',rownames(x));x}))
# c2lfl = sweep(c2lal,1,apply(c2lal,1,sum),'/')
# nmf.frq10.limb = llply(1:N,function(i){nmf(c2lfl, rank = rank)},.parallel = T)
# saveRDS(nmf.frq10.limb,'data.nfs/my/visium/nmf.frq.limb.r10.n100.rds')
# 
# nmf.abs10.limb = llply(1:N,function(i){nmf(c2lal, rank = rank)},.parallel = T)
# saveRDS(nmf.abs10.limb,'data.nfs/my/visium/nmf.abs.limb.r10.n100.rds')

# redo old nmf
# c2l_old_sd0 = read.csv('../2202.c2l.service/processed/2203.limb/c2l/pred/denorm.no.erythroid.v2.20/predmodel/q05_cell_abundance_w_sf.csv',row.names = 1,check.names = F)
# colnames(c2l_old_sd0)=sub('q05cell_abundance_w_sf_','',colnames(c2l_old_sd0))
# rownames(c2l_old_sd0) = sapply(strsplit(rownames(c2l_old_sd0),'_'),function(x)paste0(x[2],'_',x[3],'|',x[1]))
# 
# c2l_old_sd0
# nmf.abs10.limb_old = llply(1:N,function(i){nmf(c2l_old_sd0, rank = rank)},.parallel = T)
# saveRDS(nmf.abs10.limb_old,'data.nfs/my/visium/nmf.abs.limb_old.r10.n100.rds')


# choose run that is the most similar to consensus
nmfbst = nmfGetBest(nmf.frq10,nmfNormFs$max)
cols = RColorBrewer::brewer.pal(rank,'Set3')
cols = setNames(cols,1:rank)

# comp with previous version #######
cellRenamer = c('Postcapillary venule'='Venules',
  'Capillary/postcapillary venule'='Postcapillary venules',
  'Capillary (venular tip)'='Capillaries',
  'Tip cell (arterial)'='Capillary arterioles',
  'Arterial'='Arterioles',
  'Early endothelial cell'='Early endothelial cells')
nmf = readRDS('../2202.c2l.service/processed/2203.limb/rds/nmf.a20.seed4321.rds')
par(las=2,mar=c(4,20,1,1))
# this looks right
colfun2=function(x)num2col(x,c('yellow','orange','violet','black'),minx = 0,maxx = 1)
oldme=apply(nmf$coef.n2,2,which.max)
dotPlot(t(nmf$coef.n2),ylab.cex = 0.8,colColours = cols,rowColours = cols[oldme],max.cex = 2,scaleWM=F,colfun=colfun2,grid = F)

# rename 
oldme_ = oldme
f = names(oldme_) %in% names(cellRenamer) 
names(oldme_)[f] = cellRenamer[names(oldme_)[f]]
oldorder = names(oldme_)
oldme_ = oldme_[names(nmfbst$cl)]
t = table(nmfbst$cl,oldme_)
imageWithText(t)
clrenamer = getBestMatchByContigencyTable(t,useHungarian = TRUE)
newme_ = nmfbst$cl
newme_ = setNames(clrenamer$col[match(newme_,rownames(t))],names(newme_))

t = table(oldme_,newme_)
imageWithText(t)
clcmp = cbind(old=oldme_,new=newme_)[oldorder[oldorder %in% names(oldme_)],]
clcmp = clcmp[nrow(clcmp):1,]
clcmp = clcmp[order(clcmp[,1],clcmp[,2],decreasing = T),]
par(las=2,mar=c(4,20,1,1),cex=0.8)
imageWithText(t(clcmp),col=cols)

colfun2=function(x)num2col(x,c('yellow','orange','violet','black'),minx = 0,maxx = 1)
o = order(nmfbst$cl)
dotPlot(t(nmfbst$coefn[,o]),ylab.cex = 0.8,colColours = cols,rowColours = cols[nmfbst$cl[o]],max.cex = 2,scaleWM=F,colfun=colfun2,grid = F)

#dir.create('figures/c2l/v02.merg_spots.reydKC.rie20/nmf')
pdf('figures/visium/nmf.summary.dots.pdf',w=12,h=9)
plotNMFCons(coefs=nmfbst$coefn,clcols = cols,cons = nmfbst$cons/N,max.cex = 1.7,ylab.cex = 0.6)
dev.off()

sid = sapply(strsplit(rownames(nmfbst$basisn),'|',fixed = TRUE),'[',1)

me = lapply(split.data.frame(nmfbst$basisn,sid),function(x){
  rownames(x) = sapply(strsplit(rownames(x),'|',fixed = TRUE),'[',2)
  x
})
me = me[names(c2l_sd0)]



pdf('figures/visium/nmf.on.visiums.pdf',w=9*3,h=7*3)
par(mar=c(0.1,0.1,1.2,0),bty='n',oma=c(0,0,0,0))
for(i in m$id){
  if(i %in% names(vsf)){
    cex = scaleTo(sqrt(vsf[[i]]$nspots),1,sqrt(8),minx = 1,maxx = sqrt(7))*0.9
    plotVisium(vsf[[i]],pie.fracs=me[[i]], pie.cols=cols,main=paste0(i,', ',m[i,'File.ID']),he.img.width=400,he.grayscale =TRUE,cex=cex)
  }else
    plot.new()
}
plot.new()
legend('topleft',bty='n',legend = 1:N,fill=cols,border=NA,ncol=1)
dev.off()

