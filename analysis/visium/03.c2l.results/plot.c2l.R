library(Seurat)
library(visutils) # devtools::install_github("iaaka/visutils")
library(NMF)
library(plyr)
library(doMC)
source('src/my/visium/c2l.utils.R')

# load data ########
meta = read.csv('src/my/visium/samples.csv')
rownames(meta) = meta$Sanger_id
meta = meta[order(meta$PCD),]

vsf = readRDS('data.nfs/my/visium/viss.all.filtered.rds')
#c2l_sd1 = loadC2L('data.nfs/my/visium/skin.c2l/pred/fetal_skin.norm.maternal_removed.20220202.filtered.denorm.20/predmodel/',1)
c2l_sd0 = loadC2L('data.nfs/my/visium/skin.c2l/pred/fetal_skin.norm.maternal_removed.20220202.filtered.denorm.20/predmodel/',0)

vsf = vsf[meta$Sanger_id]
#c2l_sd1 = c2l_sd1[meta$Sanger_id]
c2l_sd0 = c2l_sd0[meta$Sanger_id]

for(n in names(vsf))
  vsf[[n]] = vsf[[n]][,rownames(c2l_sd0[[n]])]

c2l = c2l_sd0

# plot summary ############
ctcols = getColoursByDistance(1-cor(do.call(rbind,c2l)),orderBySim = T,use3D = F)
# to order in legend
cells.leg.order = names(ctcols)
ctcols = ctcols[colnames(c2l[[1]])]


ncol=5
width = ncol*3
height = 3*3
he.img.width = 400

pdf('figures/visium/sd1.meancol.summary.pdf',w=width,h=height)
par(mfrow=c(3,ncol),mar=c(0.1,0.1,1.2,4),bty='n',oma=c(0,0,0,0))
for(n in names(vsf)){
  plotVisiumMultyColours(vsf[[n]],c2l[[n]],cols = ctcols,mode = 'mean',
                         zfun = function(x)x^2,he.img.width=he.img.width,
                         scale.per.colour = T,legend.ncol=0,min.opacity = 250,
                         main=paste0(n,' (',meta[n,'GA'],')'))
}
plot.new()
legend('topleft',xpd=NA,fill=ctcols[cells.leg.order],border=NA,legend=cells.leg.order,ncol=4,bty='n')
dev.off()


# 2F Fb ########
cts = c('WNT2+ fibroblast','HOXC5+ early fibroblast')
cols.fb = setNames(c('#2cb8c9','#bb7784'),cts)
img.alpha = 0.6


pdf('figures/visium/WNT2-HOXC5-Fbs.visium.all.samples.pdf',w=4*5,h=2*3+0.2)
par(mfcol=c(2,5),mar=c(0.1,0.1,1.5,10),bty='n',tcl=-0.2,mgp=c(1.3,0.3,0),oma=c(0,0,1.3,0))
for(n in meta$Sanger_id){
    plotVisiumMultyColours(vsf[[n]],c2l[[n]][,cts],cols = cols.fb[cts],scale.per.colour = T,zfun=function(x)x,mode = 'mean',bg = '#FFFFFFFF',min.opacity = 255,
                           border=NA,img.alpha=img.alpha,title.adj = c(0,-1),main=paste0(n,' (',meta[n,'PCD'],'pcd)'))
    d = c2l[[n]][,cts]
    #d = sweep(d,2,apply(d,2,max),'/')
    
    par(mar=c(3,3,1.5,10))
    plotCellAbundanceOnDist2Border(vsf[[n]]$dist2border.graph,d,log=F,cols=cols.fb[cts],xlim=range(vsf[[n]]$dist2border.graph))
    par(mar=c(0.1,0.1,1.5,10))  
}
dev.off()


# _with PEAR1+ ###########
cts = c('WNT2+ fibroblast','HOXC5+ early fibroblast','PEAR1+ fibroblast')
cols.fb = setNames(c('#2cb8c9','#bb7784','#4ad367'),cts)
img.alpha = 0.6


pdf('figures/visium/WNT2-HOXC5-PEAR1-Fbs.visium.all.samples.pdf',w=4*5,h=2*3+0.2)
par(mfcol=c(2,5),mar=c(0.1,0.1,1.5,10),bty='n',tcl=-0.2,mgp=c(1.3,0.3,0),oma=c(0,0,1.3,0))
for(n in meta$Sanger_id){
  plotVisiumMultyColours(vsf[[n]],c2l[[n]][,cts],cols = cols.fb[cts],scale.per.colour = T,zfun=function(x)x,mode = 'mean',bg = '#FFFFFFFF',min.opacity = 255,
                         border=NA,img.alpha=img.alpha,title.adj = c(0,-1),main=paste0(n,' (',meta[n,'PCD'],'pcd)'))
  d = c2l[[n]][,cts]
  d = sweep(d,2,apply(d,2,max),'/')
  par(mar=c(3,3,1.5,10))
  plotCellAbundanceOnDist2Border(vsf[[n]]$dist2border.graph,d,log=F,cols=cols.fb[cts],xlim=range(vsf[[n]]$dist2border.graph))
  par(mar=c(0.1,0.1,1.5,10))  
}
dev.off()


# pair corr (1E) #######
# From Hudaa:
# LYVE1++ macrophage/WNT2+ fibroblast
# TREM2+ macrophage/WNT2+ fibroblast
# ILC2/Pre-dermal condensate
# LTi/Pre-dermal condensate
# Iron−recycling macrophage/All endothelial cells
# LYVE1+ macrophage/All endothelial cells
# TREM2+ macrophage/All endothelial cells
# MHCII+ macropahage/All endothelial cells
# Iron−recycling macrophage/All neuronal cells
# LYVE1+ macrophage/All neuronal cells
# TREM2+ macrophage/All neuronal cells
# MHCII+ macropahage/All neuronal cells
end = c('Early endothelial cells','Capillary arterioles','Arterioles','Capillaries','Postcapillary venules')#,'Venules') there is no Venules <= 10pcw
setdiff(end,colnames(c2l_sd0$HCA_rFSKI13460601))
neu = c('Neuron progenitors', 'Schwann/Schwann precursors', 'PID1+ schwann cellls', 'Myelinating Schwann cells')
setdiff(neu,colnames(c2l_sd0$HCA_rFSKI13460601))

ctpairs = rbind(data.frame(ct1='LYVE1++ macrophage',ct2='WNT2+ fibroblast'),
                data.frame(ct1='TREM2+ macrophage',ct2='WNT2+ fibroblast'),
                data.frame(ct1='ILC2',ct2='Pre-dermal condensate'),
                data.frame(ct1='LTi',ct2='Pre-dermal condensate'),
                data.frame(ct1='Iron-recycling macrophage',ct2=end),
                data.frame(ct1='LYVE1++ macrophage',ct2=end),
                data.frame(ct1='TREM2+ macrophage',ct2=end),
                data.frame(ct1='MHCII+ macrophage',ct2=end),
                data.frame(ct1='Iron-recycling macrophage',ct2=neu),
                data.frame(ct1='LYVE1++ macrophage',ct2=neu),
                data.frame(ct1='TREM2+ macrophage',ct2=neu),
                data.frame(ct1='MHCII+ macrophage',ct2=neu))

# look on Pre−dermal condensate
ctpairs = data.frame(ct1='Pre-dermal condensate',ct2=colnames(c2l_sd0$HCA_rFSKI13460601))
                
ctpairs = ctpairs[nrow(ctpairs):1,]
setdiff(ctpairs$ct2,colnames(c2l_sd0$HCA_rFSKI13460601))
setdiff(ctpairs$ct1,colnames(c2l_sd0$HCA_rFSKI13460601))

ctpcor.per.sam = sapply(c2l,function(x){x=sweep(x,1,apply(x,1,sum),'/');apply(ctpairs,1,function(ct)cor(x[,ct[1]],x[,ct[2]]))})
cols=char2col(meta$PCD,palette = T)
ord = nrow(meta):1
pdf('figures/visium/celltype.pairs.cor.per.sample.pdf',w=8,h=14)
par(mar=c(4,12,1,1))
b=barplot(t(ctpcor.per.sam[,ord]),beside = T,border = F,col=cols[as.character(meta$PCD)[ord]],names.arg = paste0(ctpairs$ct1,'\n',ctpairs$ct2),las=2,space = c(0,5),cex.names=0.6,horiz = TRUE,xlab='Pearson correlation coefficient')
legend('topleft',fill=cols,legend=names(cols),border=NA,bty='n',title='PCD')
x = apply(b,2,mean)
#segments(x,-0.01,x,-1,lty=2,col='gray')
segments(-0.01,x,-1,x,lty=2,col='gray')
dev.off()



# look on Pre−dermal condensate
ctpairs = data.frame(ct1='Pre-dermal condensate',ct2=setdiff(colnames(c2l_sd0$HCA_rFSKI13460601),'Pre-dermal condensate'))
ctpcor.per.sam = sapply(c2l,function(x){x=sweep(x,1,apply(x,1,sum),'/');apply(ctpairs,1,function(ct)cor(x[,ct[1]],x[,ct[2]]))})
o = order(apply(ctpcor.per.sam,1,mean))
ctpairs = ctpairs[o,]
ctpcor.per.sam = ctpcor.per.sam[o,]

setdiff(ctpairs$ct2,colnames(c2l_sd0$HCA_rFSKI13460601))
setdiff(ctpairs$ct1,colnames(c2l_sd0$HCA_rFSKI13460601))

cols=char2col(meta$PCD,palette = T)
ord = nrow(meta):1

pdf('figures/visium/Pre-dermal_condensate_vs_all_pcc.pdf',w=8,h=14)
par(mar=c(4,12,1,1))
b=barplot(t(ctpcor.per.sam[,ord]),beside = T,border = F,col=cols[as.character(meta$PCD)[ord]],names.arg = ctpairs$ct2,las=2,space = c(0,5),cex.names=0.6,horiz = TRUE,xlab='Pearson correlation coefficient')
legend('topleft',fill=cols,legend=names(cols),border=NA,bty='n',title='PCD')
x = apply(b,2,mean)
#segments(x,-0.01,x,-1,lty=2,col='gray')
segments(-0.01,x,-1,x,lty=2,col='gray')
dev.off()

# _all pairs ######
nmf.frq10 = readRDS('data.nfs/my/visium/nmf.frq.r10.n100.rds')
nmfbst = nmfGetBest(nmf.frq10,nmfNormFs$max)
cols = RColorBrewer::brewer.pal(10,'Set3')
cols = setNames(cols,1:10)
colfun=function(x)num2col(x,c('blue','cyan','gray','orange','red'),minx = -1,maxx = 1)

#ord = c(1,10,3,4,7,9,2,5,6,8)
ord = 1:10
o = names(nmfbst$cl)[order(match(nmfbst$cl,ord),names(nmfbst$cl))]
c2l.cor = lapply(c2l,function(x)cor(x[,o]))


pdf('figures/visium/per.sample.pcc.dotplot.pdf',w=4*6,h=3*5.5)
par(mfrow=c(3,4),mar=c(13,13,1,6),las=2,bty='n')
for(n in rownames(meta)){
  x=c2l.cor[[n]]
  dotPlot(x,rfun = function(x)(1+x),scaleWM = F,colfun = colfun,colColours = cols[nmfbst$cl[o]],rowColours = cols[nmfbst$cl[o]],max.cex = 0.5,ylab.cex = 0.3,xlab.cex = 0.3,
          main=paste0(meta[n,'site'],' ',round(meta[n,'PCD']/7),'pcw (',n,')'),legend.cex.at = c(-0.5,0,0.5,1),grid = F,legend.cex.title = 'PCC')
}
dev.off()
