library(Seurat)
library(visutils) # devtools::install_github("iaaka/visutils")
library(NMF)
source('src/my/visium/c2l.utils.R')

# load data ########
meta = read.csv('src/my/visium/samples.csv')
rownames(meta) = meta$Sanger_id
meta = meta[order(meta$PCD),]
meta$PCWclean = round(meta$PCD/7)

nmf.frq10 = readRDS('data.nfs/my/visium/nmf.frq.r10.n100.rds')
nmfbst = nmfGetBest(nmf.frq10,nmfNormFs$max)


vsf = readRDS('data.nfs/my/visium/viss.all.filtered.rds')
c2l_sd0 = loadC2L('data.nfs/my/visium/skin.c2l/pred/fetal_skin.norm.maternal_removed.20220202.filtered.denorm.20/predmodel/',0)

vsf = vsf[meta$Sanger_id]
c2l_sd0 = c2l_sd0[meta$Sanger_id]

c2l = c2l_sd0
# 2g #########
sids = c('WSSS_THYst9383362','WSSS_THYst9383359','HCA_rFSKI13460603')
cts = c('WNT2+ fibroblast','HOXC5+ early fibroblast')
cols.fb = setNames(c('#2cb8c9','#bb7784'),cts)
img.alpha = 0.6


vs2plot = list(rotateVisium(vsf[[sids[1]]]),rotateVisium(vsf[[sids[2]]]),rotateVisium(vsf[[sids[3]]],n=3))
names(vs2plot) = sids


cairo_pdf('figures/paper.figures/2g.WNT2-HOXC5-Fbs.pdf',w=4*3,h=2*3+0.2)
par(mfcol=c(2,3),mar=c(0.1,0.1,1.5,10),bty='n',tcl=-0.2,mgp=c(1.3,0.3,0),oma=c(0,0,1.3,0))
for(n in names(vs2plot)){
    plotVisiumMultyColours(vs2plot[[n]],c2l[[n]][,cts],cols = cols.fb[cts],scale.per.colour = T,zfun=function(x)x,mode = 'mean',bg = '#FFFFFFFF',min.opacity = 255,
                           border=NA,img.alpha=img.alpha,title.adj = c(0,-1),main=paste0(n,' (',meta[n,'PCD'],'pcd)'))
    d = c2l[[n]][,cts]
    
    par(mar=c(3,3,1.5,10))
    plotCellAbundanceOnDist2Border(vsf[[n]]$dist2border.graph,d,log=F,cols=cols.fb[cts],xlim=range(vsf[[n]]$dist2border.graph))
    par(mar=c(0.1,0.1,1.5,10))  
}
dev.off()


# celltype pairs #######
cs = colnames(c2l$WSSS_THYst9699523)
endo  = c('Early endothelial cells','Capillary arterioles','Arterioles','Capillaries','Postcapillary venules')#,'Venules') there is no Venules <= 10pcw
neuro = c('Neuron progenitors', 'Schwann/Schwann precursors', 'PID1+ schwann cellls', 'Myelinating Schwann cells')
macro = sort(cs[grep('macro',cs,ignore.case = T)],decreasing = F)[c(1,2,4,3)]
fibro = sort(cs[grep('fibro',cs,ignore.case = T)],decreasing = T)


c2l.cor = array(NA,c(ncol(c2l$WSSS_THYst9699523),ncol(c2l$WSSS_THYst9699523),length(c2l)),
                dimnames = list(colnames(c2l$WSSS_THYst9699523),colnames(c2l$WSSS_THYst9699523),names(c2l)))
for(n in names(c2l))
  c2l.cor[,,n] = cor(sweep(c2l[[n]],1,apply(c2l[[n]],1,sum),'/'))
 
hist(c2l.cor[c2l.cor<1],main='Celltype pair PCC',xlab='Celltype pair PCC')


plotPairCorr = function(cor,pairs,pcw,pcw2col,space=c(0,2),Nempty=0,...){
  d = do.call(cbind,lapply(1:nrow(pairs),function(i)cor[pairs[[1]][i],pairs[[2]][i],]))
  if(Nempty>0)
    for(i in 1:Nempty){
      d = cbind(d,0)
    }
  pcwo = order(pcw,decreasing = TRUE)
  d = d[pcwo,]
  pcw = pcw[pcwo]
  colnames(d) = c(paste0(pairs[[1]],'\n',pairs[[2]]),rep('',Nempty))
  b = barplot(d,beside = T,border = F,space = space,horiz = TRUE,las=1,
          col=cols[as.character(pcw)],...)
  x = apply(b,2,mean)[colnames(d)!='']
  y = pmin(0,apply(d,2,min))[colnames(d)!='']
  segments(y,x,-1,x,lty=2,col='gray')
  invisible(list(d=d,b=b))
}

plotPairCorrHM = function(cor,ct1,ct2,main='mean PCC across visium samples [min,max]',...){
  cor = cor[ct1,ct2,]
  cor  = apply(cor,1:2,mean)
  corl = apply(cor,1:2,min)
  corh = apply(cor,1:2,max)
  o1 = order(apply(cor,1,mean))
  o2 = order(apply(cor,2,mean))
  cor = cor[o1,o2]
  corl = corl[o1,o2]
  corh = corh[o1,o2]
  cort = paste0(round(cor,2),'\n[',round(corl,2),',',round(corh,2),']')
  imageWithText(cor,cort,centerColors0 = TRUE,main=main)
}

plotPairCorrHM(c2l.cor,endo,macro)
plotPairCorrHM(c2l.cor,fibro,macro)
plotPairCorrHM(c2l.cor,fibro,macro)
plotPairCorrHM(c2l.cor,neuro,macro)

pcw = meta$PCWclean
pcw2col = cols=char2col(pcw,palette = T)

z=plotPairCorr(c2l.cor,pairs=data.frame(ct1=macro,ct2=c("WNT2+ fibroblast")),pcw,pcw2col)
z
cs = colnames(c2l$WSSS_THYst9699523)

par(mar=c(15,15,1,1))

nmfbst$cl[macro]
nmfbst$cl[fibro]
# _current #####
pcw = meta$PCWclean
pcw2col = cols=char2col(pcw,palette = T)

ma2en = data.frame(ct1=rep(macro,each=length(endo)),
                   ct2=rep(endo,times=length(macro)))

epairs = list('main 3d'  = data.frame(ct1='WNT2+ fibroblast',ct2=macro[c(1,4,2,3)]),#c("LYVE1++ macrophage","TREM2+ macrophage")),
              'main 3f'  = data.frame(ct1="TREM2+ macrophage",ct2=neuro),
              'ED2'= data.frame(ct1='Pre-dermal condensate',ct2=rev(c('DC1','DC2','LTi','ILC3'))))

pdf('figures/paper.figures/1e_3d_3f_spatial_pcc.pdf',w=20,h=10)
par(mfrow=c(1,4),mar=c(3,8,1,1),tcl=-0.2,mgp=c(1.3,0.3,0),yaxs='i')
z1=plotPairCorr(c2l.cor,ma2en,pcw,pcw2col,cex.names=0.7,xlim=c(-0.4,1),xlab='Pearson correlation coefficient',ylim=c(0,260),main='main 1e')
legend(0.9,max(z1$b),fill=pcw2col,bty='n',border=NA,legend=names(pcw2col),title='PCW',xpd=NA)
mar_ = grconvertY(0:1,'user','lines')
mar_ = (mar_[2]-mar_[1])


for(figno in names(epairs)){
  n = nrow(ma2en)-nrow(epairs[[figno]])
  mar = mar_*n*(length(c2l)+2)
  par(mar=c(3+mar,8,1,1))
  z2=plotPairCorr(c2l.cor,epairs[[figno]],pcw,pcw2col,cex.names=0.7,xlim=c(-0.4,1),xlab='Pearson correlation coefficient',ylim=c(0,nrow(epairs[[figno]])*13),main=figno)
}
dev.off()

sid='WSSS_THYst9383359'
sid='WSSS_THYst9383362'
cts =c("TREM2+ macrophage","Capillaries")
plotVisiumMultyColours(vsf[[sid]],c2l[[sid]][,cts],c('red','blue'),mode = 'mean',zfun = function(x)x^2,he.grayscale=TRUE,min.opacity = 150)

# _v2 #######
# __pre-dermal condensate (fig 1e) ########
# select immune cells (most lymphoid and myeloid)
ctpairs = data.frame(ct1='Pre-dermal condensate',ct2=c('LTi','ILC3','DC2','DC1'))
ctpcor.per.sam = sapply(c2l,function(x){x=sweep(x,1,apply(x,1,sum),'/');apply(ctpairs,1,function(ct)cor(x[,ct[1]],x[,ct[2]]))})
o = order(apply(ctpcor.per.sam,1,mean))
ctpairs = ctpairs[o,]
ctpcor.per.sam = ctpcor.per.sam[o,]

cols=char2col(meta$PCD,palette = T)
ord = nrow(meta):1

pdf('figures/paper.figures/1e.Pre-dermal_condensate_сolocalization.pdf',w=5,h=5)
par(mar=c(4,3,1,1))
b=barplot(t(ctpcor.per.sam[,ord]),beside = T,border = F,col=cols[as.character(meta$PCD)[ord]],names.arg = ctpairs$ct2,las=2,
          space = c(0,5),cex.names=1,horiz = TRUE,xlab='Pearson correlation coefficient',xlim=c(-0.1,1),
          main='Colocalization with Pre-dermal condensate')
legend('bottomright',fill=cols,legend=names(cols),border=NA,bty='n',title='PCD')
# x = apply(b,2,mean)
# y = pmin(0,apply(t(ctpcor.per.sam[,ord]),2,min))
# segments(y,x,-1,x,lty=2,col='gray')
dev.off()

# __macro-endo (fig 1e (or 4?) ########
end = c('Early endothelial cells','Capillary arterioles','Arterioles','Capillaries','Postcapillary venules')
ctpairs = rbind(data.frame(ct1='Iron-recycling macrophage',ct2=end),
                data.frame(ct1='LYVE1++ macrophage',ct2=end),
                data.frame(ct1='TREM2+ macrophage',ct2=end),
                data.frame(ct1='MHCII+ macrophage',ct2=end))
ctpcor.per.sam = sapply(c2l,function(x){x=sweep(x,1,apply(x,1,sum),'/');apply(ctpairs,1,function(ct)cor(x[,ct[1]],x[,ct[2]]))})
#o = order(apply(ctpcor.per.sam,1,mean))
#ctpairs = ctpairs[o,]
#ctpcor.per.sam = ctpcor.per.sam[o,]

cols=char2col(meta$PCD,palette = T)
ord = nrow(meta):1

pdf('figures/paper.figures/1e.Macro-vs-endo.pdf',w=7,h=8)
par(mar=c(4,13,1,1))
b=barplot(t(ctpcor.per.sam[,ord]),beside = T,border = F,col=cols[as.character(meta$PCD)[ord]],names.arg = paste0(ctpairs$ct1,'\n',ctpairs$ct2),las=2,
          space = c(0,5),cex.names=0.7,horiz = TRUE,xlab='Pearson correlation coefficient',xlim=c(-0.3,1),
          main='')
legend(grconvertX(0,'ndc','user'),grconvertY(1,'ndc','user'),fill=cols,legend=names(cols),border=NA,bty='n',title='PCD',xpd=T)
x = apply(b,2,mean)
y = pmin(0,apply(t(ctpcor.per.sam[,ord]),2,min))
segments(y,x,-1,x,lty=2,col='gray')
dev.off()

# _old ####
# __v1 ####
# manually selected pairs of interest
setdiff(neu,colnames(c2l_sd0$HCA_rFSKI13460601))
c('LYVE1++ macrophage',)
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

# _1e ##################
ctpcor.per.sam = sapply(c2l,function(x){x=sweep(x,1,apply(x,1,sum),'/');apply(ctpairs,1,function(ct)cor(x[,ct[1]],x[,ct[2]]))})
cols=char2col(meta$PCD,palette = T)


c2lm = as.matrix(do.call(rbind,lapply(names(c2l),function(n){x=c2l[[n]];rownames(x)=paste0(n,'|',rownames(x));x})))
cp.cor = NULL
spotm = data.frame(sid=sapply(strsplit(rownames(c2lm),'|',T),'[',1),
                   barcode=sapply(strsplit(rownames(c2lm),'|',T),'[',2))
spotm$pcd = meta[spotm$sid,'PCD']

filters = list('6 PCW'=spotm$pcd==43,'8 PCW'=spotm$pcd==57,'10 PCW'=spotm$pcd==70,all=T)
for(norm in c(TRUE,FALSE)){
  for(f in names(filters)){
    t = c2lm
    if(norm){
      t = sweep(t,1,apply(t,1,sum),'/')
    }
    t = t[filters[[f]],]
    for(p in 1:nrow(ctpairs)){
      cr = cor.test(t[,ctpairs$ct1[p]],t[,ctpairs$ct2[p]])
      cp.cor = rbind(cp.cor,data.frame(cell1=ctpairs$ct1[p],cell2=ctpairs$ct2[p],samples=f,norm=norm,pcc=cr$estimate,cil=cr$conf.int[1],cih=cr$conf.int[2]))
    }
  }
}

cp.cor[1:2,]
n = T

castpcc = function(d,v){
  r = reshape::cast(d,pair ~ samples,value = v)
  rownames(r) = r[,'pair']
  r = r[paste0(ctpairs$ct1,'\n',ctpairs$ct2),names(filters)]
  t(as.matrix(as.data.frame(r)))
}

cp.cor$pair = paste0(cp.cor$cell1,'\n',cp.cor$cell2)
pccs = castpcc(cp.cor[cp.cor$norm==n,],'pcc')
pcc1 = castpcc(cp.cor[cp.cor$norm==n,],'cil')
pcc2 = castpcc(cp.cor[cp.cor$norm==n,],'cih')



pdf('figures/paper.figures/1e_draft.pdf',w=16,h=16)
# horiszont
ord = nrow(meta):1
ord2 = 4:1
pord = nrow(ctpairs):1
cols2=c(cols,'gray')

par(mfrow=c(1,2),mar=c(3,10,1,1),tcl=-0.2,mgp=c(1.3,0.3,0))
b=barplot(t(ctpcor.per.sam[pord,ord]),
          col=cols[as.character(meta$PCD)[ord]],
          names.arg = paste0(ctpairs$ct1,'\n',ctpairs$ct2)[pord],las=2,space = c(0,5),cex.names=0.6,
          horiz = TRUE,beside = T,border = F,xlab='Pearson correlation coefficient')
legend('topright',fill=cols,legend=paste0(round(as.numeric(names(cols))/7),' PCW'),border=NA,bty='n',title='PCD')
x = apply(b,2,mean)
segments(pmin(0,apply(ctpcor.per.sam[pord,ord],1,min)),x,-1,x,lty=2,col='gray')
#
legend.text = TRUE
b=barplot(pccs[ord2,pord],beside = T,horiz = T,xlab='Pearson correlation coefficient',
          las=2,legend.text = legend.text,col=cols2[ord2],border = NA,cex.names=0.6,
          args.legend = list(x='topright',bty='n',border=NA,title='PCD'),xlim=range(pcc1,pcc2))
segments(pcc1[ord2,pord],b,pcc2[ord2,pord],b)
x=apply(b,2,mean)
segments(pmin(0,apply(pccs[ord2,pord],2,min)),x,-1,x,lty=2,col='gray')

# verical
ord = 1:nrow(meta)
ord2 = 1:4
pord = 1:nrow(ctpairs)

par(mfrow=c(2,1),mar=c(10,3,1,1),tcl=-0.2,mgp=c(1.3,0.3,0))
b=barplot(t(ctpcor.per.sam[pord,ord]),
          col=cols[as.character(meta$PCD)[ord]],
          names.arg = paste0(ctpairs$ct1,'\n',ctpairs$ct2)[pord],las=2,space = c(0,5),cex.names=0.6,
          horiz = F,beside = T,border = F,ylab='Pearson correlation coefficient')
legend('topright',fill=cols,legend=paste0(round(as.numeric(names(cols))/7),' PCW'),border=NA,bty='n',title='PCD')
x = apply(b,2,mean)
segments(x,pmin(0,apply(ctpcor.per.sam[pord,ord],1,min)),x,-1,lty=2,col='gray')
#
legend.text = TRUE
b=barplot(pccs[ord2,pord],beside = T,horiz = F,ylab='Pearson correlation coefficient',
          las=2,legend.text = legend.text,col=cols2[ord2],border = NA,cex.names=0.6,
          args.legend = list(x='topright',bty='n',border=NA,title='PCD'),ylim=range(pcc1,pcc2))
segments(b,pcc1[ord2,pord],b,pcc2[ord2,pord])
x=apply(b,2,mean)
segments(x,pmin(0,apply(pcc1[ord2,pord],2,min)),x,-1,lty=2,col='gray')

dev.off()


# NMF #############
cols = RColorBrewer::brewer.pal(10,'Set3')
cols = setNames(cols,1:10)
colfun=function(x)num2col(x,c('yellow','orange','violet','black'),minx = 0,maxx = 1)
# _1d ###############
# not used any more
pdf('figures/paper.figures/1d_nmf.pdf',w=6,h=5)
o = order(nmfbst$cl,names(nmfbst$cl))
o = o[names(nmfbst$cl)[o] %in% unlist(ctpairs)]
par(las=1,mar=c(3,13,2,5),bty='n',tcl=-0.2,mgp=c(1.4,0.3,0))
lat= c(0.25,0.5,0.75,1)
dotPlot(t(nmfbst$coefn[,o]),grid = F,max.cex = 2.5,ylab.cex = 0.8,xlab.cex = 0.8,colfun = colfun,
        rowColours = cbind(microenv=cols[nmfbst$cl[o]]),
        colColours = cols,legend.cex.at =lat,legend.cex.title='fraction',
        xlab='Microenvironment')
dev.off()


# _ext1e ###########
pdf('figures/paper.figures/1e_extended.pdf',w=6.7,h=10)
o = order(nmfbst$cl,names(nmfbst$cl))
par(las=1,mar=c(3,13,2,5),bty='n',tcl=-0.2,mgp=c(1.4,0.3,0))
lat= c(0.25,0.5,0.75,1)
dotPlot(t(nmfbst$coefn[,o]),grid = F,max.cex = 2,ylab.cex = 0.8,xlab.cex = 0.8,colfun = colfun,
        rowColours = cbind(microenv=cols[nmfbst$cl[o]]),
        colColours = cols,legend.cex.at =lat,legend.cex.title='fraction',
        xlab='Microenvironment')
dev.off()
