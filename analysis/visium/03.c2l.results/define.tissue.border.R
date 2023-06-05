library(Seurat)
library(visutils) # devtools::install_github("iaaka/visutils")

meta = read.csv('src/my/visium/samples.csv')

# limb ######
# here skin regions were manually subsetted, so I'll use whole samples to define tissue border
vs1 = lapply(meta$Sanger_id[meta$site=='limb'],function(s)myLoad10X_Spatial(paste0('data.lustre/visium/spaceranger/',s,'/outs'),filter.matrix = FALSE))
names(vs1) = meta$Sanger_id[meta$site=='limb']

ans = lapply(names(vs1),function(f)read.csv(paste0('data.nfs/my/visium/man.ann/',f,'.csv')))
names(ans) = names(vs1)


for(n in names(ans)){
  a = ans[[n]]
  a = colnames(vs1[[n]]) %in% a$Barcode[a$annot=='skin']
  vs1[[n]]$is.skin = a + 0
  vs1[[n]]@images$slice1@coordinates$tissue = vs1[[n]]$is.tissue = pmax(vs1[[n]]$is.tissue,vs1[[n]]$is.skin)
}

# manual examination shown that there are some 10 misannotated non-tissue spots in sample WSSS_THYst9383359 that preclude correct determination of tissue border
# lets remove them based on H&E brightness
n = 'WSSS_THYst9383359'
sf = jsonlite::read_json(paste0('data.lustre/visium/spaceranger/',n,'/outs/spatial/scalefactors_json.json'))
br = visutils::getMeanSpotColor(vs1[[n]],sf)
br = apply(br,1,mean)
hist(br[vs1[[n]]$is.tissue==0],0:100/100,freq = F,border = F,xlim=c(0.6,1))
hist(br[vs1[[n]]$is.tissue==1],0:100/100,freq = F,add=T,border = F,col='#FF000090')

# brigtness only is not enaugh to select these spots
toremove =  br>0.9 & vs1[[n]]$is.skin==0 & vs1[[n]]@images$slice1@coordinates$imagerow*vs1[[n]]@images$slice1@scale.factors$lowres < 200
plotVisium(vs1[[n]],as.character(toremove),spot.filter = vs1[[n]]$is.tissue==1,cex=0.5,
           xaxt='s',yaxt='s',xlim=c(380,530),ylim=c(260,440),border = vs1[[n]]$is.skin==1)
vs1[[n]]@images$slice1@coordinates$tissue[toremove] = vs1[[n]]$is.tissue[toremove] = 0

for(n in names(ans)){
  brd = findTissueBorder(vs1[[n]])
  dist = calcDistance2border(brd$rc,brd$nj)
  vs1[[n]]$is.border = dist$is.border
  vs1[[n]]$dist2border.graph = dist$dist2border.graph
}

par(mfrow=c(2,4),mar=c(0,0,1,1))
for(n in names(vs1)){
  plotVisium(vs1[[n]],as.character(vs1[[n]]$is.border),spot.filter = vs1[[n]]$is.skin==1)
}


par(mfrow=c(2,4),mar=c(0,0,1,1))
for(n in names(vs1)){
  plotVisium(vs1[[n]],vs1[[n]]$dist2border.graph,spot.filter = vs1[[n]]$is.skin==1,z2col = function(x)num2col(x,col=c('blue','red')))
}

# body-n-face #####
# here all spots are skin, however some manual filtering plus border selection was performed
vs2 = lapply(meta$Sanger_id[meta$site!='limb'],function(s)myLoad10X_Spatial(paste0('data.lustre/visium/spaceranger/',s,'/outs'),filter.matrix = FALSE))
names(vs2) = meta$Sanger_id[meta$site!='limb']

ans2 = lapply(meta$Sanger_id[meta$site!='limb'],function(f)read.csv(paste0('data.nfs/my/visium/man.ann/',f,'.csv')))
names(ans2) = meta$Sanger_id[meta$site!='limb']


for(n in names(ans2)){
  a = ans2[[n]]
  f = colnames(vs2[[n]]) %in% a$Barcode[a$ann=='not_tissue']
  vs2[[n]]$tmp = ifelse(vs2[[n]]$is.tissue==0,'not_tissue','tissue')
  vs2[[n]]@images$slice1@coordinates$tissue[f] = vs2[[n]]$is.tissue[f] = 0
  vs2[[n]]$is.border = colnames(vs2[[n]]) %in% a$Barcode[a$ann=='border']
  vs2[[n]]$tmp[vs2[[n]]$is.border] = 'border'
  vs2[[n]]$tmp[vs2[[n]]$tmp!='not_tissue' & vs2[[n]]$is.tissue==0] = 'to_remove'
}

# vs2_ = vs2
# for(n in names(vs2))
#   vs2_[[n]]@images$slice1@image = enhanceImage(vs2[[n]]@images$slice1@image,qs = c(0.001,0.2))
# 
# pdf('figures/visium/face.body/clean.pdf',w=3.5*4,h=6)
# cols = setNames(RColorBrewer::brewer.pal(3,'Set1'),c('to_remove','border','tissue'))
# par(mfcol=c(2,4),mar=c(0,0,1,6),bty='n')
# for(n in names(vs2)){
#   plotVisium(vs2_[[n]],vs2[[n]]$tmp,cex=0,spot.filter = vs2[[n]]$tmp!='not_tissue',main=n,plot.legend = F)
#   plotVisium(vs2_[[n]],vs2[[n]]$tmp,cex=0.6,z2col=cols,spot.filter = vs2[[n]]$tmp!='not_tissue',main=n)
# }
# dev.off()


for(n in names(vs2)){
  vs2[[n]]$tmp = NULL
  brd = findTissueBorder(vs2[[n]])
  brd$rc$is.border = vs2[[n]]$is.border == 1
  brd$rc$border.inx[vs2[[n]]$is.border != 1] = NA
  dist = calcDistance2border(brd$rc,brd$nj)
  vs2[[n]]$is.border = dist$is.border
  vs2[[n]]$dist2border.graph = dist$dist2border.graph
  vs2[[n]]$is.skin = vs2[[n]]$is.tissue
}


vs = c(vs1,vs2)
vsf = lapply(vs,function(v)v[,v$is.skin==1])

zlim = range(unlist(lapply(vsf,function(v)v$dist2border.graph)))

pdf('figures/visium/distance2border.pdf',w=3.5*4,h=9)
par(mfrow=c(3,4),mar=c(0,0,1,6),bty='n')
for(n in names(vsf)){
  plotVisium(vsf[[n]],vsf[[n]]$dist2border.graph,z2col = function(x)num2col(x,col=c('red','orange','yellow','green','blue','violet','black')),
             he.grayscale=TRUE,main=n,border = ifelse(vsf[[n]]$is.border,'black',NA),zlim = zlim,legend.args = list(title='Dist to border\n(spots)'))
}
dev.off()

saveRDS(vsf,'data.nfs/my/visium/viss.all.filtered.rds')

