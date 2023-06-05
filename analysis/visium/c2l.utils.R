plotCellAbundanceOnDist2Border = function(d2b,x,log=FALSE,cols=c('#00FFFF','#FF00FF','#FFFF00'),...){
  cilim = c(0,Inf)
  if(log){
    x = log10(x)
    cilim = NULL
  }
  
  areas = lapply(1:ncol(x),function(i){
    r = t(sapply(split(x[,i],d2b),function(y)c(mean(y),sd(y)/sqrt(length(y)),length(y))))
    r[r[,3]>2,1:2]
  })
  
  
  t = do.call(rbind,areas)
  if(log)
    ylim = range(t[,1]-2*t[,2],t[,1]+2*t[,2])
  else
    ylim = c(0,max(t[,1]+2*t[,2]))
  
  
  plotArea(as.numeric(rownames(areas[[1]])),areas[[1]],col=cols[1],sd.mult = 2,lwd=3,cilim = cilim,
           new = T,ylim=ylim,xlab='distance to tissue border (spots)',ylab='cell abundance',...)
  if(length(areas)>1)
    for(i in 2:length(areas))
      plotArea(as.numeric(rownames(areas[[i]])),areas[[i]],col=cols[i],new = F,sd.mult = 2,lwd=2,cilim = cilim)  
}


#' @param d path to folder (predmodel) with csv outputs of c2l
#' @param sd2zero values < sd2zero*sd will be set to zero. 
loadC2L = function(d,sd2zero=0){
  r = read.csv(paste0(d,'/q05_cell_abundance_w_sf.csv'),row.names = 1,check.names = FALSE)
  sd = read.csv(paste0(d,'/stds_cell_abundance_w_sf.csv'),row.names = 1,check.names = FALSE)
  colnames(r) = sub('q05cell_abundance_w_sf_','',colnames(r))
  colnames(sd) = sub('stdscell_abundance_w_sf_','',colnames(sd))
  sd = sd [rownames(r),colnames(r)]
  r[r<sd*sd2zero] = 0
  mm = do.call(rbind,strsplit(rownames(r),'|',T))
  r$barcode = mm[,2]
  r = split(r,mm[,1])
  
  for(i in 1:length(r)){
    rownames(r[[i]]) = r[[i]]$barcode
    r[[i]]$barcode = NULL
    r[[i]] = as.matrix(r[[i]])
  }
  r
}

getBestMatchByContigencyTable = function(t,useHungarian=FALSE){
  if(useHungarian){
    o = as.data.frame(RcppHungarian::HungarianSolver(max(t)-t)$pairs)
    colnames(o) = c('row','col')
    return(o[o$row>0 & o$col>0,])
  }
  o = order(t,decreasing = T) -1
  o = data.frame(row=o %% nrow(t) + 1,col=o %/% nrow(t)+1)
  o$good = FALSE
  r = NULL
  rs = 1:nrow(t)
  cs = 1:ncol(t)
  for(i in 1:nrow(o)){
    if(all(is.na(rs)) || all(is.na(cs))) break
    if(!is.na(rs[o$row[i]]) & !is.na(cs[o$col[i]])){
      cs[o$col[i]] = NA
      rs[o$row[i]] = NA
      o$good[i] = TRUE
    }
  }
  o[o$good,1:2]
}


getCoincedence = function(cs){
  c = cs[[1]]
  r = matrix(0,ncol=length(c),nrow=length(c),dimnames=list(names(c),names(c)))
  for(c in cs){
    for(i in 1:length(c))
      for(j in 1:length(c))
        r[i,j] = r[i,j] + (c[i]==c[j])
  }
  r
}


nmfNormFs =list(
  no  = function(n){rep(1,ncol(basis(n)))},
  max = function(n){apply(basis(n),2,max)}
) 


getNMFcl = function(ns,normF){
  lapply(ns,function(n){
    c = coef(n)
    f = normF(n)
    c = sweep(c,1,f,'*')
    apply(c,2,which.max)
  })
}

nmf2cl = function(n,normF){
  c = coef(n)
  b = basis(n)
  f = normF(n)
  c = sweep(c,1,f,'*')
  b = sweep(b,2,f,'/')
  cl = apply(c,2,which.max)
  c = sweep(c,2,apply(c,2,sum),'/')
  list(nmf=n,coefn = c,basisn=b,cl=cl)
}


getSilhouetteStat = function(cls,d){
  do.call(rbind,lapply(cls,function(cl){
    x = cluster::silhouette(cl,d)[,3]
    data.frame(min=min(x),mean=mean(x),median=median(x),n0 = sum(x<0))
  }))
}


nmfGetBest = function(ns,normF){
  cls = getNMFcl(ns,normF)
  cns = getCoincedence(cls)
  slh = getSilhouetteStat(cls,1-cns/100)
  slh$dev = sapply(ns,deviance)
  i = order(-slh$mean,slh$dev)[1]
  n = ns[[i]]
  r = nmf2cl(n,normF)
  r$cons=cns
  r$slh.stat=slh
  r$besti=i
  r
}

reorderContigencyTable = function(t){
  t = t[order(apply(t,1,sum),decreasing = T),]
  tn = sweep(t,1,apply(t,1,sum),'/')
  t[,order(apply(tn,2,function(x)which.max(x)))]
}