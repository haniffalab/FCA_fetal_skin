face_body = read.csv('data.nfs/my/visium/face.body/FS_pilot_MV_metadata - Visium tracker.csv')
face_body = face_body[,c('sanger.ID','gestation.age','PCW','Tissue.Type')]
face_body$PCD = 7*as.numeric(sub('pcw','',face_body$PCW))
face_body$Tissue.Type = sapply(strsplit(face_body$Tissue.Type,' '),'[',2)

colnames(face_body) = c("Sanger_id","GA","PCW","site","PCD")

limb = openxlsx::read.xlsx('data.nfs/my/visium/201015limb_samples.xlsx')
limb = limb[29:36,c('Sample.ID','Sample.stage','Norm..Sample.Stage')]
limb$site = 'limb'
limb$Sample.stage = gsub('GA| ','',limb$Sample.stage)
limb$Norm..Sample.Stage = gsub('||Unknown','',limb$Norm..Sample.Stage,fixed = T)
# use PCW if avaliabe
limb$PCD = round(7*as.numeric(sub('pcw','',limb$Norm..Sample.Stage)))

f = !is.na(limb$Norm..Sample.Stage)
limb$Norm..Sample.Stage[f] = paste0(sub('pcw','',limb$Norm..Sample.Stage),'pcw')[f]
# all remaining have GA=10w1, set them manuall
limb$PCD[is.na(limb$PCD)] = 10*7+1 - 2*7
colnames(limb) = c("Sanger_id","GA","PCW","site","PCD")
meta = rbind(limb,face_body)
#  no skin in this sample, so I exclude it
meta = meta[meta$Sanger_id!='WSSS_THYst9383360',]
write.csv(meta,'src/my/visium/samples.csv',row.names = FALSE)
