# file names should conform 10x requrements:
# [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
fqpath='data.lustre/visium/fq/fetal.face.body/'
files = list.files(fqpath,pattern = '*fastq.gz',recursive = TRUE)
dirs = dirname(files)
basenames = basename(files)

suffix = sapply(strsplit(basenames,'_'),function(x){paste(x[(length(x)-2):length(x)],collapse = '_')})
newbasenames = paste0(dirs,'_S1_',suffix)

from = paste0(fqpath,'/',dirs,'/',basenames)
to = paste0(fqpath,'/',dirs,'/',newbasenames)

file.rename(from,to)

samples = rbind(data.frame(dir='limb',sample=list.dirs('data.lustre/visium/fq/limb',recursive = F,full.names = F)),
                data.frame(dir='fetal.face.body',sample=list.dirs('data.lustre/visium/fq/fetal.face.body',recursive = F,full.names = F)))
write.table(samples,row.names = F,col.names = F,quote = F,'data.lustre/visium/spaceranger/samples.txt')
