
library(edgeR)

rawData = read.table('regen.sense.filter20.outlierRemoved.sense.counts', row.names = 1, skip=1)

piece = factor(c('head','head','head','head','head','head','head','head','head','head','head','head','head','head','tail','tail','tail','tail','tail','tail','tail','tail','tail','tail','tail','tail','tail','tail'))
time = factor(c('00','00','06','06','12','12','24','24','36','36','48','48','72','72','00','00','06','06','12','12','24','24','36','36','48','48','72','72'))

df = data.frame(piece,time)

grouping = factor(paste(df$piece,df$time,sep='.'))

design = model.matrix(~0+grouping)

listGroup = rep(levels(grouping),each=2)
colnames(design) = levels(grouping)

DGEData = DGEList(counts=rawData,group = listGroup)

DGEData = calcNormFactors(DGEData)

normData = cpm(DGEData,normalized.lib.sizes=FALSE)

write.table(normData,'edgeR.normalized.counts',sep='\t',quote=FALSE)

DGEData = estimateGLMCommonDisp(DGEData,design)
DGEData = estimateGLMTrendedDisp(DGEData,design)
DGEData = estimateGLMTagwiseDisp(DGEData,design)

fit = glmFit(DGEData, design)

comps = makeContrasts(head.00.tail.00 = head.00-tail.00,head.00.head.06 = head.00-head.06,head.06.tail.06 = head.06-tail.06,head.06.head.12 = head.06-head.12,head.12.head.24 = head.12-head.24,head.12.tail.12 = head.12-tail.12,head.24.tail.24 = head.24-tail.24,head.24.head.36 = head.24-head.36,head.36.tail.36 = head.36-tail.36,head.36.head.48 = head.36-head.48,head.48.head.72 = head.48-head.72,head.48.tail.48 = head.48-tail.48,head.72.tail.72 = head.72-tail.72,tail.00.tail.06 = tail.00-tail.06,tail.06.tail.12 = tail.06-tail.12,tail.12.tail.24 = tail.12-tail.24,tail.24.tail.36 = tail.24-tail.36,tail.36.tail.48 = tail.36-tail.48,tail.48.tail.72 = tail.48-tail.72,levels=design)

head.00_tail.00 = glmLRT(fit, contrast=comps[,'head.00.tail.00'])
head.00_head.06 = glmLRT(fit, contrast=comps[,'head.00.head.06'])
head.06_tail.06 = glmLRT(fit, contrast=comps[,'head.06.tail.06'])
head.06_head.12 = glmLRT(fit, contrast=comps[,'head.06.head.12'])
head.12_head.24 = glmLRT(fit, contrast=comps[,'head.12.head.24'])
head.12_tail.12 = glmLRT(fit, contrast=comps[,'head.12.tail.12'])
head.24_tail.24 = glmLRT(fit, contrast=comps[,'head.24.tail.24'])
head.24_head.36 = glmLRT(fit, contrast=comps[,'head.24.head.36'])
head.36_tail.36 = glmLRT(fit, contrast=comps[,'head.36.tail.36'])
head.36_head.48 = glmLRT(fit, contrast=comps[,'head.36.head.48'])
head.48_head.72 = glmLRT(fit, contrast=comps[,'head.48.head.72'])
head.48_tail.48 = glmLRT(fit, contrast=comps[,'head.48.tail.48'])
head.72_tail.72 = glmLRT(fit, contrast=comps[,'head.72.tail.72'])
tail.00_tail.06 = glmLRT(fit, contrast=comps[,'tail.00.tail.06'])
tail.06_tail.12 = glmLRT(fit, contrast=comps[,'tail.06.tail.12'])
tail.12_tail.24 = glmLRT(fit, contrast=comps[,'tail.12.tail.24'])
tail.24_tail.36 = glmLRT(fit, contrast=comps[,'tail.24.tail.36'])
tail.36_tail.48 = glmLRT(fit, contrast=comps[,'tail.36.tail.48'])
tail.48_tail.72 = glmLRT(fit, contrast=comps[,'tail.48.tail.72'])
write.table(topTags(head.00_tail.00,n=18000,adjust.method = 'BH', sort.by='p.value'),'head.00.tail.00.de',sep='\t',quote=FALSE)
write.table(topTags(head.00_head.06,n=18000,adjust.method = 'BH', sort.by='p.value'),'head.00.head.06.de',sep='\t',quote=FALSE)
write.table(topTags(head.06_tail.06,n=18000,adjust.method = 'BH', sort.by='p.value'),'head.06.tail.06.de',sep='\t',quote=FALSE)
write.table(topTags(head.06_head.12,n=18000,adjust.method = 'BH', sort.by='p.value'),'head.06.head.12.de',sep='\t',quote=FALSE)
write.table(topTags(head.12_head.24,n=18000,adjust.method = 'BH', sort.by='p.value'),'head.12.head.24.de',sep='\t',quote=FALSE)
write.table(topTags(head.12_tail.12,n=18000,adjust.method = 'BH', sort.by='p.value'),'head.12.tail.12.de',sep='\t',quote=FALSE)
write.table(topTags(head.24_tail.24,n=18000,adjust.method = 'BH', sort.by='p.value'),'head.24.tail.24.de',sep='\t',quote=FALSE)
write.table(topTags(head.24_head.36,n=18000,adjust.method = 'BH', sort.by='p.value'),'head.24.head.36.de',sep='\t',quote=FALSE)
write.table(topTags(head.36_tail.36,n=18000,adjust.method = 'BH', sort.by='p.value'),'head.36.tail.36.de',sep='\t',quote=FALSE)
write.table(topTags(head.36_head.48,n=18000,adjust.method = 'BH', sort.by='p.value'),'head.36.head.48.de',sep='\t',quote=FALSE)
write.table(topTags(head.48_head.72,n=18000,adjust.method = 'BH', sort.by='p.value'),'head.48.head.72.de',sep='\t',quote=FALSE)
write.table(topTags(head.48_tail.48,n=18000,adjust.method = 'BH', sort.by='p.value'),'head.48.tail.48.de',sep='\t',quote=FALSE)
write.table(topTags(head.72_tail.72,n=18000,adjust.method = 'BH', sort.by='p.value'),'head.72.tail.72.de',sep='\t',quote=FALSE)
write.table(topTags(tail.00_tail.06,n=18000,adjust.method = 'BH', sort.by='p.value'),'tail.00.tail.06.de',sep='\t',quote=FALSE)
write.table(topTags(tail.06_tail.12,n=18000,adjust.method = 'BH', sort.by='p.value'),'tail.06.tail.12.de',sep='\t',quote=FALSE)
write.table(topTags(tail.12_tail.24,n=18000,adjust.method = 'BH', sort.by='p.value'),'tail.12.tail.24.de',sep='\t',quote=FALSE)
write.table(topTags(tail.24_tail.36,n=18000,adjust.method = 'BH', sort.by='p.value'),'tail.24.tail.36.de',sep='\t',quote=FALSE)
write.table(topTags(tail.36_tail.48,n=18000,adjust.method = 'BH', sort.by='p.value'),'tail.36.tail.48.de',sep='\t',quote=FALSE)
write.table(topTags(tail.48_tail.72,n=18000,adjust.method = 'BH', sort.by='p.value'),'tail.48.tail.72.de',sep='\t',quote=FALSE)
