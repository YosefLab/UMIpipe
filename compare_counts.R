count_files=list.files()
count_files=count_files[grep('.exon.counts.tsv',count_files)]

dge=read.table('/data/yosef2/users/chenling/YosefCode/packages/UMI/Macosko2015/SRR1853178.star_gene_exon_tagged.dge.txt',header=T,row.names=1,as.is=T)
prefix=sapply(strsplit(count_files,".",fixed=T),function(X){X[1]})

ht_gene=lapply(prefix,function(X){
	gene=read.table(paste(X,'.gene.counts.tsv',sep=''), row.names = 1, header= TRUE,as.is=T)
	return(gene)
})
ht_exon=lapply(prefix,function(X){
	exon=read.table(paste(X,'.exon.counts.tsv',sep=''), row.names = 1, header= TRUE,as.is=T)
	return(exon)
})
qc=lapply(prefix,function(X){
	qc=read.table(paste(X,'.deduplicated_edit_distance.tsv',sep=''), header= TRUE,as.is=T)
	return(qc)
})

shared=intersect(rownames(ht_gene[[1]]),rownames(dge))
ht_gene=lapply(ht_gene,function(X){X[match(shared,rownames(X)),1]})
ht_exon=lapply(ht_exon,function(X){X[match(shared,rownames(X)),1]})
ht_gene=do.call(cbind,ht_gene)
ht_exon=do.call(cbind,ht_exon)

dge=dge[match(shared,rownames(dge)),]

ht_names=sapply(strsplit(prefix,'_',fixed=T),function(X){X[2]})
cor=sapply(c(1:length(ht_exon[1,])),function(i){
	cor(ht_exon[,i],dge[,colnames(dge)==ht_names[i]])
})

write.table(cbind(ht_names,cor),file='ht_count.compare.tsv',quote=F,sep="\t",row.names=F,col.names=F)

qc=lapply(prefix,function(X){
	qc=read.table(paste(X,'.deduplicated_edit_distance.tsv',sep=''), header= TRUE,as.is=T)
	return(qc)
})

diff=sapply(qc,function(X){sum(X[c(2:11),5]-X[c(2:11),1])/sum(X[c(2:11),5])})
total_reads=sapply(qc,function(X){sum(X[c(2:11),5])})
total_UMI=sapply(qc,function(X){sum(X[c(2:11),1])})

cor(cor,diff,method='spearman')
cor(cor,total_reads,method='spearman')
cor(cor,total_UMI,method='spearman')

exon=read.csv('processed/counts_exon.csv',row.names=1,as.is=T)
exon=t(exon)
dge=read.table('/data/yosef2/users/chenling/YosefCode/packages/UMI/Macosko2015/SRR1853178.star_gene_exon_tagged.dge.txt',header=T,row.names=1,as.is=T)
shared=intersect(rownames(exon),rownames(dge))
dge=dge[match(shared,rownames(dge)),]
exon=exon[match(shared,rownames(exon)),]
cells=intersect(colnames(exon),colnames(dge))
dge=dge[,match(cells,colnames(dge))]
exon=exon[,match(cells,colnames(exon))]
cor=sapply(c(1:length(exon[1,])),function(i){
	exon[is.na(exon[,i]),i]=0
	cor(exon[,i],dge[,i])
})


gene=read.csv('processed/counts_gene.csv',row.names=1,as.is=T)
gene=t(gene)
shared=intersect(rownames(gene),rownames(dge))
dge=dge[match(shared,rownames(dge)),]
gene=gene[match(shared,rownames(gene)),]
cells=intersect(colnames(gene),colnames(dge))
dge=dge[,match(cells,colnames(dge))]
gene=gene[,match(cells,colnames(gene))]
cor=sapply(c(1:length(gene[1,])),function(i){
	gene[is.na(gene[,i]),i]=0
	cor(gene[,i],dge[,i])
})
