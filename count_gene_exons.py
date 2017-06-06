import HTSeq
import glob
from itertools import *
import collections
import time
import pandas

files=glob.glob("/data/yosef2/users/chenling/YosefCode/packages/UMI/Nir_David/processed/*.dedup.bam")
# gtf_file = HTSeq.GFF_Reader( "/data/yosef2/users/chenling/YosefCode/packages/UMI/Mus_musculus.GRCm38.82.chr.gtf",end_included=True )
gtf_file=HTSeq.GFF_Reader("/data/yosef2/users/chenling/ref/Homo_sapiens.GRCh38.88.chr.gtf")
exons = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
genes = HTSeq.GenomicArrayOfSets( "auto", stranded=True )

for feature in gtf_file:
    if feature.type == "exon":
       exons[ feature.iv ] += feature.attr["gene_name"]
    if feature.type == "gene":
       genes[ feature.iv ] += feature.attr["gene_name"]



barcode=[]
for filename in files:
	# filename = filename.split("_")[1]
	filename = filename.split("/")[10]
	barcode.append(filename.split(".")[0])

counts_exon =[]
counts_gene =[]
counts_filename =[]

start = time.clock()
for file in files:
	almnt_file = HTSeq.BAM_Reader( file )
	exon_count = collections.Counter( )
	gene_count = collections.Counter( )
	for almnt in almnt_file:
		if not almnt.aligned:
			count[ "_unmapped" ] += 1
			continue
		exonid = set()
		geneid = set()
		for cigop in almnt.cigar:
			if cigop.type != "M":
			  continue
			for iv, val in exons[ cigop.ref_iv ].steps():
			  exonid |= val
			for iv, val in genes[ cigop.ref_iv ].steps():
			  geneid |= val
		if len(exonid) == 1:
			exonid = list(exonid)[0]
			exon_count[ exonid ] += 1
		elif len(exonid) == 0:
			exon_count[ "_no_feature" ] += 1
		else:
			exon_count[ "_ambiguous" ] += 1
		if len(geneid) == 1:
			geneid = list(geneid)[0]
			gene_count[ geneid ] += 1
		elif len(geneid) == 0:
			gene_count[ "_no_feature" ] += 1
		else:
			gene_count[ "_ambiguous" ] += 1
	counts_exon.append(exon_count)
	counts_gene.append(gene_count)
	counts_filename.append(file)


print time.clock()-start

counts_barcode=[x.split('/')[10].split('.')[0].split('_')[1] for x in counts_filename]
# exon=pandas.DataFrame(counts_exon,index=counts_barcode)
# gene=pandas.DataFrame(counts_gene,index=counts_barcode)
exon=pandas.DataFrame(counts_exon,index=barcode)
gene=pandas.DataFrame(counts_gene,index=barcode)
exon=exon.transpose()
gene=gene.transpose()
exon.to_csv('counts_exon.csv',index=True)
gene.to_csv('counts_gene.csv',index=True)


################################# R ##############################
################################# R ##############################
################################# R ##############################
