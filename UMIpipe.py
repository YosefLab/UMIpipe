#!/usr/bin/python
import argparse
import sys
import HTSeq
import subprocess
from UMIpipe_funcs import *

parser = argparse.ArgumentParser()
parser.add_argument("--fq1", help="file name of the first read",type=str)
parser.add_argument("--fq2", help="file name of the second read",type=str)
parser.add_argument("--ref", help="name of organism and reference version",type=str)
parser.add_argument("--picard_path", help = "path to picard tools bin", default="/opt/pkg/picard-tools-2.5.0/bin/" ,type=str)
parser.add_argument("--dropseq_jar", help = "path to dropseq tools bin", default="/data/yosef2/users/chenling/YosefCode/packages/UMI/Drop-seq_tools-1.12/jar/dropseq.jar" ,type=str)
parser.add_argument("--star_path", help = "path to star bin", default="/opt/pkg/STAR-2.5.3a/bin/Linux_x86_64_static/STAR" ,type=str)
parser.add_argument("--samplename",help = "the name of the output",type=str)
parser.add_argument("--cstart",help = "start of cell barcode",default=1,type=int)
parser.add_argument("--cend",help = "end of cell barcode",default=12,type=int)
parser.add_argument("--mstart",help = "start of molecular barcode",default=13,type=int)
parser.add_argument("--mend",help = "end of molecular barcode",default=20,type=int)
parser.add_argument("--tag_read",help="1 or 2, the read that contains the barcode",default=1,type=int)
parser.add_argument("--tag_quals",help = "minimum quality score for the tag basepairs",default=10,type=int)
parser.add_argument("--adapter",help = "adapter sequence for trimming",default='AAGCAGTGGTATCAACGCAGAGTGAATGGG',type=str)
parser.add_argument("--ncells",help = "estimated number of cells in sample",type=int)
# parser.add_argument("--secondary",help = "T for keeping secondary mapping reads, F for not keeping them",type=str)

args = parser.parse_args()

if args.ref=='mm10':
	fasta="/data/yosef/index_files/mm10/Gencode/GRCm38.primary_assembly.genome.fa"
	star_index="/data/yosef/index_files/mm10/Gencode/STAR_index"
	gtf="/data/yosef/index_files/mm10/Gencode/gencode.vM10.annotation.gtf"
elif args.ref=="hg38":
	fasta="/data/yosef2/users/chenling/ref/STARindex/hg38/hg38.fa"
	star_index="/data/yosef2/users/chenling/ref/STARindex/hg38/"
	gtf="/data/yosef2/users/chenling/ref/Homo_sapiens.GRCh38.88.chr.gtf" 

#build picard dict
/opt/genomics/bin/picard CreateSequenceDictionary \
R=/data/yosef2/users/chenling/ref/STARindex/hg38/hg38.fa \
O=/data/yosef2/users/chenling/ref/STARindex/hg38/hg38.dict

# #######################################################
# #convert fastq to sam
# #######################################################
print 'convert '+args.fq1+' and '+args.fq2+' sam file'
returnCode=Fq2Sam(args.picard_path,args.fq1,args.fq2,args.samplename)

# #######################################################
# #tag barcode
# #######################################################
print 'tag '+args.samplename+' with cell barcode'
returnCode=TagSam("Ctag",args.dropseq_jar,args.samplename,args.cstart,args.cend,args.tag_quals,args.tag_read,'False','XC')
args.samplename = args.samplename+".Ctag"

# command = "/opt/genomics/bin/samtools view -h "+args.samplename+".bam |  paste -d';' - - |grep -v 'XC:Z:NNNNNNNN'" 
# command = command + '|gawk \'{if (length($10)==16)print $0}\' '
# command = command + '|gawk \'BEGIN { FS = ";" };{print $1"\\n"$2}\' >' +  args.samplename + '.filter.sam'
# returnCode = subprocess.call(command,shell=True)
# print command
# command = "/opt/genomics/bin/samtools view -H "+args.samplename+".bam > header.txt"
# returnCode = subprocess.call(command,shell=True)
# command = 'cat header.txt ' + args.samplename + '.filter.sam' + '| /opt/genomics/bin/samtools view -b >'   +  args.samplename + '.filter.bam'
# returnCode = subprocess.call(command,shell=True)
# if(returnCode != 0):
#       raise Exception("filtering bad reads");
# args.samplename = args.samplename + '.filter'

print 'tag '+args.samplename+' with molecular barcode'
#some previous versions have 'MC' for the molecular barcode
# correct with this 
# samtools view -h $f.bam|sed 's/MC:Z:/XM:Z:/' |samtools view -b > $f.corrected.bam
returnCode=TagSam("Mtag",args.dropseq_jar,args.samplename,args.mstart,args.mend,args.tag_quals,args.tag_read,'True','XM')
args.samplename = args.samplename + '.Mtag'
# #######################################################
# #filtering and trimming
# #######################################################
returnCode = FilterSam(args.samplename,args.dropseq_jar)
args.samplename = args.samplename + '.filter'
returnCode = TrimSam(args.samplename,args.dropseq_jar,args.adapter)
args.samplename = args.samplename + '.trimmed'
returnCode = TrimPolyA(args.samplename,args.dropseq_jar)
args.samplename = args.samplename + '.polyA'
unmapped = args.samplename
# #######################################################
# #Mapping 
# #######################################################
returnCode = Bam2Fastq(args.samplename,args.picard_path)
returnCode = StarMap(args.samplename,args.star_path,star_index)
returnCode = SortBam(args.samplename,args.picard_path)
args.samplename = args.samplename + '.aligned.sorted'
# #######################################################
# #Merge bams 
# #######################################################
returnCode=MergeBams(args.samplename,unmapped,args.picard_path,"T",fasta)
returnCode=MergeBams(args.samplename,unmapped,args.picard_path,"F",fasta)
withsec = args.samplename+".sec.merged"
nosec = args.samplename+".nosec.merged"
# #######################################################
# #Get exrpession 1 
# #######################################################
returnCode=TagExons(nosec,args.dropseq_jar,args.gtf)
nosec=nosec+'.exon'
returnCode=CheckMixing(nosec,args.dropseq_jar,args.ncells)
nosec=nosec+'.cleaned'
returnCode=GetExpression(nosec,args.dropseq_jar,args.ncells)
returnCode=GetNCells(nosec,args.dropseq_jar)