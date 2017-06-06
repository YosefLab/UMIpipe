#!/usr/bin/python
import argparse
import sys
import HTSeq
import subprocess

def Fq2Sam(picard_path,fq1name,fq2name,samplename):
	command=picard_path+'FastqToSam'+ ' FASTQ='+fq1name + ' FASTQ2='+ fq2name
	command = command + ' OUTPUT=' + samplename+".bam"
	command = command + ' SM=' + samplename
	command = command + ' READ_GROUP_NAME=NA' + " LIBRARY_NAME=UNKOWN" 
	command = command + " PLATFORM_UNIT=NA" + " PLATFORM=illumina" + " SEQUENCING_CENTER=NA" 

	print(command)
	returnCode = subprocess.call(command,shell=True)
	if(returnCode != 0):
		raise Exception("converting fq to sam failed")
	return returnCode

def TagSam(taskname,dropseq_jar,samplename,cstart,cend,tag_quals,tag_read,discard,tag_name):
	command='java -jar '+ dropseq_jar +' TagBamWithReadSequenceExtended' + ' INPUT='+samplename+'.bam'
	command = command + ' OUTPUT=' + samplename+"."+taskname+".bam"
	command = command + ' SUMMARY=' + samplename+"."+taskname+".summary"
	command = command + ' BASE_RANGE=' + str(cstart)+'-'+str(cend)
	command = command + ' BASE_QUALITY=' + str(tag_quals)
	command = command + ' BARCODED_READ=' + str(tag_read) 
	command = command + ' DISCARD_READ=' + discard + ' TAG_NAME=' + tag_name + ' NUM_BASES_BELOW_QUALITY=1'
	print(command)
	returnCode = subprocess.call(command,shell=True)
	if(returnCode != 0):
		raise Exception("tagging barcode failed");
	return returnCode

def FilterSam(samplename,dropseq_jar):
	command = 'java -jar '+dropseq_jar+' FilterBAM ' + 'TAG_REJECT=XQ '
	command = command + " INPUT=" + samplename+".bam" + " OUTPUT=" + samplename+".filter.bam"
	print(command)
	returnCode = subprocess.call(command,shell=True)
	if(returnCode != 0):
		raise Exception("filter sam failed");
	return returnCode

def TrimSam(samplename,dropseq_jar,adapter):
	command = 'java -jar '+dropseq_jar+' TrimStartingSequence '
	command = command + " INPUT=" + samplename+".bam" + " OUTPUT=" + samplename+".trimmed.bam"
	command = command + " OUTPUT_SUMMARY=" + samplename+".trim.summary" 
	command = command + " SEQUENCE=" + adapter + " MISMATCHES=0 NUM_BASES=5" 
	print(command)
	returnCode = subprocess.call(command,shell=True)
	if(returnCode != 0):
		raise Exception("trim sam failed");
	return returnCode

def TrimPolyA(samplename,dropseq_jar):
	command = 'java -jar '+dropseq_jar+' PolyATrimmer '
	command = command + " INPUT=" + samplename+".bam" + " OUTPUT=" + samplename+".polyA.bam"
	command = command + " OUTPUT_SUMMARY=" + samplename+".polyA.summary" 
	command = command + " MISMATCHES=0 NUM_BASES=6" 
	print(command)
	returnCode = subprocess.call(command,shell=True)
	if(returnCode != 0):
		raise Exception("trim polyA failed");
	return returnCode

def Bam2Fastq(samplename,picard_path):
	command=picard_path+'SamToFastq'
	command = command + ' INPUT=' + samplename + ".bam"
	command = command + ' FASTQ=' + samplename + ".fastq"
	print(command)
	returnCode = subprocess.call(command,shell=True)
	if(returnCode != 0):
		raise Exception("converting sam to faq failed")
	return returnCode

def StarMap(samplename,star_path,star_index):
	command=star_path + ' --genomeDir ' +  star_index
	command = command + ' --readFilesIn ' + samplename + ".fastq"
	command = command + ' --outFileNamePrefix ' + samplename + "."
	print(command)
	returnCode = subprocess.call(command,shell=True)
	if(returnCode != 0):
		raise Exception("STAR mapping failed")
	return returnCode

def SortBam(samplename,picard_path):
	command=picard_path+'SortSam'
	command = command + ' INPUT=' + samplename + ".Aligned.out.sam"
	command = command + ' OUTPUT=' + samplename + ".aligned.sorted.bam" + ' SO=queryname'
	print(command)
	returnCode = subprocess.call(command,shell=True)
	if(returnCode != 0):
		raise Exception("soring aligned sam failed")
	return returnCode

def MergeBams(samplename,unmapped,picard_path,secondary,fasta):
	command=picard_path+'MergeBamAlignment'
	command = command + ' UNMAPPED_BAM=' + unmapped + ".bam"
	command = command + ' ALIGNED_BAM=' + samplename + ".bam"
	command = command + ' R=' + fasta
	if secondary=='T':
		command = command + ' INCLUDE_SECONDARY_ALIGNMENTS=true'
		secondary='sec'
	elif secondary=='F':
		command = command + ' INCLUDE_SECONDARY_ALIGNMENTS=false'
		secondary='nosec'
	else: 
		raise Exception("--secondary only takes T or F as input")
	command = command + ' OUTPUT=' + samplename + "."+secondary+".merged.bam"
	command = command + ' PAIRED_RUN=false'
	print(command)
	returnCode = subprocess.call(command,shell=True)
	if(returnCode != 0):
		raise Exception("merging unmapped and mapped bam failed")
	return returnCode

def CheckMixing(samplename,dropseq_jar,ncells):
	command = 'java -jar '+dropseq_jar+' DetectBeadSynthesisErrors '
	command = command + " INPUT=" + samplename+".bam" + " OUTPUT=" + samplename+".cleaned.bam"
	command = command + " OUTPUT_STATS=" + samplename+".synthesis_stats.txt" 
	command = command + " SUMMARY=" + samplename+".summary.txt" 
	command = command + " NUM_BARCODES=" + str(ncells*3)
	print(command)
	returnCode = subprocess.call(command,shell=True)
	if(returnCode != 0):
		raise Exception("check incomplete mixing in bead synthesis failed");
	return returnCode

def TagExons(samplename,dropseq_jar,gtf):
	command = 'java -jar '+dropseq_jar+' TagReadWithGeneExon '
	command = command + " I=" + samplename+".bam" + " O=" + samplename+".exon.bam"
	command = command + " ANNOTATIONS_FILE=" + gtf + " TAG=GE "
	print(command)
	returnCode = subprocess.call(command,shell=True)
	if(returnCode != 0):
		raise Exception("tagging file with exon failed");
	return returnCode

def GetExpression(samplename,dropseq_jar,ncells):
	command = 'java -jar '+dropseq_jar+' DigitalExpression '
	command = command + " INPUT=" + samplename+".bam" + " OUTPUT=" + samplename+".dge.txt.gz"
	command = command + " SUMMARY=" + samplename+".dge.summary.txt" 
	command = command + " NUM_CORE_BARCODES=" + str(ncells)
	print(command)
	returnCode = subprocess.call(command,shell=True)
	if(returnCode != 0):
		raise Exception("digital exrpession failed");
	return returnCode

def GetNCells(samplename,dropseq_jar):
	command = 'java -jar '+dropseq_jar+' BAMTagHistogram '
	command = command + " INPUT=" + samplename+".bam" + " OUTPUT=" + samplename+".cellcounts.txt.gz"
	command = command + " TAG=XC" 
	print(command)
	returnCode = subprocess.call(command,shell=True)
	if(returnCode != 0):
		raise Exception("estimate number of cells failed");
	return returnCode
