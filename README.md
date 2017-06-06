# UMIpipe.py

* This python script converts fastq file of dropseq sequencing data to an expression matrix where each column correspond to a cell, and each row correspond to a gene. 
* The dependencies of this pipeline include
  * picard-tools
  * dropseq-tools (included in the depository) 
  * STAR
  * samtools
  * references files: fasta, a picard .dict file in the same path as the fasta file, STAR index, a gtf file. The paths to the file is pre-specified and the user only need to specify the species using the --ref option. For now mm10 and hg38 are supported. 
* The script take a number of arguments at the beginning (see using -h), most of which have default values adapted to running on the yosef2 queue. An example of command is in runUMIpipe.sh. It has 5 basic parts:
  * Convert fastq to sam: requires 3 arguments 
    * --fq1: fastq read 1
    * --fq2: fastq read 2
    * --samplename: output name 
  * Tag barcode: this lets tag in bam files. The cell barcode tag is attached to an optional field in the sam file with the non-barcode read XC. The molecular barcode is attached to field XM. It assumes by default that read 1 contains the barcode sequence, and that cell barcode is base 1-12, and the molecular barcode is 13-20. At the end of the tagging, the first read is discarded. The default 
  
  
