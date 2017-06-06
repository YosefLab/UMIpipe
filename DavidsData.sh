while file= read -r var; do echo "$var"; done < filenames.txt

while IFS= read -r var
do
   samplename=$(echo "$var" |cut -f 3 -d '/'); sed "s/SAMPLE/$samplename/" runUMIpipe.sh > runUMIpipe_$samplename.sh
done < filenames.txt 

grep '_R1_' filenames.txt > read1names.txt
grep '_R2_' filenames.txt > read2names.txt

while IFS= read -r var
do
   samplename=$(echo "$var" |cut -f 3 -d '/')
   sed -i "s|FQ1|$var|" runUMIpipe_$samplename.sh 
   echo $var
   echo $samplename
done < read1names.txt


while IFS= read -r var
do
   samplename=$(echo "$var" |cut -f 3 -d '/')
   sed -i "s|FQ2|$var|" runUMIpipe_$samplename.sh 
   echo $var
   echo $samplename
done < read2names.txt


for f in *.sec.merged.bam
do
samtools index $f
#tag could also be XM
umi_tools dedup -I $f -S $g.dedup.bam -L $g.dedup.log --extract-umi-method=tag --umi-tag='MC' \
--output-stats=$g.dedup
done

