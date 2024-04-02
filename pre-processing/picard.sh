#!/usr/bin/env bash                                                                                                                                                                 
i=1
while read line; do
((i++))
varname="var$i"
printf -v $varname "$line"
done < sample_list/sample_names28
for j in `seq 2 $i`; do
curr_var=var$j
eval curr_var=\$$curr_var
if [ "$curr_var" != "" ]; then
id_name=`echo $curr_var | awk 'END {print $1}'`

ml java

#time samtools flagstat /dcs04/mathias/data/CAAPA2_Transcriptomics/alignment/BAM/${id_name}_hisat2_hits_adaptertrimmed.bam 
java -jar /users/gshankar/tools/picard.jar CollectRnaSeqMetrics I=/dcs04/mathias/data/CAAPA2_Transcriptomics/alignment/BAM/${id_name}_hisat2_hits_adaptertrimmed.bam O=/dcs04/mathias/data/CAAPA2_Transcriptomics/picard_rRNA/${id_name}.RNAmetrics.txt REF_FLAT=/dcs04/mathias/data/CAAPA2_Transcriptomics/refFlat.txt STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=/dcs04/mathias/data/CAAPA2_Transcriptomics/hg38_ribosomal_interval.list

echo $?
sleep 1 # pause                                                                                                                                                             
fi
done
