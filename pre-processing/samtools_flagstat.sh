#!/usr/bin/env bash                                                                                                                                                                 
i=1
while read line; do
((i++))
varname="var$i"
printf -v $varname "$line"
done < sample_list/sample_list_flagstat
for j in `seq 2 $i`; do
curr_var=var$j
eval curr_var=\$$curr_var
if [ "$curr_var" != "" ]; then
id_name=`echo $curr_var | awk 'END {print $1}'`

ml python/3.7.3
ml samtools
export PATH=$PATH:/users/gshankar/tools
export PATH=$PATH:/users/gshankar/tools/subread-2.0.2-Linux-x86_64/bin

time samtools flagstat /dcs04/mathias/data/CAAPA2_Transcriptomics/alignment/BAM/${id_name}_hisat2_hits_adaptertrimmed.bam 

echo $?
sleep 1 # pause                                                                                                                                                             
fi
done
