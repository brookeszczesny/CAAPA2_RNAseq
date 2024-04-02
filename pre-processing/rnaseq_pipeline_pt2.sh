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

#sort SAM file and convert to BAM
#time samtools sort -o /dcs04/mathias/data/CAAPA2_Transcriptomics/alignment/BAM/${id_name}_hisat2_hits_adaptertrimmed.bam /dcs04/mathias/data/CAAPA2_Transcriptomics/alignment/${id_name}_hisat2_hits_adaptertrimmed.sam

#generate read counts
#time /users/gshankar/tools/coco/bin/coco correct_count /dcs04/mathias/data/CAAPA2_Transcriptomics/gtf/Homo_sapiens.GRCh38.104.correct_annotation.gtf /dcs04/mathias/data/CAAPA2_Transcriptomics/alignment/BAM/${id_name}_hisat2_hits_adaptertrimmed.bam /dcs04/mathias/data/CAAPA2_Transcriptomics/read_counts/${id_name}_CoCo.txt --paired

#echo "gene_id,${id_name}" > /dcs04/mathias/data/CAAPA2_Transcriptomics/read_counts/${id_name}_raw_CoCo_counts.txt
#awk 'NR > 1 {print $1 "-" $2 "," $3}' /dcs04/mathias/data/CAAPA2_Transcriptomics/read_counts/${id_name}_CoCo.txt >> /dcs04/mathias/data/CAAPA2_Transcriptomics/read_counts/${id_name}_raw_CoCo_counts.txt

#echo "gene_id,${id_name}" > /dcs04/mathias/data/CAAPA2_Transcriptomics/read_counts/${id_name}_TPM_CoCo_counts.txt
#awk 'NR > 1 {print $1 "-" $2 "," $5}' /dcs04/mathias/data/CAAPA2_Transcriptomics/read_counts//${id_name}_CoCo.txt >> /dcs04/mathias/data/CAAPA2_Transcriptomics/read_counts/${id_name}_TPM_CoCo_counts.txt

echo "gene_id,${id_name}" > /dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/read_counts/${id_name}_CPM_CoCo_counts.txt
awk 'NR > 1 {print $1 "-" $2 "," $4}' /dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/read_counts//${id_name}_CoCo.txt >> /dcs04/mathias/data/CAAPA2_Transcriptomics_analysis/read_counts/${id_name}_CPM_CoCo_counts.txt

echo $?
sleep 1 # pause                                                                                                                                                             
fi
done
