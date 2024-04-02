#!/usr/bin/env bash                                                                                                                                                                 
i=1
while read line; do
((i++))
varname="var$i"
printf -v $varname "$line"
done < sample_list6.txt
for j in `seq 2 $i`; do
curr_var=var$j
eval curr_var=\$$curr_var
if [ "$curr_var" != "" ]; then
id_name=`echo $curr_var | awk 'END {print $1}'`

ml python/3.6.9

#FastQC
time /users/gshankar/tools/FastQC/fastqc  /dcs04/mathias/data/CAAPA2_Transcriptomics/2102UNHS-0288/${id_name}/${id_name}_1.fastq.gz -o /dcs04/mathias/data/CAAPA2_Transcriptomics/fastQC_reports/
time /users/gshankar/tools/FastQC/fastqc  /dcs04/mathias/data/CAAPA2_Transcriptomics/2102UNHS-0288/${id_name}/${id_name}_2.fastq.gz -o /dcs04/mathias/data/CAAPA2_Transcriptomics/fastQC_reports/


#adapter trimming
time /users/gshankar/tools/bbmap/bbduk.sh -Xmx1g \
in1="/dcs04/mathias/data/CAAPA2_Transcriptomics/2102UNHS-0288/${id_name}/${id_name}_1.fastq.gz" \
in2="/dcs04/mathias/data/CAAPA2_Transcriptomics/2102UNHS-0288/${id_name}/${id_name}_2.fastq.gz" \
out1="/dcs04/mathias/data/CAAPA2_Transcriptomics/trimmed_fastq/${id_name}_1_trimmed.fastq.gz" \
out2="/dcs04/mathias/data/CAAPA2_Transcriptomics/trimmed_fastq/${id_name}_2_trimmed.fastq.gz" \
overwrite=t literal=AGATCGGAAGAGCACACGTCT rcomp=t ktrim=r k=20 mink=11 \
minlen=20 qtrim=rl ftm=5 trimq=10 hdist=1 tpe tbo

#FastQC og trimmed fastq files
time /users/gshankar/tools/FastQC/fastqc /dcs04/mathias/data/CAAPA2_Transcriptomics/trimmed_fastq/${id_name}_1_trimmed.fastq.gz -o /dcs04/mathias/data/CAAPA2_Transcriptomics/fastQC_reports/
time /users/gshankar/tools/FastQC/fastqc /dcs04/mathias/data/CAAPA2_Transcriptomics/trimmed_fastq/${id_name}_2_trimmed.fastq.gz -o /dcs04/mathias/data/CAAPA2_Transcriptomics/fastQC_reports/

#Hisat2 alignment
/users/gshankar/tools/hisat2-2.2.1/hisat2 --time --dta -p 8 \
-x /dcl01/mathias1/data/HISAT2_hg38/genome \
-1 /dcs04/mathias/data/CAAPA2_Transcriptomics/trimmed_fastq/${id_name}_1_trimmed.fastq.gz \
-2 /dcs04/mathias/data/CAAPA2_Transcriptomics/trimmed_fastq/${id_name}_2_trimmed.fastq.gz \
-S /dcs04/mathias/data/CAAPA2_Transcriptomics/alignment/${id_name}_hisat2_hits_adaptertrimmed_redo.sam

echo $?
sleep 1 # pause                                                                                                                                                             
fi
done
