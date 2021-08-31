#!/bin/bash
module purge
module load gcc
module load intel
module load jdk

java -jar tools_newcovid19/picard-tools-2.4.1/picard.jar BedToIntervalList \
      I=ref_dat/ref_iontorrent/S04380110_Regions.bed \
      O=ref_dat/ref_iontorrent/SureSelect_Human_All_Exon_V5.interval_list \
      SD=ref_dat/ref_gatk_newcovid19/Homo_sapiens_assembly38.dict

# 50
tools_newcovid19/gatk-4.2.0.0/gatk SplitIntervals \
  -R ref_dat/ref_gatk_newcovid19/Homo_sapiens_assembly38.fasta \
  -L ref_dat/ref_iontorrent/SureSelect_Human_All_Exon_V5.interval_list \
  --scatter-count 50 \
  -O ref_dat/SureSelect_Human_All_Exon_V5_50intervals


for i in `seq -f %04g 0 49`; do
    echo ref_dat/SureSelect_Human_All_Exon_V5_50intervals/$i-scattered.interval_list >> \
      ref_dat/SureSelect_Human_All_Exon_V5_50intervals_index.txt;
done

# 400
tools_newcovid19/gatk-4.2.0.0/gatk SplitIntervals \
  -R ref_dat/ref_gatk_newcovid19/Homo_sapiens_assembly38.fasta \
  -L ref_dat/ref_iontorrent/SureSelect_Human_All_Exon_V5.interval_list \
  --scatter-count 400 \
  -O ref_dat/SureSelect_Human_All_Exon_V5_400intervals


for i in `seq -f %04g 0 399`; do
    echo ref_dat/SureSelect_Human_All_Exon_V5_400intervals/$i-scattered.interval_list >> \
      ref_dat/SureSelect_Human_All_Exon_V5_400intervals_index.txt;
done
