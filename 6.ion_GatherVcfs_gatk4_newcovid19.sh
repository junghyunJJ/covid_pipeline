#!/bin/bash
#for i in {1..50}; do echo "-I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_${i}.vcf.gz \ "; done

module purge
module load gcc
module load jdk
module load bcftools

tools_newcovid19/gatk-4.2.0.0/gatk GatherVcfs \
    -R ref_dat/ref_gatk_newcovid19/Homo_sapiens_assembly38.fasta \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_1.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_2.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_3.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_4.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_5.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_6.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_7.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_8.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_9.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_10.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_11.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_12.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_13.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_14.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_15.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_16.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_17.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_18.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_19.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_20.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_21.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_22.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_23.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_24.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_25.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_26.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_27.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_28.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_29.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_30.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_31.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_32.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_33.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_34.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_35.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_36.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_37.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_38.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_39.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_40.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_41.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_42.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_43.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_44.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_45.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_46.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_47.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_48.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_49.vcf.gz \
    -I ion_res_gatk4_newcovid19/ion_res_vcf_gatk4_newcovid19/ion_gatk4_newcovid19_50.vcf.gz \
    -O ion_res_gatk4_newcovid19/ion_gatk4_newcovid19.vcf.gz

bcftools index -t ion_res_gatk4_newcovid19/ion_gatk4_newcovid19.vcf.gz
