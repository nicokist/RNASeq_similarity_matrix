#!/usr/bin/env python
import os
import subprocess

def call_and_check(command):
    print('Running: "' + command + '"')
    ret = subprocess.call(command, shell=True)
    if ret == 0:
        print("Success")
    if ret != 0:
        raise ValueError("non-zero return code")


call_and_check(
    "/app/samtools faidx /app/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
)
call_and_check(
    "/app/gatk-4.1.0.0/gatk CreateSequenceDictionary -R /app/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
)
call_and_check(
    "/app/gatk-4.1.0.0/gatk IndexFeatureFile -F /app/00-common_all.vcf.gz "
)


bam_files = [i for i in os.listdir('bam_files') if i.endswith('.bam')]
roots=[filename[0:-4] for filename in bam_files]
try:
    os.mkdir('intermediate_files')
except OSError:
    pass

for root in roots:
    call_and_check(
        "/app/samtools index bam_files/{0}.bam".format(root)
    )


    call_and_check("/app/samtools view -b bam_files/{0}.bam 1 > intermediate_files/{0}.chr1.bam".format(root))


    call_and_check("java -jar /app/picard.jar AddOrReplaceReadGroups I=intermediate_files/{0}.chr1.bam O=intermediate_files/{0}.chr1.RGs.bam CREATE_INDEX=true RGID='XX' RGLB='XX' RGPL='XX' RGPU='XX' RGSM='${0}'".format(root))

    call_and_check("java -jar /app/picard.jar MarkDuplicates I=intermediate_files/{0}.chr1.RGs.bam O=intermediate_files/{0}.chr1.RGs.no_duplicates.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT REMOVE_SEQUENCING_DUPLICATES=true M=intermediate_files/{0}.metrics".format(root))

    call_and_check("/app/gatk-4.1.0.0/gatk SplitNCigarReads -R /app/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -L 1 -I intermediate_files/{0}.chr1.RGs.no_duplicates.bam -O intermediate_files/{0}.chr1.RGs.no_duplicates.split_cigar.bam".format(root))

    call_and_check("/app/gatk-4.1.0.0/gatk BaseRecalibrator -R /app/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -L 1 --known-sites /app/00-common_all.vcf.gz -I intermediate_files/{0}.chr1.RGs.no_duplicates.split_cigar.bam -O intermediate_files/{0}.recalibration.table".format(root))

    call_and_check("/app/gatk-4.1.0.0/gatk ApplyBQSR -R /app/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -L 1 --bqsr-recal-file intermediate_files/{0}.recalibration.table -I intermediate_files/{0}.chr1.RGs.no_duplicates.split_cigar.bam -O intermediate_files/{0}.chr1.RG.NoDup.Split.BQSR.bam".format(root))

    call_and_check("/app/gatk-4.1.0.0/gatk HaplotypeCaller -R /app/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -I intermediate_files/{0}.chr1.RG.NoDup.Split.BQSR.bam -L 1 --dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20.0 --dbsnp /app/00-common_all.vcf.gz -O intermediate_files/{0}.chr1.vcf.gz -A QualByDepth -A Coverage -A ClippingRankSumTest -A DepthPerSampleHC".format(root))
      
    call_and_check('/app/gatk-4.1.0.0/gatk VariantFiltration -V intermediate_files/{0}.chr1.vcf.gz -window 35 -cluster 3 --filter-name "FS" --filter-expression "FS > 30.0" --filter-name "QD" --filter-expression "QD < 2.0" -O intermediate_files/{0}.chr1.filtered.vcf.gz'.format(root))

call_and_check("bcftools merge intermediate_files/*.chr1.filtered.vcf.gz > intermediate_files/merged.chr1.vcf.gz")

call_and_check("R CMD BATCH /app/vcf_to_similarity_matrix.R")
