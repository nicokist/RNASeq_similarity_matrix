#!/usr/bin/env python
##Insert code to move bam files to bam_files.
import os
import subprocess


def call_and_check(command):
    print('Running: "' + command + '"')
    ret = subprocess.call(command, shell=True)
    if ret == 0:
        print("Success")
    if ret != 0:
        raise ValueError("non-zero return code")


call_and_check("find bam_files/*.bam | xargs -P `nproc` -n 1 samtools index")
call_and_check(
    "samtools merge --threads `nproc` -r -R 1 bam_files.merged_chr1.bam bam_files/*.bam"
)
call_and_check(
    "samtools view -H bam_files.merged_chr1.bam | grep -v '^@RG' > bam_files.merged_chr1.new_header"
)
call_and_check(
    "find bam_files/*.bam | generate_RGs.py >> bam_files.merged_chr1.new_header"
)
call_and_check(
    "samtools reheader bam_files.merged_chr1.new_header bam_files.merged_chr1.bam > bam_files.merged_chr1.header_withRG.bam"
)
call_and_check(
    "java -jar /app/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar MarkDuplicates --INPUT bam_files.merged_chr1.header_withRG.bam --OUTPUT bam_files.merged_chr1.header_withRG.MarkDuplicates.bam --CREATE_INDEX -M MarkDuplicates.metrics --VALIDATION_STRINGENCY LENIENT"
)

# ~/tools/freebayes/scripts/fasta_generate_regions.py ~/temp/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa 100000 |  grep '^1\:'  > ~/temp/
# presupplied as chr1_regions in resources directory.

call_and_check(
    "samtools faidx /app/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
)



call_and_check(
    """ulimit -n 160000;
    cd /app/freebayes/scripts;
    ./freebayes-parallel /app/chr1_regions `nproc` --use-best-n-alleles 4 -f /app/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -b /data/bam_files.merged_chr1.header_withRG.MarkDuplicates.bam > /data/bam_files.merged_chr1.header_withRG.MarkDuplicates.freebayes_best_4_alleles.vcf;
    """
)


call_and_check(
    "vcffilter -f 'QUAL > 20' bam_files.merged_chr1.header_withRG.MarkDuplicates.freebayes_best_4_alleles.vcf > bam_files.merged_chr1.header_withRG.MarkDuplicates.freebayes_best_4_alleles.QUAL_GT_20.vcf"
)
call_and_check(
    "grep '^#' bam_files.merged_chr1.header_withRG.MarkDuplicates.freebayes_best_4_alleles.QUAL_GT_20.vcf > bam_files.merged_chr1.header_withRG.MarkDuplicates.freebayes_best_4_alleles.QUAL_GT_20.common_snps_only.vcf"
)
call_and_check(
    "bedtools intersect -a bam_files.merged_chr1.header_withRG.MarkDuplicates.freebayes_best_4_alleles.QUAL_GT_20.vcf -b /app/00-common_all.vcf.gz -wa >> bam_files.merged_chr1.header_withRG.MarkDuplicates.freebayes_best_4_alleles.QUAL_GT_20.common_snps_only.vcf"
)
call_and_check("R CMD BATCH /app/vcf_to_similarity_matrix.R")
