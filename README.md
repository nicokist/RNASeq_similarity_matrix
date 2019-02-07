# Generate a sample relatedness matrices from RNASeq data
Sample confusion is a common laboratory problem. In RNASeq this is frequently tested for by checking whether sex-specific genes (y chromosomal, Xist) are congruent with the sex listed for that sample in the metadata. However, this method cannot be used to detect sample confusion between patients of the same sex. Here we leverage RNASeq reads to call genomic SNPs, and use that to generate a relatedness matrix between all samples. Samples from the same patient should be highly similar, samples from different patients should not be similar (assuming patients are unrelated).



RNASeq data is often used to detetermine transcript abundances for each gene after which the original data is discarded. However, the SNP data obtained by RNASeq can be used to generate relatedness.

This can be done as follows:
- Align the RNASeq data to the genome (gcHR37). We recommend the STAR aligner (https://github.com/alexdobin/STAR). The rest of the pipeline assumes you've generated one bam file per sample.
- Generate indices for the alignment bam files and merge them.

```
find bam_files/*.bam | xargs -P 30 -n 1 /samtools index
/samtools merge --threads 30 -r -R 1 bam_files.merged_chr1.bam bam_files/*.bam
```

- Samtools merge does not keep readgroup information, which links the reads to the sample. We regenerate this using the bam filenames.

```
/samtools view -H bam_files.merged_chr1.bam | grep -v '^@RG' > bam_files.merged_chr1.new_header
find bam_files/*.bam | /generate_RGs.py >> bam_files.merged_chr1.new_header"
/samtools reheader bam_files.merged_chr1.new_header bam_files.merged_chr1.bam > bam_files.merged_chr1.header_withRG.bam
```
 - Use gatk MarkDuplicates to keep freebayes from crashing due to deep piles of highly abundant transcripts.

```
java -jar /gatk-package-4.0.11.0-local.jar MarkDuplicates --INPUT bam_files.merged_chr1.header_withRG.bam --OUTPUT bam_files.merged_chr1.header_withRG.MarkDuplicates.bam --CREATE_INDEX -M MarkDuplicates.metrics --VALIDATION_STRINGENCY LENIENT
```
 
- Call SNPs using freebayes

```
# ~/tools/freebayes/scripts/fasta_generate_regions.py ~/temp/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa 100000 |  grep '^1\:'  > ~/temp/chr1_regions
# presupplied as chr1_regions in resources directory.

ulimit -n 160000; 
        cd /freebayes/scripts; 
        ./freebayes-parallel /chr1_regions 30 --use-best-n-alleles 4 -f /Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -b ~/bam_files.merged_chr1.header_withRG.MarkDuplicates.bam > ~/bam_files.merged_chr1.header_withRG.MarkDuplicates.freebayes_best_4_alleles.vcf;
```

- Filter the RNASeq SNP calls by quality.

```
/vcffilter -f 'QUAL > 20' vcf_file > vcf_file.QUAL_GT_20.vcf
```

- Filter the high quality RNASeq SNPs, keeping only common genomic SNPs.

```
grep '^#' vcf_file.QUAL_GT_20.vcf > vcf_file.QUAL_GT_20.common_snps_only.vcf
/bedtools intersect -a vcf_file.QUAL_GT_20.vcf -b /00-common_all.vcf.gz -wa >> vcf_file.QUAL_GT_20.common_snps_only.vcf")
```

- Finally we can generate the similarity matrix.

```
R CMD BATCH vcf_to_similarity_matrix.R
```



# Limitations
By default this analysis is limited to chromosome 1 as this gives abundant signal and significantly reduces the compututation required.

If the bam file entrypoint is used we require a different BAM file be provided for each sample, as the readgroups will be overwritten using the filename after merging the files.

In spite of our stringent filtering we do not fully recover the genomic SNPs, resulting in some leftover difference between different samples from the same patients, however the resulting identity-by-sequence should still be around 98%
