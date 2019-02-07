# Generating sample relatedness matrices from RNASeq data
RNASeq data is often used to detetermine transcript abundances for each gene after which the original data is discarded. However, the SNP data obtained by RNASeq can be used to generate relatedness.

This can be done as follows:
- Align the RNASeq data to the genome (gcHR37).
- Generate indices for the alignment bam files and merge them.
`find bam_files/*.bam | xargs -P 30 -n 1 /samtools index`
`/samtools merge --threads 30 -r -R 1 bam_files.merged_chr1.bam bam_files/*.bam`
- Samtools merge does not keep readgroup information, which links the reads to the sample. We regenerate this using the bam filenames.
```
/samtools view -H bam_files.merged_chr1.bam | grep -v '^@RG' > bam_files.merged_chr1.new_header
find bam_files/*.bam | /generate_RGs.py >> bam_files.merged_chr1.new_header"
/samtools reheader bam_files.merged_chr1.new_header bam_files.merged_chr1.bam > bam_files.merged_chr1.header_withRG.bam
```
 - Use gatk MarkDuplicates to keep freebayes from crashing and burning.
 `java -jar /gatk-package-4.0.11.0-local.jar MarkDuplicates --INPUT bam_files.merged_chr1.header_withRG.bam --OUTPUT bam_files.merged_chr1.header_withRG.MarkDuplicates.bam --CREATE_INDEX -M MarkDuplicates.metrics --VALIDATION_STRINGENCY LENIENT`
 
- Call SNPs using freebayes
```
# ~/tools/freebayes/scripts/fasta_generate_regions.py ~/temp/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa 100000 |  grep '^1\:'  > ~/temp/chr1_regions
# presupplied as chr1_regions in resources directory.

ulimit -n 160000; 
        cd /freebayes/scripts; 
        ./freebayes-parallel /chr1_regions 30 --use-best-n-alleles 4 -f /Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -b ~/bam_files.merged_chr1.header_withRG.MarkDuplicates.bam > ~/bam_files.merged_chr1.header_withRG.MarkDuplicates.freebayes_best_4_alleles.vcf;
```

#Limitations
By default this analysis is limited to chromosome 1 as this gives abundant signal and significantly reduces the compututation required.

If the bam file entrypoint is used we require a different BAM file be provided for each sample, as the readgroups will be overwritten using the filename after merging the files.

