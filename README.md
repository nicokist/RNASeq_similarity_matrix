# Generate a sample similarity matrices from RNASeq data
Sample confusion is a common laboratory problem. In RNASeq this is frequently tested for by checking whether sex-specific genes (y chromosomal, Xist) are congruent with the sex listed for that sample in the metadata. However, this method cannot be used to detect sample confusion between patients of the same sex. Here we leverage RNASeq reads to call genomic SNPs, and use that to generate a similarity matrix between all samples. 

Samples from the same patient should be highly similar, samples from different patients should not be similar (assuming patients are unrelated).

RNASeq data is often used to detetermine transcript abundances for each gene after which the original data is discarded. However, the SNP data obtained by RNASeq can be used to generate a similarity matrix.

Doing this requires a number of tools and commands, we have wrapped these in a docker image. 

## Instructions
- Align the RNASeq data to the genome (gcHR37). We recommend the STAR aligner (https://github.com/alexdobin/STAR). The rest of the pipeline assumes you've generated one sorted bam file per sample, place these in an otherwise empty directory called `bam_files`. 
- Download and build the docker image.
```
git clone https://github.com/nicokist/RNASeq_similarity_matrix
docker build RNASeq_similarity_matrix -t rnaseq_similarity_matrix
```
- Run the following docker invocation on a machine with sufficient memory and cores. This may take a few days if you have many samples. It is recommended to start with just three bam files to see if the pipeline completes.
```
docker run --user `id -u`:`id -g` -it -v `pwd`:/data rnaseq_similarity_matrix
```
- The sequence similarity matrix and a visualization thereof will be left in the current directory if everything finishes successfully.

# Citation
Link to preprint to follow.

# Troubleshooting
## If you don't have bam files
If you don't have GRCh38 bam files it is recommended to use STAR to align your fastq files to the reference genome, creating a single bam file for each sample.

First generate a genome index:
```
wget 'ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz'
gzip -d Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz


wget 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz'
gzip -d gencode.v29.primary_assembly.annotation.gtf.gz


docker build https://raw.githubusercontent.com/alexdobin/STAR/master/extras/docker/Dockerfile -t star
mkdir Homo_sapiens.GRCh38.dna_sm.primary_assembly.star_genome
docker run -v `pwd`:/data star STAR --runThreadN `nproc` --runMode genomeGenerate --genomeDir /data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.star_genome --genomeFastaFiles /data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile /data/gencode.v29.primary_assembly.annotation.gtf
```

Then generate a bam file by running the following command, which will align each fastq file to the reference separately (single unpaired reads):
```
find fasta_files/*.gz | xargs -t -I {} -n 1 docker run -v `pwd`:/data star STAR --runThreadN `nproc` --genomeDir /data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.star_genome  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn /data/{} --outFileNamePrefix /data/{}

```

If you have paired reads you will need to write a quick script rather than use xargs as above. See STAR manual for details.


## Error: 'Region "1" specifies an unkown reference name'
If you get the following output you need to convert the IDs of the bam file.

```
Running: "find bam_files/*.bam | xargs -P `nproc` -n 1 /app/samtools index"
Success
Running: "/app/samtools merge --threads `nproc` -r -R 1 bam_files.merged_chr1.bam bam_files/*.bam"
[bam_merge_core2] Region "1" specifies an unknown reference name
Traceback (most recent call last):
  File "/app/RNASeq_sample_confusion.py", line 18, in <module>
    "/app/samtools merge --threads `nproc` -r -R 1 bam_files.merged_chr1.bam bam_files/*.bam"
  File "/app/RNASeq_sample_confusion.py", line 13, in call_and_check
    raise ValueError("non-zero return code")
ValueError: non-zero return code
```

Rename the bam_files directory to original_bam_files and run the following, which will generate a new bam_files directory with appropriate sequence IDs in the bam files.

```
docker run --user `id -u`:`id -g` -it -v `pwd`:/data rnaseq_similarity_matrix remove_chr_from_bam.sh
```

## The program fails unexpectedly
Perhaps you ran out of memory? Run `dmesg` and check if the oom-killer ended the process.

# Technical Details
We use Docker to deliver a linux image with everything needed for the pipeline pre-installed. This includes a script which executes all steps in turn. If for some reason you wish to use one of the supplied programs manually you can use the following invocation, with the desired command in quotation marks after `bash -c` at the end:

```
sudo docker run  -it -v `pwd`:/data rnaseq_similarity_matrix bash -c "find bam_files/*.bam | xargs -P 30 -n 1 /samtools index"
```

## What the docker invocation does.
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

In spite of our stringent filtering we do not fully recover the genomic SNPs, resulting in some leftover difference between different samples from the same patients.
