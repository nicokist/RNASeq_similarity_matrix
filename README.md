# Generate a sample relatedness matrices from RNASeq data
Sample confusion is a common laboratory problem. In RNASeq this is frequently tested for by checking whether sex-specific genes (y chromosomal, Xist) are congruent with the sex listed for that sample in the metadata. However, this method cannot be used to detect sample confusion between patients of the same sex. Here we leverage RNASeq reads to call genomic SNPs, and use that to generate a relatedness matrix between all samples. 

Samples from the same patient should be highly similar, samples from different patients should not be similar (assuming patients are unrelated).

RNASeq data is often used to detetermine transcript abundances for each gene after which the original data is discarded. However, the SNP data obtained by RNASeq can be used to generate relatedness.

Doing this requires a number of tools and commands, we have wrapped these in a docker image. 

## Instructions
- Align the RNASeq data to the genome (gcHR37). We recommend the STAR aligner (https://github.com/alexdobin/STAR). The rest of the pipeline assumes you've generated one sorted bam file per sample, place these in an otherwise empty directory called `bam_files`. 
- Download and build the docker image.
```
git clone https://github.com/nicokist/RNASeq_sample_relatedness_matrix
docker build RNASeq_sample_relatedness_matrix -t rnaseq_relatedness
```
- Run the following docker invocation on a machine with sufficient memory and cores. This may take a few days if you have many samples. It is recommended to start with just three bam files to see if the pipeline completes.
```
sudo docker run --user `id -u`:`id -g` -it -v `pwd`:/data rnaseq_relatedness
```
- The sequence similarity matrix and a visualization thereof will be left in the current directory if everything finishes successfully.


#If you don't have bam files
If you don't have GRCh38 bam files it is recommended to use STAR to align your fastq files to the reference genome, creating a single bam file for each sample.

First generate a genome index:
```
wget 'ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz'
gzip -d Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz


wget 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz'
gzip -d gencode.v29.annotation.gtf.gz


docker build https://raw.githubusercontent.com/alexdobin/STAR/master/extras/docker/Dockerfile -t star
mkdir Homo_sapiens.GRCh38.dna_sm.primary_assembly.star_genome
docker run -v `pwd`:/data star STAR --runThreadN `nproc` --runMode genomeGenerate --genomeDir /data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.star_genome --genomeFastaFiles /data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile /data/gencode.v29.primary_assembly.annotation.gtf
```

Then generate a bam file by running the following command for each sample:
```

```

# Technical Details
We use Docker to deliver a linux image with everything needed for the pipeline pre-installed. This includes a script which executes all steps in turn. If for some reason you wish to use one of the supplied programs manually you can use the following invocation, with the desired command in quotation marks after `bash -c` at the end:

```
sudo docker run  -it -v `pwd`:/data rnaseq_relatedness bash -c "find bam_files/*.bam | xargs -P 30 -n 1 /samtools index"
```

# Limitations
By default this analysis is limited to chromosome 1 as this gives abundant signal and significantly reduces the compututation required.

If the bam file entrypoint is used we require a different BAM file be provided for each sample, as the readgroups will be overwritten using the filename after merging the files.

In spite of our stringent filtering we do not fully recover the genomic SNPs, resulting in some leftover difference between different samples from the same patients, however the resulting identity-by-sequence should still be around 98%
