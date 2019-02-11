FROM ubuntu:16.04
RUN apt-get update
RUN apt-get install -y python parallel default-jre
COPY app/* /app/
WORKDIR /app

##Samtools is included as a binary in the repository, gatk and freebayes were too large for github so are obtained properly.

RUN git clone --recursive git://github.com/ekg/freebayes.git
WORKDIR /app/freebayes
RUN make

WORKDIR /app
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip
RUN unzip gatk-4.1.0.0.zip

RUN wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz


WORKDIR /data/
CMD ["python", "/app/RNASeq_sample_confusion.py"]
