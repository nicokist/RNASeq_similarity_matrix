FROM ubuntu:16.04

RUN apt-get update
RUN apt-get install -y tzdata
RUN apt-get install -y python parallel default-jre git build-essential zlib1g-dev libbz2-dev libhts-dev liblzma-dev wget unzip lzma git r-base libcurl4-openssl-dev

WORKDIR /app
RUN wget -q ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
RUN gzip -d /app/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip
RUN unzip gatk-4.1.0.0.zip


#R
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/"
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN apt-get update
RUN apt-get install -y r-base r-base-dev


##Samtools is included as a binary in the repository, gatk and freebayes were too large for github so are obtained properly.
RUN git clone --recursive git://github.com/ekg/freebayes.git
WORKDIR /app/freebayes
RUN make
RUN cp /app/freebayes/bin/freebayes /usr/local/bin

WORKDIR /app/freebayes/vcflib
RUN make
RUN cp bin/* /usr/local/bin/

COPY app/* /app/
WORKDIR /data/
CMD ["python", "/app/RNASeq_sample_confusion.py"]