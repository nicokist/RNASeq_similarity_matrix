FROM ubuntu:16.04

RUN apt-get update
RUN apt-get install -y tzdata
RUN apt-get install -y python parallel default-jre git build-essential zlib1g-dev libbz2-dev libhts-dev liblzma-dev wget unzip lzma git libcurl4-openssl-dev

RUN mkdir /.parallel
RUN touch /.parallel/will-cite

WORKDIR /app
RUN wget -q ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
RUN gzip -d /app/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

RUN wget -q https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip
RUN unzip gatk-4.1.0.0.zip

RUN wget -q https://github.com/broadinstitute/picard/releases/download/2.20.4/picard.jar

RUN wget -q ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz

#R
RUN apt-get -y install software-properties-common
RUN add-apt-repository "deb http://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/"
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 
RUN apt-get update
RUN apt-get install -y r-base r-base-dev bcftools

COPY app/* /app/
RUN chmod -R 777 /app

WORKDIR /usr/lib/R
RUN tar -xjf /app/usr_lib_R_site-library.tar.bz2


WORKDIR /data/
ENV PATH="/app/:${PATH}"
CMD ["python", "/app/RNASeq_sample_confusion.py"]
