FROM ubuntu:16.04
RUN apt-get update
RUN apt-get install -y python
COPY resources/* /app/
WORKDIR /data/
CMD ["python", "/app/RNASeq_sample_confusion.py"]
