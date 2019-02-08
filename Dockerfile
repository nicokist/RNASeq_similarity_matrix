FROM ubuntu:16.04
RUN apt-get update
RUN apt-get install -y python
COPY . /
CMD ["python", "RNASeq_sequence_similarity_matrix.py"]
