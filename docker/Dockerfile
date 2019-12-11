FROM continuumio/miniconda3:4.7.12

RUN apt-get update && apt-get install -y procps
RUN conda create -n preon -c conda-forge openjdk==8.0.192
RUN echo "source activate preon" > ~/.bashrc
ENV PATH /opt/conda/envs/preon/bin:$PATH

RUN wget https://igor.sbgenomics.com/downloads/cli/sbg-uploader.tgz && tar xzf sbg-uploader.tgz && rm sbg-uploader.tgz && cp sbg-uploader/bin/* /opt/conda/envs/preon/bin && cp sbg-uploader/lib/* /opt/conda/envs/preon/lib
