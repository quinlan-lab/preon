FROM continuumio/miniconda3:4.8.2

RUN apt-get update && apt-get install -y procps
RUN conda create -n preon -c conda-forge -c bioconda openjdk==8.0.192 snpeff==4.3.1t snpsift==4.3.1t bcftools==1.10.2 bedtools==2.29.2 fastp==0.20.1 alignstats==0.9.1 && conda clean -a
RUN echo "source activate preon" > ~/.bashrc
ENV PATH /opt/conda/envs/preon/bin:$PATH

# Setup default ARG variables
ARG GENOME=GRCh38
ARG SNPEFF_CACHE_VERSION=86

# Download Genome
RUN snpEff download -v ${GENOME}.${SNPEFF_CACHE_VERSION}
RUN wget https://igor.sbgenomics.com/downloads/cli/sbg-uploader.tgz && tar xzf sbg-uploader.tgz && rm sbg-uploader.tgz && cp sbg-uploader/bin/* /opt/conda/envs/preon/bin && cp sbg-uploader/lib/* /opt/conda/envs/preon/lib
