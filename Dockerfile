FROM debian:trixie

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       ca-certificates \
       curl \
       git \
    && rm -rf /var/lib/apt/lists/* \
    # Install conda to /conda
    && mkdir /conda \
    && cd /conda \
    && curl -L 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh' -O \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /conda/install \
    && /conda/install/bin/conda init \
    && echo '' > /conda/install/.condarc \
    && git clone https://github.com/Xinglab/rmats2sashimiplot.git /rmats2sashimiplot \
    && cd /rmats2sashimiplot \
    # && git checkout {commit} \
    && /conda/install/bin/conda install -c conda-forge -c bioconda --file conda_requirements.txt \
    && /conda/install/bin/python -m pip install .

# Make conda installed programs available on PATH
ENV PATH /conda/install/bin:${PATH}

# Set defaults for running the image.
# The ENTRYPOINT AND CMD are empty to be compatible with
# CWL and WDL implementations that cannot override those values
WORKDIR /rmats2sashimiplot
ENTRYPOINT []
CMD []
