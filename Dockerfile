FROM debian:buster

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       2to3 \
       bedtools \
       ca-certificates \
       curl \
       git \
       python3 \
       python3-matplotlib \
       python3-numpy \
       python3-pysam \
       python3-scipy \
       samtools \
    && rm -rf /var/lib/apt/lists/* \
    && git clone https://github.com/Xinglab/rmats2sashimiplot.git /rmats2sashimiplot \
    && cd /rmats2sashimiplot \
    # && git checkout {commit} \
    && ./2to3.sh

# Set defaults for running the image.
# The ENTRYPOINT AND CMD are empty to be compatible with
# CWL and WDL implementations that cannot override those values
WORKDIR /rmats2sashimiplot
ENTRYPOINT []
CMD []
