# Combined test image for sceQTL-flow
# Includes packages from all pipeline Dockerfiles + tensorqtl + development tools

FROM r-base:4.4.0

# Install necessary system dependencies
RUN apt-get update && apt-get install -y \
    software-properties-common \
    && apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    liblapack-dev \
    libopenblas-dev \
    liblapacke-dev \
    procps \
    pkg-config \
    wget \
    make \
    gcc \
    g++ \
    gfortran \
    binutils \
    bzip2 \
    texlive-full \
    pandoc \
    git \
    vim \
    curl \
    ca-certificates \
    htop \
    python3 \
    python3-pip \
    python3-dev \
    zlib1g-dev \
    && apt-get upgrade -y gcc g++ binutils \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install PLINK 1.9
RUN wget -qO plink1.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip && \
    unzip plink1.zip -d /usr/local/bin/ && \
    rm plink1.zip

# Install PLINK 2.0
RUN wget -qO plink2.zip https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241222.zip && \
    unzip plink2.zip -d /usr/local/bin/ && \
    rm plink2.zip

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.23/bcftools-1.23.tar.bz2 && \
    tar -xvjf bcftools-1.23.tar.bz2 && \
    cd bcftools-1.23 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf bcftools-1.23 bcftools-1.23.tar.bz2

# Install Python packages and plink-ng Python bindings for tensorqtl
RUN pip3 install --no-cache-dir --break-system-packages \
    numpy \
    pandas \
    scipy \
    scikit-learn \
    Cython \
    rpy2 \
    Pgenlib \
    pandas-plink \
    setuptools \
    wheel

# Pgenlib Python bindings (for PLINK2 pgen/pvar/psam support)
# Installed via pip above; avoid building plink-ng from source to prevent BLAS/ATLAS link issues.

# Install tensorqtl after dependencies
RUN pip3 install --no-cache-dir --break-system-packages tensorqtl

# Install BiocManager
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'

# Install qvalue for tensorqtl
RUN R -e 'BiocManager::install("qvalue")'

# Install expression analysis packages
RUN R -e 'install.packages(c("dplyr","data.table","Seurat","remotes","presto","glmGamPoi"))'

# Install genotype analysis packages
RUN R -e 'BiocManager::install(c("SeqArray","SNPRelate","GenomicRanges","BSgenome"))'
RUN R -e 'BiocManager::install(c("SNPlocs.Hsapiens.dbSNP155.GRCh38","SNPlocs.Hsapiens.dbSNP155.GRCh37","GenomeInfoDb"))'
RUN R -e 'BiocManager::install(c("TxDb.Hsapiens.UCSC.hg38.knownGene","org.Hs.eg.db","edgeR","EnsDb.Hsapiens.v86"))'

# Install QTL and statistical packages
RUN R -e 'install.packages(c("argparse","MatrixEQTL","coloc","MendelianRandomization","ieugwasr","genetics.binaRies","lobstr","arrow"))'

# Install reporting packages
RUN R -e 'install.packages(c("ggplot2","cowplot","patchwork","ggrepel","tidyverse","gridExtra","ggbeeswarm","rmarkdown","tinytex"), repos="http://cran.rstudio.com/")'

# Install BPCells from GitHub
RUN R -e 'remotes::install_github("bnprks/BPCells/r")'

# Download the chain file for liftOver
RUN wget -P /usr/local/src https://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
RUN gunzip /usr/local/src/hg19ToHg38.over.chain.gz

# Create vscode user for development
RUN useradd -m -s /bin/bash vscode

# Set the working directory
WORKDIR /workspace

# Set the default command to R
CMD ["R"]
