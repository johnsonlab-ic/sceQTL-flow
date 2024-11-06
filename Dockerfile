# Use the official R base image
FROM r-base:latest

# Install necessary system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

# Install BiocManager
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'

RUN R -e 'install.packages(c("dplyr","data.table"))'

# Install the required R packages
RUN R -e 'BiocManager::install(c("seqArray", "SNPRelate", "GenomicRanges", "BSgenome", "SNPlocs.Hsapiens.dbSNP155.GRCh38","SNPlocs.Hsapiens.dbSNP155.GRCh37","GenomeInfoDb"))'

# Set the working directory
WORKDIR /usr/local/src

# Copy any additional scripts or files if needed
# COPY . .

# Set the default command to R
CMD ["R"]