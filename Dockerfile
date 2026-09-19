FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/usr/local/bin:${PATH}"

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    ncbi-blast+ \
    hmmer \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /usr/local/amrfinder && cd /usr/local/amrfinder && \
    URL=$(curl -s https://api.github.com/repos/ncbi/amr/releases/latest | grep "browser_download_url.*amrfinder_binaries" | cut -d '"' -f 4) && \
    curl -sL "$URL" -o amrfinder_bin.tar.gz && \
    tar -xzf amrfinder_bin.tar.gz && \
    rm -f amrfinder_bin.tar.gz && \
    ln -sf /usr/local/amrfinder/amrfinder* /usr/local/bin/ && \
    ln -sf /usr/local/amrfinder/dna_mutation* /usr/local/bin/ 2>/dev/null || true && \
    ln -sf /usr/local/amrfinder/fasta_check* /usr/local/bin/ 2>/dev/null || true

RUN amrfinder -u

WORKDIR /app

COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

COPY . .

# Render defaults to port 10000 (or uses $PORT environment variable)
ENV PORT=10000
EXPOSE 10000

CMD ["sh", "-c", "uvicorn app:app --host 0.0.0.0 --port ${PORT:-10000}"]
