FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="tecap"
LABEL org.opencontainers.image.description="3' terminal exon capture diagnostics for long-read scRNA-seq"
LABEL org.opencontainers.image.source="https://github.com/FullLengthFanatic/tecap"
LABEL org.opencontainers.image.licenses="MIT"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

USER mambauser
ENV MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install -n base -y -c conda-forge -c bioconda \
        python=3.12 \
        pysam">=0.22" \
        intervaltree">=3.1" \
        matplotlib-base">=3.5" \
        numpy">=1.22" \
        pip \
    && micromamba clean --all --yes

WORKDIR /opt/tecap
COPY --chown=mambauser:mambauser . /opt/tecap

RUN pip install --no-deps --no-build-isolation .

ENTRYPOINT ["tecap"]
CMD ["--help"]
