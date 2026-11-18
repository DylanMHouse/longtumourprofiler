FROM docker.io/library/rust:1-bookworm AS build

ARG MAFSMITH_REV=e009e984560671c941f709f7c16dbe6505ee1f29
ARG FASTVEP_REV=4ca2a6e36c6785ff4002c3f75a62b4b57a8d6a04

RUN cargo install --root /opt/tools \
        --git https://github.com/nf-osi/mafsmith --rev ${MAFSMITH_REV} mafsmith \
 && cargo install --locked --root /opt/tools \
        --git https://github.com/Huang-lab/fastVEP --rev ${FASTVEP_REV} fastvep-cli \
 && /opt/tools/bin/mafsmith --version \
 && /opt/tools/bin/fastvep --version

FROM docker.io/library/debian:bookworm-slim

RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates procps \
 && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/tools/bin/mafsmith /opt/tools/bin/fastvep /usr/local/bin/

LABEL org.opencontainers.image.source="https://github.com/nf-osi/mafsmith" \
      org.opencontainers.image.description="mafsmith + fastVEP for nf-core/longtumourprofiler"
