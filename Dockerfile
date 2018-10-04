# Dockerfile for the sequence refinement pipeline implemented in Nextflow

FROM nfcore/base

MAINTAINER Julia F. Soellner <julia.soellner@boehringer-ingelheim.com>

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/ao-tool_v1/bin:$PATH

CMD ["/bin/bash"]

