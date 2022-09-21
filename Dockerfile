FROM condaforge/mambaforge
MAINTAINER Xiaofang Jiang
WORKDIR /workdir
COPY environment.yml .
RUN mamba env create -f environment.yml && mamba clean --all --yes && rm environment.yml
ENV PATH /opt/conda/envs/wwbe/bin:$PATH
RUN /bin/bash -c "source activate wwbe"
COPY data/NC_045512.2/snpEffectPredictor.bin /opt/conda/envs/wwbe/share/snpeff-5.0-1/data/NC_045512.2/snpEffectPredictor.bin
RUN chmod 755 /opt/conda/envs/wwbe/share/snpeff-5.0-1/data/NC_045512.2/snpEffectPredictor.bin
