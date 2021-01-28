
# this is the base image
FROM continuumio/miniconda3

# set the working directory that we can copy files into
WORKDIR /usr/src/app

# add the conda-forge channel in case we need to install a package from there
RUN conda config --add channels conda-forge

# copy the yml file into our workdir
COPY environment.yml .
# Install biology environment into base conda environment
RUN /opt/conda/bin/conda env update -n base -f environment.yml

# Test installation
RUN /bin/bash -c "echo -e \"import shapely.geometry\" | python"