FROM continuumio/miniconda3

WORKDIR /app

# Install system dependencies for Qt/X11
RUN apt-get update && apt-get install -y \
    libgl1-mesa-glx \
    libxcb-xinerama0 \
    x11-utils \
    && rm -rf /var/lib/apt/lists/*

COPY analyze.py /app/
COPY requirements.txt /app/

# Install pymol using conda
RUN conda install -y -c conda-forge -c schrodinger pymol-bundle && \
    pip install -r requirements.txt

ENTRYPOINT ["/bin/bash"]
