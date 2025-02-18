# Protein Structure Analysis Tools

A toolkit for analyzing protein structures with focus on cancer mutations. This project provides tools to:
- Download and analyze PDB structures for specific genes
- Calculate structure quality scores
- Identify oncogenic hotspots using GENIE cohort data
- Align multiple structures using PyMOL

## Quick Start

1. Clone the repository:
```bash
git clone https://github.com/campandtell/protein-structure-tools
cd protein-structure-tools
```

2. Build and run Docker container:
```bash
docker build -t protein-tools .
docker run -it \
    -v $(pwd):/app \
    -e DISPLAY=$DISPLAY \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    --net=host \
    protein-tools
```

## Tutorial: Analyzing EGFR Structures

### 1. Analyze Available Structures
First, let's analyze all available EGFR structures:
```bash
python analyze.py analyze EGFR
```
This creates `results/egfr_structures.csv` with quality scores for each structure.

### 2. Download High-Quality Structures
Download structures with quality score ≥ 30:
```bash
python analyze.py download results/egfr_structures.csv --min-quality 30
```
Structures will be downloaded to the `pdbs` directory.

### 3. Calculate Oncogenic Hotspots
Using GENIE cohort data, identify mutation hotspots:
```bash
python get_hotspots.py --csv ./genie_data/EGFR.csv
```

### 4. Align Structures in PyMOL
Navigate to the PDB directory and align all structures:
```bash
cd pdbs
pymol ../align_all.py
```

## Command Reference

### analyze.py
```bash
# Analyze structures for a gene
python analyze.py analyze GENE_NAME

# Download structures above quality threshold
python analyze.py download RESULTS_CSV --min-quality SCORE

# Fetch specific PDB files
python analyze.py fetch PDB_ID [PDB_ID ...]
```

### get_hotspots.py
```bash
python get_hotspots.py --csv PATH_TO_GENIE_CSV
```

### align_all.py
Run within PyMOL after loading PDB files.

## Requirements
- Docker
- X11 for PyMOL visualization

## Directory Structure
```
protein-structure-tools/
├── analyze.py          # Main analysis script
├── align_all.py        # Structure alignment script
├── get_hotspots.py     # Mutation hotspot analysis
├── genie_data/         # GENIE cohort data
│   └── EGFR.csv
├── results/            # Analysis results (created)
└── pdbs/              # Downloaded structures (created)
```

## Note
Make sure to allow X11 connections before running PyMOL:
```bash
xhost +local:docker
```
