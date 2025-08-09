# PancancerCP
R code for figures in the pan-cancer cellular programs paper

## Overview
This repository contains R scripts and example data for processing, integration, visualization, and analysis of single-cell RNA-seq and spatial transcriptomic data.
The workflow is based on Seurat and related R packages, and reproduces the main computational analyses and figures reported in our manuscript.

## Software dependencies
R ≥ 4.2.0

R packages:

Seurat (v4.2.0 or above)

DoubletFinder (v2.0.3)

ComplexHeatmap (v2.6.2)

clusterProfiler (v4.3.2.991)

UCell (v1.99.7)

monocle (v2.18.0)

CytoTRACE (v0.3.3)

CellPhoneDB (v5.0)

ggplot2 (v3.4.4)

ggsankey (v0.0.9)

pheatmap (v1.0.12)

dendextend (v1.15.2)

cowplot (v1.1.2)

patchwork (v1.1.3)

ggrepel (v0.9.3)

Non-standard hardware
None required; tested on standard laptops/desktops with ≥8 GB RAM (16 GB recommended for large datasets)

## Code Functionality Description
The provided R scripts implement a complete workflow for:

Pre-processing and quality control of scRNA-seq data

Integration of multiple datasets to remove batch effects

Identification of cell types using clustering and marker gene analysis

Visualization of cell populations in reduced dimensions

Mapping scRNA-seq annotations to spatial transcriptomics data

Functional analysis via pseudotime inference, CNV detection, enrichment analysis

Inference of cell–cell communication networks

The workflow is modular, allowing users to run the full pipeline or specific modules independently.

## Data Availability for Demo
The demo dataset (pbmc_demo.rds) is included in the data/ folder.
Full raw data used in the manuscript can be downloaded from public databases:

GEO: GSE139829, GSE158803, GSE178318, GSE140312, GSE149614, GSE138709, GSE179994, GSE123902, GSE169246, GSE176078, GSE197177, GSE154778, GSE162708, GSE163558, GSE225857, GSE272362

NGDC: OMIX001073

## Citation:
In preparation
