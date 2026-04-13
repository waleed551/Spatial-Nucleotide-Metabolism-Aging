# Spatial-Nucleotide-Metabolism-Aging

This repository contains the computational workflows and analysis scripts used in our study:

**“Integrative multi-omics reveals spatial nucleotide metabolic reprogramming in the ageing brain”**

The project integrates spatial transcriptomics, metabolic flux inference, spatial metabolomics, machine learning, and entropy-based systems analysis to characterize age-dependent nucleotide metabolic remodeling across the mouse brain.

---

## Overview

Brain ageing is associated with region-specific metabolic reprogramming.  
In this study, we combine:

- Spatial transcriptomics
- scFEA-based metabolic flux inference
- AFADESI-MSI spatial metabolomics
- Random forest machine learning
- Entropy-based multi-omics integration

to investigate nucleotide metabolic network remodeling in the ageing mouse brain.

---

## Data Sources

### Public datasets
- Spatial transcriptomics: GEO accession **GSE212903**
- Bulk RNA-seq differential expression: GEO accession **GSE212336**

### Experimental datasets
- AFADESI-MSI spatial metabolomics from mouse brain sections:
  - 6 months
  - 15 months
  - 21 months

---

## Repository Structure

```bash
Spatial-Nucleotide-Metabolism-Aging/
│
├── Figure1/        # Transcriptomics and DEG-metabolite network analysis
├── Figure2%4/        # scFEA metabolic flux analysis & Thalamus multi-omics integration analysis
├── Figure3/        # Spatial metabolomics + machine learning classification
├── Figure5/        # Entropy-based systems analysis
│
├── data/           # Input and processed datasets
├── results/        # Output figures and tables
├── docs/           # Workflow diagrams and supplementary documents
│
├── README.md
├── requirements.txt
├── environment.yml
└── .gitignore
