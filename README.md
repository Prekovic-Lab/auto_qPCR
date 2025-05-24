# Prekovic Lab qPCR Analysis Portal

This repository provides the **Prekovic Lab qPCR Analysis App**, an intuitive and interactive online tool designed to streamline quantitative PCR (qPCR) data analysis. The app enables users to rapidly normalize expression data using custom-selected housekeeping genes and methods (individual or geometric mean normalization), visualize amplification and melt curves for quality control, and generate publication-quality bar plots with customizable colors and replicate-level details. Each plot can be downloaded as vectorized PDFs, and normalized data can be exported as CSV for further analysis.

## Key Features:
- **Interactive data normalization:** Choose housekeeping genes and normalization methods.
- **Customizable visualization:** Select colors for experimental conditions and inspect individual replicate points interactively.
- **Quality control tools:** Visualize melt curves and amplification data per well.
- **Export functionality:** Download publication-quality vectorized plots (PDF) and normalized data tables (CSV).

## Input requirements:
- Excel file (.xlsx format) from qPCR runs (compatible with standard Applied Biosystems export structure).

## 📌 Naming Conditions & Replicates:
For accurate replicate grouping, samples should follow the naming convention:

- Condtion_name floowed by space and the replicate number; example -> "MYC-OE 1"

## Access the App:
[👉 autoqpcr.streamlit.app](https://autoqpcr.streamlit.app/)

---
Maintained by the [Prekovic Lab](https://www.prekovic-lab.org).  
Contact: [s.prekovic@umcutrecht.nl](mailto:s.prekovic@umcutrecht.nl)
