# Volcano Plot Generator & Cross-Dataset DEG Analysis Tool

A comprehensive Streamlit application for visualizing and comparing differential gene expression (DEG) data across multiple datasets.

## Core Question This Tool Answers

> "Is gene X differentially expressed across multiple experimental conditions, and if so, is it regulated in the same direction?"

## Features

### Data Import
- **Multi-sheet Excel support**: Load multiple datasets from a single Excel file
- **Flexible column detection**: Automatically recognizes various column naming conventions (Gene, Gene Symbol, Protein, log2FC, FC, pvalue, padj, FDR, etc.)
- **Gene ID conversion**: Convert UniProt, Ensembl, or Entrez IDs to gene symbols via MyGene API
- **FC to log2FC conversion**: Automatically converts linear fold change to log2 scale

### Visualizations
- **Volcano Plots**: Interactive plots with customizable thresholds, gene highlighting, and multi-dataset grid view
- **UpSet Plots**: Visualize overlaps between multiple datasets with gene annotations for small intersections
- **Heatmaps**: Clustered expression heatmaps with optional Z-score normalization for cross-assay comparison
- **FC vs FC Scatter**: Pairwise dataset comparison with concordance statistics
- **Gene Bar Plots**: Single gene expression across all datasets
- **Pathway Analysis**: Visualize pathway gene expression patterns

### Cross-Dataset Analysis
- **Batch Comparison Matrix**: Pairwise correlation/concordance across all datasets
- **Top Genes Report**: Find genes consistently significant across multiple datasets
- **Intersection Explorer**: Quickly identify genes shared across selected datasets

### Export & Session Management
- **Figure Export**: Download plots as PNG or PDF
- **Gene List Export**: Copy/download gene lists in comma-separated format
- **Session Save/Load**: Save analysis state to JSON for later use

## Installation

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/volcano-tool.git
cd volcano-tool

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# For figure export (PNG/PDF), also install:
pip install kaleido
```

## Usage

```bash
streamlit run app.py
```

Then open your browser to `http://localhost:8501`

### Input File Format

Your Excel/CSV file should contain DEG results with columns for:
- **Gene identifier**: Gene, Gene Symbol, Protein, UniProt, etc.
- **Fold change**: log2FC, log2FoldChange, logFC, FC, FoldChange, etc.
- **P-value**: pvalue, p-value, PValue, etc.
- **Adjusted p-value** (optional): padj, FDR, qvalue, adj.P.Val, etc.

Multiple sheets in an Excel file are treated as separate datasets/conditions.

## Project Structure

```
volcano_tool/
├── app.py                 # Main Streamlit application
├── requirements.txt       # Python dependencies
├── core/
│   ├── __init__.py
│   ├── data_loader.py    # File loading and column detection
│   ├── preprocessing.py  # Data standardization
│   ├── analysis.py       # DEG analysis utilities
│   └── deseq_runner.py   # DE analysis from counts
├── plotting/
│   ├── __init__.py
│   ├── volcano.py        # Volcano plot generation
│   ├── upset.py          # UpSet plot for overlaps
│   ├── heatmap.py        # Expression heatmaps
│   ├── scatter.py        # FC vs FC scatter plots
│   └── barplot.py        # Gene/pathway bar plots
├── utils/
│   ├── __init__.py
│   ├── gene_mapping.py   # Gene ID conversion
│   ├── pathways.py       # Pathway database
│   └── export.py         # Export utilities
└── config/
    └── pathways_cache.json
```

## Dependencies

- streamlit
- pandas
- numpy
- plotly
- scipy
- openpyxl
- requests (for gene ID conversion)
- kaleido (optional, for figure export)

## Tips for Effective Use

1. **Comparing different assay types**: Use Z-score normalization in the heatmap when comparing proximity data (smaller FC range) with RNA-seq (larger FC range)

2. **Finding robust hits**: Use the Report tab to find genes significant across multiple datasets - these are more likely to be real biological effects

3. **Quick gene lists**: Click on gene lists displayed in `code blocks` to select all, then Ctrl+C to copy

4. **Dataset management**: Use checkboxes to quickly toggle datasets on/off for analysis

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions welcome! Please open an issue or pull request.
