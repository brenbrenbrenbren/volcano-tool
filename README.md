# Volcano Plot Generator & Cross-Dataset DEG Analysis Tool

A comprehensive Streamlit application for visualizing and comparing differential gene expression (DEG) data across multiple datasets.

## Core Question This Tool Answers

> "Is gene X differentially expressed across multiple experimental conditions, and if so, is it regulated in the same direction?"

## Features

### 🚀 Quick Workflow
1. **Upload** multiple Excel files (all sheets loaded by default)
2. **Dashboard** shows overview statistics and gene distribution
3. **UpSet** identifies overlapping hits across datasets
4. **Concordance** shows genes regulated in same direction (robust hits!)
5. **Report** exports top recurring genes

### ⌨️ Keyboard Shortcuts
- **`/`** - Focus gene search
- **`Ctrl+L`** - Toggle dark mode
- **`?`** - Show keyboard shortcuts help
- **`Ctrl+K`** - Command palette
- **Quick Presets** - One-click threshold adjustment (Stringent/Moderate/Lenient)

### 📊 Key Analysis Tabs
- **Dashboard**: Overview statistics, gene distribution charts, quick gene lists
- **Concordance** (NEW!): Identify genes regulated in same direction across datasets
  - Concordant Up/Down (robust biological hits)
  - Discordant (potential artifacts or context-dependent)
  - Downloadable gene lists by concordance type
- **UpSet Plots**: Visualize overlaps between multiple datasets with hover gene names
- **Report**: Find genes consistently significant across multiple datasets

### 📁 Data Import
- **Multi-file Excel support**: Upload multiple Excel files at once
- **All sheets selected by default**: No manual clicking
- **Checkbox selection**: Quick toggle with Select All/Clear All buttons
- **Flexible column detection**: Auto-recognizes Gene, log2FC, FC, pvalue, padj, FDR, etc.
- **Intensity data processing**: Auto-calculate log2FC + p-values from proteomics data
- **Gene ID conversion**: UniProt, Ensembl, Entrez → Gene symbols via MyGene API
- **FC to log2FC conversion**: Automatic transformation

### 📈 Visualizations
- **Volcano Plots**: Interactive plots with customizable thresholds and multi-dataset grid
- **Heatmaps**: Clustered expression with optional Z-score normalization
- **FC vs FC Scatter**: Pairwise dataset comparison with concordance stats
- **Batch Comparison Matrix**: Pairwise correlation/concordance across all datasets
- **Gene Bar Plots**: Single gene expression across all datasets
- **Pathway Analysis**: Visualize pathway gene expression patterns

### 💾 Export & Session
- **Figure Export**: Download plots as PNG or PDF
- **Gene List Export**: Copy/download gene lists (comma-separated or newlines)
- **Session Save/Load**: Save entire analysis state to JSON

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

### Finding Robust Biological Hits
1. **Start with Concordance Tab**: Genes that are concordant (same direction) across multiple datasets are most likely real biological effects
2. **Use Dashboard**: Quickly see which genes appear in ALL datasets or 2+ datasets
3. **Check UpSet Plot**: Visualize exact overlap patterns between datasets
4. **Validate with Report Tab**: Top recurring genes sorted by frequency

### Quick Actions
- **Press `?`** to see all keyboard shortcuts
- **Use Quick Presets** (Stringent/Moderate/Lenient) to rapidly adjust thresholds
- **Click gene lists** in code blocks to select all, then `Ctrl+C` to copy
- **Press `/`** to jump to gene search

### Data Analysis Tips
- **Z-score normalization**: Use in Heatmap when comparing different assay types (e.g., proximity vs RNA-seq)
- **Concordant genes**: Focus on these for robust, reproducible findings
- **Discordant genes**: May indicate technical artifacts or context-dependent regulation
- **Intensity data**: Tool auto-detects proteomics data and calculates log2FC + p-values for you

### Performance Tips
- Load all sheets at once using **Select All** button
- Use **threshold presets** instead of manual slider adjustment
- **Export gene lists** for downstream pathway analysis

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions welcome! Please open an issue or pull request.
