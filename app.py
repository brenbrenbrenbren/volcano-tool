"""
Volcano Plot Generator & Cross-Dataset DEG Analysis Tool

A comprehensive Streamlit application for visualizing and comparing
differential gene expression data across multiple datasets.
"""

import streamlit as st
import pandas as pd
import numpy as np
import json
import io
import base64
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple
from itertools import combinations

# Import local modules
from core import DataLoader, DataPreprocessor, DEGAnalyzer, run_deseq2, run_simple_de
from plotting import (
    create_volcano_plot, create_multi_volcano,
    create_upset_plot, get_intersection_genes,
    create_heatmap,
    create_fc_scatter, get_concordance_stats,
    create_gene_barplot, create_pathway_barplot
)
from utils import GeneMapper, PathwayDatabase, Exporter

# Page config
st.set_page_config(
    page_title="Volcano Plot Generator",
    page_icon="🌋",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Add keyboard shortcut hints via custom CSS
st.markdown("""
<style>
    /* Improve dark mode visibility */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    .stTabs [data-baseweb="tab"] {
        padding: 8px 16px;
        font-weight: 500;
    }

    /* Keyboard shortcut hints */
    .shortcut-hint {
        font-size: 10px;
        color: #888;
        margin-left: 4px;
    }

    /* Gene code blocks styling */
    .stCode {
        background-color: rgba(0,0,0,0.1);
        border-radius: 4px;
        padding: 8px;
    }

    /* Quick copy styling */
    .gene-list-box {
        font-family: monospace;
        font-size: 12px;
        user-select: all;
        cursor: pointer;
    }

    /* Checkbox grid for dataset selection */
    .dataset-checkbox-grid {
        display: grid;
        grid-template-columns: repeat(auto-fill, minmax(150px, 1fr));
        gap: 4px;
    }
</style>
""", unsafe_allow_html=True)

# Initialize session state
if 'datasets' not in st.session_state:
    st.session_state.datasets = {}
if 'raw_datasets' not in st.session_state:
    st.session_state.raw_datasets = {}  # Store original data for ID conversion
if 'loader' not in st.session_state:
    st.session_state.loader = DataLoader()
if 'preprocessor' not in st.session_state:
    st.session_state.preprocessor = DataPreprocessor()
if 'pathway_db' not in st.session_state:
    st.session_state.pathway_db = PathwayDatabase()
if 'exporter' not in st.session_state:
    st.session_state.exporter = Exporter()
if 'gene_mapper' not in st.session_state:
    st.session_state.gene_mapper = GeneMapper()
if 'dark_mode' not in st.session_state:
    st.session_state.dark_mode = False
if 'id_type_detected' not in st.session_state:
    st.session_state.id_type_detected = None


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_figure_download_buttons(fig, filename_base: str, key_suffix: str = ""):
    """Create download buttons for PNG and PDF export of a Plotly figure."""
    col1, col2 = st.columns(2)

    with col1:
        try:
            png_bytes = fig.to_image(format="png", width=1200, height=800, scale=2)
            st.download_button(
                "📥 Download PNG",
                data=png_bytes,
                file_name=f"{filename_base}.png",
                mime="image/png",
                key=f"png_{filename_base}_{key_suffix}"
            )
        except Exception as e:
            st.caption(f"PNG export requires kaleido: `pip install kaleido`")

    with col2:
        try:
            pdf_bytes = fig.to_image(format="pdf", width=1200, height=800, scale=2)
            st.download_button(
                "📥 Download PDF",
                data=pdf_bytes,
                file_name=f"{filename_base}.pdf",
                mime="application/pdf",
                key=f"pdf_{filename_base}_{key_suffix}"
            )
        except Exception as e:
            st.caption(f"PDF export requires kaleido: `pip install kaleido`")


def save_session() -> bytes:
    """Save current session state to JSON bytes."""
    session_data = {
        'version': '1.0',
        'timestamp': datetime.now().isoformat(),
        'datasets': {
            name: df.to_dict(orient='records')
            for name, df in st.session_state.datasets.items()
        },
        'settings': {
            'dark_mode': st.session_state.dark_mode,
        }
    }
    return json.dumps(session_data, indent=2).encode('utf-8')


def load_session(uploaded_file) -> bool:
    """Load session state from uploaded JSON file."""
    try:
        session_data = json.load(uploaded_file)

        # Restore datasets
        st.session_state.datasets = {
            name: pd.DataFrame(records)
            for name, records in session_data.get('datasets', {}).items()
        }

        # Restore settings
        settings = session_data.get('settings', {})
        st.session_state.dark_mode = settings.get('dark_mode', False)

        return True
    except Exception as e:
        st.error(f"Failed to load session: {str(e)}")
        return False


def create_batch_comparison_matrix(
    datasets: Dict[str, pd.DataFrame],
    log2fc_threshold: float,
    pvalue_threshold: float,
    metric: str = 'correlation'
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create pairwise comparison matrix for all dataset pairs.

    Returns:
        Tuple of (metric_matrix, details_dataframe)
    """
    dataset_names = list(datasets.keys())
    n = len(dataset_names)

    # Initialize matrix
    matrix = pd.DataFrame(
        np.eye(n),  # 1s on diagonal
        index=dataset_names,
        columns=dataset_names
    )

    details = []

    for i, name_a in enumerate(dataset_names):
        for j, name_b in enumerate(dataset_names):
            if i >= j:
                continue

            try:
                stats = get_concordance_stats(
                    datasets, name_a, name_b,
                    log2fc_threshold, pvalue_threshold
                )

                if metric == 'correlation':
                    value = stats['correlation']
                elif metric == 'concordance':
                    value = stats['concordance_rate']
                elif metric == 'overlap':
                    value = stats['sig_in_both'] / max(stats['sig_in_a'], stats['sig_in_b'], 1)
                else:
                    value = stats['correlation']

                matrix.loc[name_a, name_b] = value
                matrix.loc[name_b, name_a] = value

                details.append({
                    'Dataset A': name_a,
                    'Dataset B': name_b,
                    'Correlation': stats['correlation'],
                    'Concordant': stats['concordant'],
                    'Discordant': stats['discordant'],
                    'Concordance Rate': stats['concordance_rate'],
                    'Sig in Both': stats['sig_in_both'],
                    'Total Overlap': stats['total_genes']
                })
            except Exception:
                matrix.loc[name_a, name_b] = np.nan
                matrix.loc[name_b, name_a] = np.nan

    return matrix, pd.DataFrame(details)


def find_top_recurring_genes(
    datasets: Dict[str, pd.DataFrame],
    log2fc_threshold: float,
    pvalue_threshold: float,
    use_padj: bool = True,
    min_datasets: int = 1,
    top_n: int = 20,
    direction: str = 'both'
) -> pd.DataFrame:
    """
    Find genes that appear as significant across multiple datasets.

    Args:
        datasets: Dictionary of dataset name to DataFrame
        log2fc_threshold: Threshold for fold change significance
        pvalue_threshold: Threshold for p-value significance
        use_padj: Use adjusted p-value
        min_datasets: Minimum number of datasets gene must be significant in
        top_n: Number of top genes to return
        direction: 'both', 'up', or 'down'

    Returns:
        DataFrame with gene statistics across datasets
    """
    pval_col = 'padj' if use_padj else 'pvalue'
    gene_stats = {}

    for name, df in datasets.items():
        # Get significant genes
        if direction == 'up':
            sig_mask = (df['log2FC'] >= log2fc_threshold) & (df[pval_col] <= pvalue_threshold)
        elif direction == 'down':
            sig_mask = (df['log2FC'] <= -log2fc_threshold) & (df[pval_col] <= pvalue_threshold)
        else:
            sig_mask = (df['log2FC'].abs() >= log2fc_threshold) & (df[pval_col] <= pvalue_threshold)

        sig_genes = df[sig_mask]

        for _, row in sig_genes.iterrows():
            gene = row['Gene']
            if gene not in gene_stats:
                gene_stats[gene] = {
                    'Gene': gene,
                    'Datasets_Significant': [],
                    'log2FC_values': [],
                    'pvalue_values': [],
                    'Directions': []
                }

            gene_stats[gene]['Datasets_Significant'].append(name)
            gene_stats[gene]['log2FC_values'].append(row['log2FC'])
            gene_stats[gene]['pvalue_values'].append(row.get(pval_col, row.get('pvalue', 1)))
            gene_stats[gene]['Directions'].append('Up' if row['log2FC'] > 0 else 'Down')

    # Convert to DataFrame
    results = []
    for gene, stats in gene_stats.items():
        n_datasets = len(stats['Datasets_Significant'])
        if n_datasets >= min_datasets:
            # Check direction consistency
            up_count = stats['Directions'].count('Up')
            down_count = stats['Directions'].count('Down')

            if up_count == n_datasets:
                consistency = 'Consistently Up'
            elif down_count == n_datasets:
                consistency = 'Consistently Down'
            else:
                consistency = f'Mixed ({up_count} up, {down_count} down)'

            results.append({
                'Gene': gene,
                'Datasets_Significant': n_datasets,
                'Dataset_Names': ', '.join(stats['Datasets_Significant']),
                'Mean_log2FC': np.mean(stats['log2FC_values']),
                'Max_log2FC': max(stats['log2FC_values'], key=abs),
                'Min_pvalue': min(stats['pvalue_values']),
                'Direction_Consistency': consistency,
                'All_log2FC': stats['log2FC_values']
            })

    if not results:
        return pd.DataFrame()

    result_df = pd.DataFrame(results)

    # Sort by number of datasets (descending), then by mean fold change magnitude
    result_df = result_df.sort_values(
        by=['Datasets_Significant', 'Mean_log2FC'],
        ascending=[False, False],
        key=lambda x: x.abs() if x.name == 'Mean_log2FC' else x
    )

    return result_df.head(top_n)


# =============================================================================
# MAIN APPLICATION
# =============================================================================

def main():
    """Main application entry point."""

    # Sidebar for global settings
    with st.sidebar:
        st.title("🌋 Volcano Tool")
        st.markdown("---")

        # Theme toggle
        st.session_state.dark_mode = st.toggle("Dark Mode", value=st.session_state.dark_mode)

        st.markdown("### Global Thresholds")
        log2fc_threshold = st.slider(
            "|log2FC| threshold",
            min_value=0.0,
            max_value=5.0,
            value=1.0,
            step=0.1,
            help="Absolute log2 fold change threshold for significance"
        )

        pvalue_threshold = st.select_slider(
            "P-value threshold",
            options=[0.001, 0.01, 0.05, 0.1],
            value=0.05
        )

        use_padj = st.checkbox("Use adjusted p-value", value=True)

        direction = st.selectbox(
            "DEG direction",
            options=['both', 'up', 'down'],
            format_func=lambda x: {'both': 'Both', 'up': 'Upregulated only', 'down': 'Downregulated only'}[x]
        )

        st.markdown("---")

        # Session management
        st.markdown("### Session")

        col_save, col_load = st.columns(2)
        with col_save:
            if st.session_state.datasets:
                session_bytes = save_session()
                st.download_button(
                    "💾 Save",
                    data=session_bytes,
                    file_name=f"volcano_session_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
                    mime="application/json",
                    help="Save current session"
                )

        with col_load:
            uploaded_session = st.file_uploader(
                "Load",
                type=['json'],
                label_visibility="collapsed",
                help="Load a saved session"
            )
            if uploaded_session:
                if load_session(uploaded_session):
                    st.success("Session loaded!")
                    st.rerun()

        st.markdown("---")

        # Dataset summary
        if st.session_state.datasets:
            st.markdown("### Loaded Datasets")
            for name in st.session_state.datasets:
                df = st.session_state.datasets[name]
                n_genes = len(df)
                n_sig = len(df[(df['log2FC'].abs() >= log2fc_threshold) & (df['padj'] <= pvalue_threshold)])
                st.caption(f"**{name}**: {n_genes} genes ({n_sig} sig)")

        st.markdown("---")

        # Keyboard shortcuts help
        with st.expander("⌨️ Tips"):
            st.markdown("""
            **Quick Actions:**
            - Click gene lists to select all (then Ctrl+C)
            - Use checkboxes to toggle datasets
            - `Select All` / `Clear All` for bulk selection

            **Navigation:**
            - Use tabs to switch between views
            - Sidebar controls apply globally
            """)

    # Main content tabs
    tabs = st.tabs([
        "📁 Data Import",
        "🌋 Volcano Plots",
        "🔗 Overlap (UpSet)",
        "🗺️ Heatmap",
        "📊 FC vs FC",
        "📈 Batch Compare",
        "🔍 Gene Search",
        "🛤️ Pathways",
        "📋 Report",
        "💾 Export"
    ])

    # Tab 1: Data Import
    with tabs[0]:
        data_import_tab(log2fc_threshold, pvalue_threshold)

    # Tab 2: Volcano Plots
    with tabs[1]:
        volcano_tab(log2fc_threshold, pvalue_threshold, use_padj)

    # Tab 3: UpSet Plot
    with tabs[2]:
        upset_tab(log2fc_threshold, pvalue_threshold, use_padj, direction)

    # Tab 4: Heatmap
    with tabs[3]:
        heatmap_tab(log2fc_threshold, pvalue_threshold)

    # Tab 5: FC vs FC Scatter
    with tabs[4]:
        scatter_tab(log2fc_threshold, pvalue_threshold, use_padj)

    # Tab 6: Batch Comparison (NEW)
    with tabs[5]:
        batch_comparison_tab(log2fc_threshold, pvalue_threshold)

    # Tab 7: Gene Search
    with tabs[6]:
        gene_search_tab(log2fc_threshold, pvalue_threshold, use_padj)

    # Tab 8: Pathway Analysis
    with tabs[7]:
        pathway_tab(log2fc_threshold, pvalue_threshold)

    # Tab 9: Report
    with tabs[8]:
        report_tab(log2fc_threshold, pvalue_threshold, use_padj)

    # Tab 10: Export
    with tabs[9]:
        export_tab(log2fc_threshold, pvalue_threshold)


def data_import_tab(log2fc_threshold: float, pvalue_threshold: float):
    """Data import and configuration tab with gene ID conversion."""
    st.header("Data Import & Configuration")

    col1, col2 = st.columns([1, 1])

    with col1:
        st.subheader("Upload File")
        uploaded_file = st.file_uploader(
            "Choose Excel or CSV file",
            type=['xlsx', 'xls', 'csv'],
            help="Upload an Excel file with multiple sheets or a CSV file"
        )

        if uploaded_file:
            try:
                sheet_names = st.session_state.loader.load_file_from_buffer(uploaded_file)
                st.success(f"Loaded file with {len(sheet_names)} sheet(s)")

                # Sheet selector with search
                search_term = st.text_input("🔍 Search sheets", placeholder="Type to filter...")
                filtered_sheets = st.session_state.loader.filter_sheet_names(search_term)

                selected_sheets = st.multiselect(
                    "Select sheets to analyze",
                    options=filtered_sheets,
                    default=filtered_sheets[:min(5, len(filtered_sheets))]
                )

                if selected_sheets:
                    st.markdown("---")
                    st.subheader("Column Mapping")

                    # Get first sheet for column detection
                    sample_df = st.session_state.loader.get_sheet(selected_sheets[0])
                    col_info = st.session_state.loader.get_column_info(sample_df)
                    data_type = st.session_state.loader.detect_data_type(sample_df)

                    st.info(f"Detected data type: **{data_type.replace('_', ' ').title()}**")

                    all_cols = col_info['all_columns']

                    gene_col = st.selectbox(
                        "Gene/Protein column",
                        options=all_cols,
                        index=all_cols.index(col_info['gene_col']) if col_info['gene_col'] in all_cols else 0
                    )

                    # Gene ID type detection and conversion option
                    st.markdown("---")
                    st.subheader("🔄 Gene ID Conversion")

                    # Detect ID type from sample
                    sample_ids = sample_df[gene_col].dropna().astype(str).tolist()[:100]
                    detected_type = st.session_state.gene_mapper.detect_id_type(sample_ids)
                    st.session_state.id_type_detected = detected_type

                    st.write(f"Detected ID type: **{detected_type.upper()}**")

                    convert_ids = st.checkbox(
                        "Convert to gene symbols",
                        value=(detected_type not in ['symbol', 'unknown']),
                        help="Convert UniProt, Ensembl, or Entrez IDs to gene symbols"
                    )

                    if convert_ids:
                        species = st.selectbox(
                            "Species",
                            options=['human', 'mouse'],
                            help="Select species for ID conversion"
                        )
                    else:
                        species = 'human'

                    st.markdown("---")

                    if data_type == 'deg_results':
                        # Determine the best FC column and whether log transform is needed
                        fc_col_detected = col_info.get('log2fc_col') or col_info.get('fc_col')
                        needs_log_transform = col_info.get('needs_log_transform', False)

                        # Show info about detected columns
                        if col_info.get('log2fc_col'):
                            st.success(f"Detected log2FC column: **{col_info['log2fc_col']}**")
                        elif col_info.get('fc_col'):
                            st.warning(f"Detected FC column: **{col_info['fc_col']}** (will convert to log2)")

                        log2fc_col = st.selectbox(
                            "Fold Change column (log2FC or FC)",
                            options=all_cols,
                            index=all_cols.index(fc_col_detected) if fc_col_detected and fc_col_detected in all_cols else 0
                        )

                        # Check if selected column needs log transform
                        selected_col_lower = log2fc_col.lower().replace(' ', '').replace('_', '').replace('-', '').replace('.', '')
                        is_log_col = any(x in selected_col_lower for x in ['log2', 'log', 'lfc'])
                        needs_log_transform_final = not is_log_col

                        if needs_log_transform_final:
                            st.info("📊 Selected column appears to be linear FC - will convert to log2FC")

                        pvalue_col = st.selectbox(
                            "P-value column",
                            options=all_cols,
                            index=all_cols.index(col_info['pvalue_col']) if col_info['pvalue_col'] and col_info['pvalue_col'] in all_cols else 0
                        )

                        padj_options = ['(Use p-value)'] + all_cols
                        padj_default = 0
                        if col_info['padj_col'] and col_info['padj_col'] in all_cols:
                            padj_default = padj_options.index(col_info['padj_col'])

                        padj_col = st.selectbox(
                            "Adjusted p-value column",
                            options=padj_options,
                            index=padj_default
                        )
                        if padj_col == '(Use p-value)':
                            padj_col = None

                        if st.button("🚀 Process Datasets", type="primary"):
                            process_deg_datasets(
                                selected_sheets, gene_col, log2fc_col, pvalue_col, padj_col,
                                convert_ids=convert_ids, species=species,
                                needs_log_transform=needs_log_transform_final
                            )

                    else:  # raw_counts
                        st.markdown("**Configure sample groups for DE analysis:**")

                        numeric_cols = sample_df.select_dtypes(include=['number']).columns.tolist()

                        group_a = st.multiselect("Control samples", options=numeric_cols)
                        group_b = st.multiselect("Treatment samples", options=numeric_cols)

                        if group_a and group_b:
                            if st.button("🧬 Run Differential Expression", type="primary"):
                                run_de_analysis(selected_sheets, gene_col, group_a, group_b,
                                               convert_ids=convert_ids, species=species)

            except Exception as e:
                st.error(f"Error loading file: {str(e)}")

    with col2:
        st.subheader("Loaded Datasets Preview")

        if st.session_state.datasets:
            preview_name = st.selectbox(
                "Select dataset to preview",
                options=list(st.session_state.datasets.keys())
            )

            if preview_name:
                df = st.session_state.datasets[preview_name]

                # Stats
                n_sig = len(df[(df['log2FC'].abs() >= log2fc_threshold) & (df['padj'] <= pvalue_threshold)])
                n_up = len(df[(df['log2FC'] >= log2fc_threshold) & (df['padj'] <= pvalue_threshold)])
                n_down = len(df[(df['log2FC'] <= -log2fc_threshold) & (df['padj'] <= pvalue_threshold)])

                col_a, col_b, col_c = st.columns(3)
                col_a.metric("Total Genes", len(df))
                col_b.metric("Upregulated", n_up)
                col_c.metric("Downregulated", n_down)

                st.dataframe(df.head(100), use_container_width=True)
        else:
            st.info("No datasets loaded yet. Upload a file to get started.")


def process_deg_datasets(sheets: List[str], gene_col: str, log2fc_col: str,
                         pvalue_col: str, padj_col: Optional[str],
                         convert_ids: bool = False, species: str = 'human',
                         needs_log_transform: bool = False):
    """Process DEG result datasets with optional ID conversion and FC transformation."""
    preprocessor = st.session_state.preprocessor
    mapper = st.session_state.gene_mapper

    with st.spinner("Processing datasets..."):
        for sheet in sheets:
            try:
                df = st.session_state.loader.get_sheet(sheet)

                # Verify required columns exist
                missing_cols = []
                if gene_col not in df.columns:
                    missing_cols.append(f"gene column '{gene_col}'")
                if log2fc_col not in df.columns:
                    missing_cols.append(f"FC column '{log2fc_col}'")
                if pvalue_col not in df.columns:
                    missing_cols.append(f"p-value column '{pvalue_col}'")

                if missing_cols:
                    st.warning(f"Sheet '{sheet}' missing: {', '.join(missing_cols)}. Skipping.")
                    continue

                # Convert IDs if requested
                if convert_ids:
                    with st.spinner(f"Converting gene IDs for {sheet}..."):
                        try:
                            df = mapper.convert_to_symbols(df, gene_col, species)
                            # Use converted symbols
                            actual_gene_col = 'Gene_Symbol'
                        except Exception as e:
                            st.warning(f"ID conversion failed for {sheet}: {str(e)}. Using original IDs.")
                            actual_gene_col = gene_col
                else:
                    actual_gene_col = gene_col

                standardized = preprocessor.standardize_dataset(
                    df, actual_gene_col, log2fc_col, pvalue_col, padj_col,
                    needs_log_transform=needs_log_transform
                )
                st.session_state.datasets[sheet] = standardized
            except Exception as e:
                st.warning(f"Failed to process {sheet}: {str(e)}")

    st.success(f"Processed {len(st.session_state.datasets)} datasets")
    st.rerun()


def run_de_analysis(sheets: List[str], gene_col: str, group_a: List[str], group_b: List[str],
                   convert_ids: bool = False, species: str = 'human'):
    """Run differential expression analysis on count data."""
    mapper = st.session_state.gene_mapper

    with st.spinner("Running differential expression analysis..."):
        for sheet in sheets:
            try:
                df = st.session_state.loader.get_sheet(sheet)

                try:
                    result = run_deseq2(df, gene_col, group_a, group_b)
                except ImportError:
                    st.warning(f"PyDESeq2 not available for {sheet}, using simple DE")
                    result = run_simple_de(df, gene_col, group_a, group_b)

                # Convert IDs if requested
                if convert_ids:
                    try:
                        result = mapper.convert_to_symbols(result, 'Gene', species)
                        result['Gene'] = result['Gene_Symbol']
                        result = result.drop(columns=['Gene_Symbol'])
                    except Exception as e:
                        st.warning(f"ID conversion failed for {sheet}: {str(e)}")

                st.session_state.datasets[sheet] = result
            except Exception as e:
                st.warning(f"Failed to analyze {sheet}: {str(e)}")

    st.success(f"Analyzed {len(st.session_state.datasets)} datasets")
    st.rerun()


def volcano_tab(log2fc_threshold: float, pvalue_threshold: float, use_padj: bool):
    """Volcano plot visualization tab with export buttons."""
    st.header("Volcano Plots")

    if not st.session_state.datasets:
        st.info("Please load datasets in the Data Import tab first.")
        return

    col1, col2 = st.columns([3, 1])

    with col2:
        st.subheader("Options")

        selected_dataset = st.selectbox(
            "Select dataset",
            options=list(st.session_state.datasets.keys())
        )

        show_labels = st.checkbox("Show gene labels", value=True)
        top_n = st.slider("Top N genes to label", 0, 50, 10)

        highlight_input = st.text_area(
            "Highlight genes (comma-separated)",
            placeholder="CTNNB1, TP53, MYC"
        )
        highlight_genes = [g.strip() for g in highlight_input.split(',') if g.strip()]

    with col1:
        if selected_dataset:
            df = st.session_state.datasets[selected_dataset]

            fig = create_volcano_plot(
                df,
                title=selected_dataset,
                log2fc_threshold=log2fc_threshold,
                pvalue_threshold=pvalue_threshold,
                use_padj=use_padj,
                highlight_genes=highlight_genes if highlight_genes else None,
                dark_mode=st.session_state.dark_mode,
                show_labels=show_labels,
                top_n_labels=top_n
            )

            st.plotly_chart(fig, use_container_width=True)

            # Export buttons
            get_figure_download_buttons(fig, f"volcano_{selected_dataset}", "single")

            # Show multi-volcano option
            if len(st.session_state.datasets) > 1:
                if st.checkbox("Show all datasets in grid"):
                    multi_fig = create_multi_volcano(
                        st.session_state.datasets,
                        log2fc_threshold=log2fc_threshold,
                        pvalue_threshold=pvalue_threshold,
                        use_padj=use_padj,
                        dark_mode=st.session_state.dark_mode
                    )
                    st.plotly_chart(multi_fig, use_container_width=True)
                    get_figure_download_buttons(multi_fig, "volcano_grid", "multi")


def upset_tab(log2fc_threshold: float, pvalue_threshold: float, use_padj: bool, direction: str):
    """UpSet plot for overlap analysis with export."""
    st.header("Cross-Dataset Overlap Analysis")

    if len(st.session_state.datasets) < 2:
        st.info("Please load at least 2 datasets to analyze overlaps.")
        return

    col1, col2 = st.columns([3, 1])

    with col2:
        st.subheader("Dataset Selection")

        # Toggle all buttons
        col_all, col_none = st.columns(2)
        with col_all:
            if st.button("Select All", key="upset_select_all", use_container_width=True):
                st.session_state.upset_selected = list(st.session_state.datasets.keys())
        with col_none:
            if st.button("Clear All", key="upset_clear_all", use_container_width=True):
                st.session_state.upset_selected = []

        # Initialize selection state
        if 'upset_selected' not in st.session_state:
            st.session_state.upset_selected = list(st.session_state.datasets.keys())

        # Dataset checkboxes
        st.caption(f"**{len(st.session_state.upset_selected)}/{len(st.session_state.datasets)}** selected")

        selected_datasets = []
        for ds_name in st.session_state.datasets.keys():
            is_selected = st.checkbox(
                ds_name,
                value=ds_name in st.session_state.upset_selected,
                key=f"upset_cb_{ds_name}"
            )
            if is_selected:
                selected_datasets.append(ds_name)

        st.session_state.upset_selected = selected_datasets

        st.markdown("---")
        st.subheader("Plot Options")
        min_size = st.slider("Min intersection size", 1, 50, 1, key="upset_min_size")
        max_subsets = st.slider("Max intersections", 10, 100, 40, key="upset_max_subsets")

    with col1:
        if len(selected_datasets) >= 2:
            # Get DEG sets
            analyzer = DEGAnalyzer(log2fc_threshold, pvalue_threshold, use_padj)
            filtered_datasets = {k: st.session_state.datasets[k] for k in selected_datasets}
            deg_sets = analyzer.get_deg_sets(filtered_datasets, direction)

            # Create UpSet plot (now returns intersection data too)
            fig, intersection_data = create_upset_plot(
                deg_sets,
                min_subset_size=min_size,
                max_subsets=max_subsets,
                dark_mode=st.session_state.dark_mode,
                show_gene_labels_threshold=10
            )

            st.plotly_chart(fig, use_container_width=True)
            get_figure_download_buttons(fig, "upset_plot", "upset")

            # Quick summary of key intersections
            st.markdown("---")
            st.subheader("Key Intersections")

            # Show genes in ALL datasets first
            all_datasets_combo = tuple(sorted(selected_datasets))
            genes_in_all = get_intersection_genes(deg_sets, selected_datasets, exclusive=False)

            if genes_in_all:
                st.success(f"**{len(genes_in_all)} genes** significant in ALL {len(selected_datasets)} datasets")
                genes_list = sorted(list(genes_in_all))
                if len(genes_list) <= 20:
                    st.write(", ".join(genes_list))
                else:
                    st.write(", ".join(genes_list[:20]) + f"... (+{len(genes_list)-20} more)")

                # Copy button
                col_copy, col_download = st.columns(2)
                with col_copy:
                    st.code(", ".join(genes_list), language=None)
                with col_download:
                    st.download_button(
                        "📋 Download Gene List",
                        "\n".join(genes_list),
                        f"genes_in_all_{len(selected_datasets)}_datasets.txt",
                        "text/plain",
                        key="download_all_genes"
                    )
            else:
                st.info("No genes significant in all selected datasets")

            # Show small intersections prominently
            small_intersections = [d for d in intersection_data if d['count'] <= 10 and d['count'] > 0]
            if small_intersections:
                st.markdown("---")
                st.subheader("Small Intersections (≤10 genes)")
                for inter in small_intersections[:5]:
                    with st.expander(f"**{inter['count']} genes** in {len(inter['datasets'])} datasets: {', '.join(inter['datasets'][:3])}{'...' if len(inter['datasets']) > 3 else ''}"):
                        st.write(", ".join(inter['genes']))
                        st.code(", ".join(inter['genes']), language=None)

            # Intersection explorer
            st.markdown("---")
            st.subheader("Custom Intersection Explorer")

            intersection_sets = st.multiselect(
                "Select datasets to find intersection",
                options=selected_datasets,
                key="intersection_selector"
            )

            if intersection_sets:
                exclusive = st.checkbox("Exclusive intersection only", value=True)
                genes = get_intersection_genes(deg_sets, intersection_sets, exclusive)

                st.write(f"**{len(genes)} genes** in this intersection:")
                if genes:
                    genes_sorted = sorted(list(genes))
                    st.code(", ".join(genes_sorted), language=None)
                    st.download_button(
                        "📋 Download",
                        "\n".join(genes_sorted),
                        "intersection_genes.txt",
                        "text/plain",
                        key="download_custom_intersection"
                    )
        else:
            st.warning("Please select at least 2 datasets")


def heatmap_tab(log2fc_threshold: float, pvalue_threshold: float):
    """Heatmap visualization tab with export."""
    st.header("Expression Heatmap")

    if not st.session_state.datasets:
        st.info("Please load datasets in the Data Import tab first.")
        return

    col1, col2 = st.columns([3, 1])

    with col2:
        st.subheader("Options")

        filter_sig = st.checkbox("Filter to significant genes", value=True)
        n_datasets = len(st.session_state.datasets)
        if n_datasets > 1:
            min_datasets = st.slider("Minimum datasets significant in", 1, n_datasets, 1, key="heatmap_min_ds")
        else:
            min_datasets = 1
            st.caption("Min datasets: 1 (only one dataset loaded)")

        st.markdown("---")
        st.subheader("Display")

        show_values = st.checkbox(
            "Show values on cells",
            value=False,
            help="Display log2FC values as text on each cell"
        )

        normalize_rows = st.checkbox(
            "Normalize rows (Z-score)",
            value=False,
            help="Z-score normalize each gene across datasets. Useful when comparing different assay types (e.g., proximity vs sequencing)"
        )

        if normalize_rows:
            st.info("Z-score normalization: shows relative expression pattern across datasets, not absolute FC values")

        st.markdown("---")
        st.subheader("Clustering")

        cluster_genes = st.checkbox("Cluster genes", value=True)
        cluster_datasets = st.checkbox("Cluster datasets", value=False)

        st.markdown("---")
        max_genes = st.slider("Maximum genes to show", 20, 200, 100)

        custom_genes = st.text_area(
            "Custom gene list (optional)",
            placeholder="Enter genes comma-separated"
        )
        genes_list = [g.strip() for g in custom_genes.split(',') if g.strip()] if custom_genes else None

    with col1:
        fig = create_heatmap(
            st.session_state.datasets,
            genes=genes_list,
            cluster_genes=cluster_genes,
            cluster_datasets=cluster_datasets,
            log2fc_threshold=log2fc_threshold,
            pvalue_threshold=pvalue_threshold,
            filter_significant=filter_sig,
            min_datasets=min_datasets,
            dark_mode=st.session_state.dark_mode,
            max_genes=max_genes,
            show_values=show_values,
            normalize_rows=normalize_rows
        )

        st.plotly_chart(fig, use_container_width=True)
        get_figure_download_buttons(fig, "heatmap", "heatmap")


def scatter_tab(log2fc_threshold: float, pvalue_threshold: float, use_padj: bool):
    """FC vs FC scatter plot tab with export."""
    st.header("Fold Change Comparison")

    if len(st.session_state.datasets) < 2:
        st.info("Please load at least 2 datasets to compare.")
        return

    col1, col2 = st.columns([3, 1])

    dataset_names = list(st.session_state.datasets.keys())

    with col2:
        st.subheader("Options")

        dataset_a = st.selectbox("Dataset A (x-axis)", options=dataset_names, index=0)
        dataset_b = st.selectbox(
            "Dataset B (y-axis)",
            options=dataset_names,
            index=min(1, len(dataset_names)-1)
        )

        show_diagonal = st.checkbox("Show diagonal (y=x)", value=True)

        highlight_input = st.text_area(
            "Highlight genes",
            placeholder="CTNNB1, TP53"
        )
        highlight_genes = [g.strip() for g in highlight_input.split(',') if g.strip()]

    with col1:
        if dataset_a != dataset_b:
            fig = create_fc_scatter(
                st.session_state.datasets,
                dataset_a,
                dataset_b,
                log2fc_threshold=log2fc_threshold,
                pvalue_threshold=pvalue_threshold,
                use_padj=use_padj,
                highlight_genes=highlight_genes if highlight_genes else None,
                show_diagonal=show_diagonal,
                dark_mode=st.session_state.dark_mode
            )

            st.plotly_chart(fig, use_container_width=True)
            get_figure_download_buttons(fig, f"scatter_{dataset_a}_vs_{dataset_b}", "scatter")

            # Show concordance stats
            stats = get_concordance_stats(
                st.session_state.datasets,
                dataset_a,
                dataset_b,
                log2fc_threshold,
                pvalue_threshold
            )

            st.markdown("### Concordance Statistics")
            col_a, col_b, col_c, col_d = st.columns(4)
            col_a.metric("Overlap Genes", stats['total_genes'])
            col_b.metric("Concordant", stats['concordant'])
            col_c.metric("Discordant", stats['discordant'])
            col_d.metric("Correlation", f"{stats['correlation']:.3f}")
        else:
            st.warning("Please select two different datasets to compare.")


def batch_comparison_tab(log2fc_threshold: float, pvalue_threshold: float):
    """Batch comparison matrix view for all dataset pairs."""
    st.header("Batch Dataset Comparison")

    if len(st.session_state.datasets) < 2:
        st.info("Please load at least 2 datasets to compare.")
        return

    col1, col2 = st.columns([3, 1])

    with col2:
        st.subheader("Options")

        metric = st.selectbox(
            "Comparison metric",
            options=['correlation', 'concordance', 'overlap'],
            format_func=lambda x: {
                'correlation': 'Pearson Correlation',
                'concordance': 'Concordance Rate',
                'overlap': 'DEG Overlap Ratio'
            }[x]
        )

        selected_datasets = st.multiselect(
            "Select datasets",
            options=list(st.session_state.datasets.keys()),
            default=list(st.session_state.datasets.keys())
        )

    with col1:
        if len(selected_datasets) >= 2:
            filtered_datasets = {k: st.session_state.datasets[k] for k in selected_datasets}

            with st.spinner("Computing pairwise comparisons..."):
                matrix, details = create_batch_comparison_matrix(
                    filtered_datasets,
                    log2fc_threshold,
                    pvalue_threshold,
                    metric
                )

            # Create heatmap of comparison matrix
            import plotly.graph_objects as go

            if st.session_state.dark_mode:
                colorscale = [[0, '#4ECDC4'], [0.5, '#2D2D2D'], [1, '#FF6B6B']]
                paper_bgcolor = '#1E1E1E'
                font_color = '#FFFFFF'
            else:
                colorscale = [[0, '#3498DB'], [0.5, '#FFFFFF'], [1, '#E74C3C']]
                paper_bgcolor = '#FFFFFF'
                font_color = '#2C3E50'

            fig = go.Figure(data=go.Heatmap(
                z=matrix.values,
                x=matrix.columns.tolist(),
                y=matrix.index.tolist(),
                colorscale=colorscale,
                zmin=-1 if metric == 'correlation' else 0,
                zmax=1,
                text=np.round(matrix.values, 3),
                texttemplate='%{text}',
                textfont={"size": 10},
                hovertemplate=(
                    '%{y} vs %{x}<br>'
                    f'{metric.title()}: %{{z:.3f}}<br>'
                    '<extra></extra>'
                )
            ))

            fig.update_layout(
                title=f"Pairwise {metric.title()} Matrix",
                paper_bgcolor=paper_bgcolor,
                font=dict(color=font_color),
                height=max(400, len(selected_datasets) * 50),
                xaxis=dict(tickangle=45),
            )

            st.plotly_chart(fig, use_container_width=True)
            get_figure_download_buttons(fig, f"batch_comparison_{metric}", "batch")

            # Show details table
            st.markdown("### Pairwise Details")
            st.dataframe(
                details.style.format({
                    'Correlation': '{:.3f}',
                    'Concordance Rate': '{:.2%}'
                }),
                use_container_width=True
            )

            # Download details
            csv = details.to_csv(index=False)
            st.download_button(
                "📥 Download Comparison Table (CSV)",
                csv,
                "batch_comparison_details.csv",
                "text/csv"
            )
        else:
            st.warning("Please select at least 2 datasets.")


def gene_search_tab(log2fc_threshold: float, pvalue_threshold: float, use_padj: bool):
    """Single gene search tab with export."""
    st.header("Gene Search")

    if not st.session_state.datasets:
        st.info("Please load datasets in the Data Import tab first.")
        return

    col1, col2 = st.columns([3, 1])

    with col2:
        st.subheader("Search")

        gene_name = st.text_input(
            "Enter gene name",
            placeholder="e.g., CTNNB1"
        ).strip()

        if gene_name:
            # Quick info
            analyzer = DEGAnalyzer(log2fc_threshold, pvalue_threshold, use_padj)
            presence = analyzer.get_gene_presence(st.session_state.datasets, gene_name)

            sig_count = presence['Significant'].sum()
            st.metric("Significant in", f"{sig_count}/{len(st.session_state.datasets)} datasets")

    with col1:
        if gene_name:
            fig = create_gene_barplot(
                st.session_state.datasets,
                gene_name,
                log2fc_threshold=log2fc_threshold,
                pvalue_threshold=pvalue_threshold,
                use_padj=use_padj,
                dark_mode=st.session_state.dark_mode
            )

            st.plotly_chart(fig, use_container_width=True)
            get_figure_download_buttons(fig, f"gene_{gene_name}", "gene")

            # Detail table
            st.markdown("### Detailed Values")
            st.dataframe(presence, use_container_width=True)
        else:
            st.info("Enter a gene name to search across all datasets.")


def pathway_tab(log2fc_threshold: float, pvalue_threshold: float):
    """Pathway analysis tab with export."""
    st.header("Pathway Analysis")

    if not st.session_state.datasets:
        st.info("Please load datasets in the Data Import tab first.")
        return

    col1, col2 = st.columns([3, 1])

    pathway_db = st.session_state.pathway_db

    with col2:
        st.subheader("Select Pathway")

        # List available pathways
        all_pathways = pathway_db.list_all_pathways()

        selected_pathway = st.selectbox(
            "Choose pathway",
            options=all_pathways,
            format_func=lambda x: x.replace('_', ' ').title()
        )

        show_all = st.checkbox("Show all genes (not just significant)", value=False)
        max_genes = st.slider("Max genes to display", 10, 100, 30)

        # Show pathway info
        if selected_pathway:
            genes = pathway_db.get_pathway_genes(selected_pathway)
            st.info(f"**{len(genes)} genes** in this pathway")

    with col1:
        if selected_pathway:
            genes = pathway_db.get_pathway_genes(selected_pathway)

            fig = create_pathway_barplot(
                st.session_state.datasets,
                genes,
                selected_pathway.replace('_', ' ').title(),
                log2fc_threshold=log2fc_threshold,
                pvalue_threshold=pvalue_threshold,
                dark_mode=st.session_state.dark_mode,
                show_all_genes=show_all,
                max_genes=max_genes
            )

            st.plotly_chart(fig, use_container_width=True)
            get_figure_download_buttons(fig, f"pathway_{selected_pathway}", "pathway")

            # Show overlap statistics
            st.markdown("### Pathway Overlap by Dataset")

            overlap_data = []
            for name, df in st.session_state.datasets.items():
                sig_genes = set(df[
                    (df['log2FC'].abs() >= log2fc_threshold) &
                    (df['padj'] <= pvalue_threshold)
                ]['Gene'].str.upper())

                overlap = pathway_db.get_pathway_overlap(
                    {g.upper() for g in genes},
                    sig_genes
                )

                overlap_data.append({
                    'Dataset': name,
                    'DEGs in Pathway': overlap['overlap_count'],
                    'Overlap %': f"{overlap['overlap_fraction']*100:.1f}%"
                })

            st.dataframe(pd.DataFrame(overlap_data), use_container_width=True)


def report_tab(log2fc_threshold: float, pvalue_threshold: float, use_padj: bool):
    """Generate a report of top recurring genes across datasets."""
    st.header("Top Genes Report")

    if not st.session_state.datasets:
        st.info("Please load datasets in the Data Import tab first.")
        return

    st.markdown("""
    This report identifies genes that are **consistently differentially expressed**
    across multiple datasets. Genes appearing in more datasets are prioritized.
    """)

    col1, col2 = st.columns([3, 1])

    with col2:
        st.subheader("Report Options")

        # Dataset selection
        selected_datasets = st.multiselect(
            "Select datasets to include",
            options=list(st.session_state.datasets.keys()),
            default=list(st.session_state.datasets.keys()),
            help="Choose which datasets to analyze"
        )

        # Number of top genes
        top_n = st.selectbox(
            "Number of top genes",
            options=[5, 10, 20, 30, 50, 100],
            index=2,
            help="How many top genes to include in the report"
        )

        # Minimum datasets
        n_selected = len(selected_datasets)
        if n_selected > 1:
            min_datasets = st.slider(
                "Minimum datasets",
                min_value=1,
                max_value=n_selected,
                value=max(1, n_selected // 2),
                help="Minimum number of datasets a gene must be significant in"
            )
        else:
            min_datasets = 1

        # Direction filter
        direction = st.selectbox(
            "Expression direction",
            options=['both', 'up', 'down'],
            format_func=lambda x: {
                'both': 'Both (up & down)',
                'up': 'Upregulated only',
                'down': 'Downregulated only'
            }[x]
        )

        generate_report = st.button("Generate Report", type="primary", use_container_width=True)

    with col1:
        if len(selected_datasets) < 1:
            st.warning("Please select at least 1 dataset.")
            return

        if generate_report or 'last_report' in st.session_state:
            # Store report parameters to maintain state
            if generate_report:
                st.session_state.report_params = {
                    'selected_datasets': selected_datasets,
                    'top_n': top_n,
                    'min_datasets': min_datasets,
                    'direction': direction
                }

            params = st.session_state.get('report_params', {
                'selected_datasets': selected_datasets,
                'top_n': top_n,
                'min_datasets': min_datasets,
                'direction': direction
            })

            filtered_datasets = {k: st.session_state.datasets[k] for k in params['selected_datasets']}

            with st.spinner("Analyzing genes across datasets..."):
                report_df = find_top_recurring_genes(
                    filtered_datasets,
                    log2fc_threshold,
                    pvalue_threshold,
                    use_padj=use_padj,
                    min_datasets=params['min_datasets'],
                    top_n=params['top_n'],
                    direction=params['direction']
                )

            if report_df.empty:
                st.warning("No genes found matching the criteria. Try adjusting thresholds or minimum datasets.")
                return

            st.session_state.last_report = report_df

            # Summary metrics
            st.subheader("Summary")
            total_datasets = len(params['selected_datasets'])

            col_m1, col_m2, col_m3, col_m4 = st.columns(4)

            # Genes in all datasets
            genes_in_all = len(report_df[report_df['Datasets_Significant'] == total_datasets])
            col_m1.metric("In All Datasets", genes_in_all)

            # Consistently up
            consistently_up = len(report_df[report_df['Direction_Consistency'] == 'Consistently Up'])
            col_m2.metric("Consistently Up", consistently_up)

            # Consistently down
            consistently_down = len(report_df[report_df['Direction_Consistency'] == 'Consistently Down'])
            col_m3.metric("Consistently Down", consistently_down)

            # Mixed direction
            mixed = len(report_df) - consistently_up - consistently_down
            col_m4.metric("Mixed Direction", mixed)

            st.markdown("---")

            # Main report table
            st.subheader(f"Top {len(report_df)} Recurring DEGs")

            # Format for display
            display_df = report_df[[
                'Gene', 'Datasets_Significant', 'Direction_Consistency',
                'Mean_log2FC', 'Max_log2FC', 'Min_pvalue', 'Dataset_Names'
            ]].copy()

            display_df.columns = [
                'Gene', '# Datasets', 'Direction', 'Mean log2FC',
                'Max log2FC', 'Min p-value', 'Datasets'
            ]

            # Style the dataframe
            def highlight_direction(val):
                if 'Up' in str(val) and 'Down' not in str(val):
                    return 'background-color: rgba(231, 76, 60, 0.3)'
                elif 'Down' in str(val) and 'Up' not in str(val):
                    return 'background-color: rgba(52, 152, 219, 0.3)'
                elif 'Mixed' in str(val):
                    return 'background-color: rgba(241, 196, 15, 0.3)'
                return ''

            styled_df = display_df.style.applymap(
                highlight_direction, subset=['Direction']
            ).format({
                'Mean log2FC': '{:.3f}',
                'Max log2FC': '{:.3f}',
                'Min p-value': '{:.2e}'
            })

            st.dataframe(styled_df, use_container_width=True, height=400)

            # Visualization
            st.markdown("---")
            st.subheader("Visualization")

            viz_type = st.radio(
                "Chart type",
                options=['bar', 'heatmap'],
                format_func=lambda x: 'Bar Chart' if x == 'bar' else 'Heatmap',
                horizontal=True
            )

            import plotly.graph_objects as go
            import plotly.express as px

            if st.session_state.dark_mode:
                paper_bgcolor = '#1E1E1E'
                plot_bgcolor = '#2D2D2D'
                font_color = '#FFFFFF'
            else:
                paper_bgcolor = '#FFFFFF'
                plot_bgcolor = '#FAFAFA'
                font_color = '#2C3E50'

            if viz_type == 'bar':
                # Bar chart showing genes and their dataset count
                fig = go.Figure()

                # Color by direction consistency
                colors = []
                for _, row in display_df.iterrows():
                    if 'Consistently Up' in row['Direction']:
                        colors.append('#E74C3C' if not st.session_state.dark_mode else '#FF6B6B')
                    elif 'Consistently Down' in row['Direction']:
                        colors.append('#3498DB' if not st.session_state.dark_mode else '#4ECDC4')
                    else:
                        colors.append('#F39C12' if not st.session_state.dark_mode else '#FFD93D')

                fig.add_trace(go.Bar(
                    x=display_df['Gene'],
                    y=display_df['# Datasets'],
                    marker_color=colors,
                    text=display_df['# Datasets'],
                    textposition='outside',
                    hovertemplate=(
                        '<b>%{x}</b><br>'
                        'Significant in: %{y} datasets<br>'
                        '<extra></extra>'
                    )
                ))

                # Add line for total datasets
                fig.add_hline(
                    y=total_datasets,
                    line_dash="dash",
                    line_color="gray",
                    annotation_text=f"All {total_datasets} datasets"
                )

                fig.update_layout(
                    title=f"Top {len(display_df)} Recurring DEGs",
                    xaxis_title="Gene",
                    yaxis_title="Number of Datasets",
                    paper_bgcolor=paper_bgcolor,
                    plot_bgcolor=plot_bgcolor,
                    font=dict(color=font_color),
                    xaxis=dict(tickangle=45),
                    height=500
                )

                st.plotly_chart(fig, use_container_width=True)
                get_figure_download_buttons(fig, "top_genes_bar", "report_bar")

            else:  # Heatmap
                # Build heatmap matrix: genes x datasets showing log2FC
                genes_to_show = report_df['Gene'].tolist()[:min(30, len(report_df))]

                heatmap_data = []
                for gene in genes_to_show:
                    row_data = {'Gene': gene}
                    for ds_name in params['selected_datasets']:
                        df = filtered_datasets[ds_name]
                        gene_row = df[df['Gene'] == gene]
                        if len(gene_row) > 0:
                            row_data[ds_name] = gene_row.iloc[0]['log2FC']
                        else:
                            row_data[ds_name] = np.nan
                    heatmap_data.append(row_data)

                heatmap_df = pd.DataFrame(heatmap_data).set_index('Gene')

                # Create heatmap
                max_val = min(5, heatmap_df.abs().max().max())

                if st.session_state.dark_mode:
                    colorscale = [[0, '#4ECDC4'], [0.5, '#2D2D2D'], [1, '#FF6B6B']]
                else:
                    colorscale = [[0, '#3498DB'], [0.5, '#FFFFFF'], [1, '#E74C3C']]

                fig = go.Figure(data=go.Heatmap(
                    z=heatmap_df.values,
                    x=heatmap_df.columns.tolist(),
                    y=heatmap_df.index.tolist(),
                    colorscale=colorscale,
                    zmin=-max_val,
                    zmax=max_val,
                    colorbar=dict(title=dict(text="log2FC", font=dict(color=font_color))),
                    hovertemplate=(
                        '<b>%{y}</b><br>'
                        'Dataset: %{x}<br>'
                        'log2FC: %{z:.3f}<br>'
                        '<extra></extra>'
                    )
                ))

                fig.update_layout(
                    title=f"Top Recurring DEGs - Expression Pattern",
                    xaxis_title="Dataset",
                    yaxis_title="Gene",
                    paper_bgcolor=paper_bgcolor,
                    plot_bgcolor=plot_bgcolor,
                    font=dict(color=font_color),
                    xaxis=dict(tickangle=45),
                    yaxis=dict(autorange="reversed"),
                    height=max(400, len(genes_to_show) * 20 + 150)
                )

                st.plotly_chart(fig, use_container_width=True)
                get_figure_download_buttons(fig, "top_genes_heatmap", "report_heatmap")

            # Download options
            st.markdown("---")
            st.subheader("Download Report")

            col_d1, col_d2 = st.columns(2)

            with col_d1:
                # CSV download
                csv_df = report_df.drop(columns=['All_log2FC'])
                csv = csv_df.to_csv(index=False)
                st.download_button(
                    "📥 Download Report (CSV)",
                    csv,
                    f"top_genes_report_{datetime.now().strftime('%Y%m%d')}.csv",
                    "text/csv",
                    use_container_width=True
                )

            with col_d2:
                # Gene list only (for pathway analysis tools)
                gene_list = '\n'.join(report_df['Gene'].tolist())
                st.download_button(
                    "📥 Download Gene List (TXT)",
                    gene_list,
                    f"top_genes_list_{datetime.now().strftime('%Y%m%d')}.txt",
                    "text/plain",
                    use_container_width=True
                )

        else:
            st.info("Configure options and click **Generate Report** to identify top recurring genes.")

            # Show quick preview
            if len(selected_datasets) >= 2:
                st.markdown("### Quick Preview")
                st.caption("Genes appearing in the most datasets:")

                quick_preview = find_top_recurring_genes(
                    {k: st.session_state.datasets[k] for k in selected_datasets},
                    log2fc_threshold,
                    pvalue_threshold,
                    use_padj=use_padj,
                    min_datasets=1,
                    top_n=5,
                    direction='both'
                )

                if not quick_preview.empty:
                    for _, row in quick_preview.iterrows():
                        st.write(f"**{row['Gene']}** - {row['Datasets_Significant']}/{len(selected_datasets)} datasets ({row['Direction_Consistency']})")


def export_tab(log2fc_threshold: float, pvalue_threshold: float):
    """Export tab with session save/load."""
    st.header("Export Data & Figures")

    if not st.session_state.datasets:
        st.info("Please load datasets in the Data Import tab first.")
        return

    exporter = st.session_state.exporter

    # Session management section
    st.subheader("💾 Session Management")

    col_s1, col_s2 = st.columns(2)

    with col_s1:
        st.markdown("**Save Current Session**")
        st.caption("Save all loaded datasets and settings to a JSON file for later use.")

        session_bytes = save_session()
        st.download_button(
            "💾 Download Session File",
            data=session_bytes,
            file_name=f"volcano_session_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
            mime="application/json",
            use_container_width=True
        )

    with col_s2:
        st.markdown("**Load Previous Session**")
        st.caption("Restore a previously saved session.")

        uploaded_session = st.file_uploader(
            "Choose session file",
            type=['json'],
            key="export_tab_session_loader"
        )

        if uploaded_session:
            if st.button("🔄 Restore Session", use_container_width=True):
                if load_session(uploaded_session):
                    st.success("Session restored successfully!")
                    st.rerun()

    st.markdown("---")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Export Summary Report")

        if st.button("Generate Summary Report"):
            report = exporter.create_deg_report(
                st.session_state.datasets,
                {'log2fc': log2fc_threshold, 'pvalue': pvalue_threshold}
            )

            st.dataframe(report, use_container_width=True)

            csv = report.to_csv(index=False)
            st.download_button(
                "📥 Download Summary (CSV)",
                csv,
                "deg_summary.csv",
                "text/csv"
            )

    with col2:
        st.subheader("Export Full Gene Table")

        if st.button("Generate Full Export"):
            full_export = exporter.create_full_export(
                st.session_state.datasets,
                {'log2fc': log2fc_threshold, 'pvalue': pvalue_threshold}
            )

            st.write(f"**{len(full_export)} genes** across all datasets")
            st.dataframe(full_export.head(100), use_container_width=True)

            csv = full_export.to_csv(index=False)
            st.download_button(
                "📥 Download Full Table (CSV)",
                csv,
                "deg_full_export.csv",
                "text/csv"
            )

    st.markdown("---")
    st.subheader("Export Individual Dataset")

    export_dataset = st.selectbox(
        "Select dataset",
        options=list(st.session_state.datasets.keys()),
        key="individual_export_selector"
    )

    if export_dataset:
        df = st.session_state.datasets[export_dataset]

        col_e1, col_e2 = st.columns(2)

        with col_e1:
            csv = df.to_csv(index=False)
            st.download_button(
                f"📥 Download {export_dataset} (CSV)",
                csv,
                f"{export_dataset}.csv",
                "text/csv",
                use_container_width=True
            )

        with col_e2:
            # Excel export
            buffer = io.BytesIO()
            df.to_excel(buffer, index=False, engine='openpyxl')
            buffer.seek(0)
            st.download_button(
                f"📥 Download {export_dataset} (Excel)",
                buffer.getvalue(),
                f"{export_dataset}.xlsx",
                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                use_container_width=True
            )


if __name__ == "__main__":
    main()
