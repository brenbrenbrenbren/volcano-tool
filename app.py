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
from core import (
    DataLoader, DataPreprocessor, DEGAnalyzer,
    run_deseq2, run_simple_de, run_intensity_de, detect_sample_groups, auto_run_de
)
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

# Add keyboard shortcuts and custom CSS
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

    /* Modal styling */
    .modal-overlay {
        position: fixed;
        top: 0;
        left: 0;
        right: 0;
        bottom: 0;
        background: rgba(0,0,0,0.7);
        z-index: 9998;
        display: none;
    }
    .modal-content {
        position: fixed;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
        background: white;
        padding: 2rem;
        border-radius: 8px;
        z-index: 9999;
        min-width: 400px;
    }
</style>

<script>
// Keyboard shortcuts handler
document.addEventListener('keydown', function(e) {
    // Ignore if user is typing in an input
    if (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA') {
        return;
    }

    const key = e.key.toLowerCase();

    // Global shortcuts
    if (key === 'd' && e.ctrlKey) {
        e.preventDefault();
        // Toggle dark mode - we'll handle this via session state
        const darkToggle = document.querySelector('[data-testid="stCheckbox"]');
        if (darkToggle) darkToggle.click();
    }

    // Help modal
    if (key === '?') {
        e.preventDefault();
        showHelpModal();
    }

    // Gene search
    if (key === '/') {
        e.preventDefault();
        focusGeneSearch();
    }
});

function showHelpModal() {
    alert(`Keyboard Shortcuts:

Navigation:
  1-9    Jump to tab

Tools:
  /      Gene search
  Ctrl+D Dark mode toggle
  ?      Show this help

Tips:
  - Click gene lists to select all (Ctrl+C to copy)
  - Use threshold sliders in sidebar for quick filtering`);
}

function focusGeneSearch() {
    // Try to find and focus gene search input
    const geneInputs = Array.from(document.querySelectorAll('input')).filter(
        input => input.placeholder && input.placeholder.toLowerCase().includes('gene')
    );
    if (geneInputs.length > 0) {
        geneInputs[0].focus();
        geneInputs[0].scrollIntoView({behavior: 'smooth', block: 'center'});
    }
}
</script>
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

        # Quick preset buttons
        st.caption("Quick Presets:")
        preset_col1, preset_col2, preset_col3 = st.columns(3)

        with preset_col1:
            if st.button("Stringent", use_container_width=True, help="log2FC≥1.5, p≤0.01"):
                st.session_state.preset_fc = 1.5
                st.session_state.preset_p = 0.01
        with preset_col2:
            if st.button("Moderate", use_container_width=True, help="log2FC≥1.0, p≤0.05"):
                st.session_state.preset_fc = 1.0
                st.session_state.preset_p = 0.05
        with preset_col3:
            if st.button("Lenient", use_container_width=True, help="log2FC≥0.5, p≤0.1"):
                st.session_state.preset_fc = 0.5
                st.session_state.preset_p = 0.1

        st.markdown("---")

        # Initialize preset values if not set
        if 'preset_fc' not in st.session_state:
            st.session_state.preset_fc = 1.0
        if 'preset_p' not in st.session_state:
            st.session_state.preset_p = 0.05

        log2fc_threshold = st.slider(
            "|log2FC| threshold",
            min_value=0.0,
            max_value=5.0,
            value=st.session_state.preset_fc,
            step=0.1,
            help="Absolute log2 fold change threshold for significance"
        )

        pvalue_threshold = st.select_slider(
            "P-value threshold",
            options=[0.001, 0.01, 0.05, 0.1],
            value=st.session_state.preset_p
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
        with st.expander("⌨️ Keyboard Shortcuts"):
            st.markdown("""
            **Navigation:**
            - `1-9` - Jump to specific tab

            **Tools:**
            - `/` - Focus gene search
            - `Ctrl+D` - Toggle dark mode
            - `?` - Show keyboard help

            **Quick Actions:**
            - Click gene lists to select all (then `Ctrl+C`)
            - Use preset buttons for quick threshold changes
            - `Select All` / `Clear All` for bulk sheet selection
            """)

    # Main content tabs - reordered for optimal workflow
    tabs = st.tabs([
        "📁 Data Import",
        "📊 Dashboard",
        "🔗 UpSet",
        "🔄 Concordance",
        "📋 Report",
        "🔍 Gene Search",
        "🌋 Volcano",
        "🗺️ Heatmap",
        "📈 FC vs FC",
        "🧮 Batch Compare",
        "🛤️ Pathways",
        "💾 Export"
    ])

    # Tab 1: Data Import
    with tabs[0]:
        data_import_tab(log2fc_threshold, pvalue_threshold)

    # Tab 2: Dashboard (NEW) - Overview of all datasets
    with tabs[1]:
        dashboard_tab(log2fc_threshold, pvalue_threshold, use_padj)

    # Tab 3: UpSet - Key overlap visualization
    with tabs[2]:
        upset_tab(log2fc_threshold, pvalue_threshold, use_padj, direction)

    # Tab 4: Concordance (NEW) - Direction consistency
    with tabs[3]:
        concordance_tab(log2fc_threshold, pvalue_threshold, use_padj)

    # Tab 5: Report - Top recurring genes
    with tabs[4]:
        report_tab(log2fc_threshold, pvalue_threshold, use_padj)

    # Tab 6: Gene Search - Quick lookup
    with tabs[5]:
        gene_search_tab(log2fc_threshold, pvalue_threshold, use_padj)

    # Tab 7: Volcano Plots
    with tabs[6]:
        volcano_tab(log2fc_threshold, pvalue_threshold, use_padj)

    # Tab 8: Heatmap
    with tabs[7]:
        heatmap_tab(log2fc_threshold, pvalue_threshold)

    # Tab 9: FC vs FC Scatter
    with tabs[8]:
        scatter_tab(log2fc_threshold, pvalue_threshold, use_padj)

    # Tab 10: Batch Comparison
    with tabs[9]:
        batch_comparison_tab(log2fc_threshold, pvalue_threshold)

    # Tab 11: Pathway Analysis
    with tabs[10]:
        pathway_tab(log2fc_threshold, pvalue_threshold)

    # Tab 12: Export
    with tabs[11]:
        export_tab(log2fc_threshold, pvalue_threshold)


def data_import_tab(log2fc_threshold: float, pvalue_threshold: float):
    """Data import and configuration tab with gene ID conversion."""
    st.header("Data Import & Configuration")

    # Initialize multi-file storage
    if 'uploaded_files_data' not in st.session_state:
        st.session_state.uploaded_files_data = {}  # {filename: {excel_file, sheets, dataframes}}
    if 'all_available_sheets' not in st.session_state:
        st.session_state.all_available_sheets = []  # [(filename, sheetname), ...]
    if 'selected_sheet_keys' not in st.session_state:
        st.session_state.selected_sheet_keys = set()

    col1, col2 = st.columns([1, 1])

    with col1:
        st.subheader("Upload Files")
        uploaded_files = st.file_uploader(
            "Choose Excel or CSV files",
            type=['xlsx', 'xls', 'csv'],
            accept_multiple_files=True,
            help="Upload one or more Excel files (each can have multiple sheets)"
        )

        if uploaded_files:
            # Process all uploaded files
            all_sheets = []
            new_files_added = False

            for uploaded_file in uploaded_files:
                filename = uploaded_file.name

                # Skip if already loaded
                if filename not in st.session_state.uploaded_files_data:
                    try:
                        suffix = Path(filename).suffix.lower()
                        if suffix in ['.xlsx', '.xls']:
                            excel_file = pd.ExcelFile(uploaded_file)
                            sheet_names = excel_file.sheet_names
                            st.session_state.uploaded_files_data[filename] = {
                                'excel_file': excel_file,
                                'sheets': sheet_names,
                                'dataframes': {}
                            }
                        elif suffix == '.csv':
                            df = pd.read_csv(uploaded_file)
                            dataset_name = Path(filename).stem
                            st.session_state.uploaded_files_data[filename] = {
                                'excel_file': None,
                                'sheets': [dataset_name],
                                'dataframes': {dataset_name: df}
                            }
                        new_files_added = True
                    except Exception as e:
                        st.error(f"Error loading {filename}: {str(e)}")
                        continue

            # Build list of all available sheets
            all_sheets = []
            for filename, file_data in st.session_state.uploaded_files_data.items():
                for sheet in file_data['sheets']:
                    all_sheets.append((filename, sheet))

            st.session_state.all_available_sheets = all_sheets

            # Default: select ALL sheets on first load
            if new_files_added or not st.session_state.selected_sheet_keys:
                st.session_state.selected_sheet_keys = {f"{fn}::{sn}" for fn, sn in all_sheets}

            # Summary
            total_files = len(st.session_state.uploaded_files_data)
            total_sheets = len(all_sheets)
            st.success(f"**{total_files} file(s)** loaded with **{total_sheets} total sheet(s)**")

            # Select All / Clear All
            st.markdown("---")
            st.markdown("**Select sheets to analyze:**")
            col_all, col_none, col_count = st.columns([1, 1, 2])
            with col_all:
                if st.button("Select All", key="import_select_all", use_container_width=True):
                    st.session_state.selected_sheet_keys = {f"{fn}::{sn}" for fn, sn in all_sheets}
                    st.rerun()
            with col_none:
                if st.button("Clear All", key="import_clear_all", use_container_width=True):
                    st.session_state.selected_sheet_keys = set()
                    st.rerun()
            with col_count:
                n_selected = len(st.session_state.selected_sheet_keys)
                st.caption(f"**{n_selected}/{total_sheets}** selected")

            # Show checkboxes grouped by file
            selected_sheets_list = []  # [(filename, sheetname), ...]

            for filename, file_data in st.session_state.uploaded_files_data.items():
                sheets = file_data['sheets']
                if len(st.session_state.uploaded_files_data) > 1:
                    st.markdown(f"**{filename}** ({len(sheets)} sheets)")

                for sheet in sheets:
                    key = f"{filename}::{sheet}"
                    is_selected = st.checkbox(
                        sheet if len(st.session_state.uploaded_files_data) == 1 else f"  {sheet}",
                        value=key in st.session_state.selected_sheet_keys,
                        key=f"import_cb_{key}"
                    )
                    if is_selected:
                        selected_sheets_list.append((filename, sheet))
                        st.session_state.selected_sheet_keys.add(key)
                    else:
                        st.session_state.selected_sheet_keys.discard(key)

            if selected_sheets_list:
                st.markdown("---")
                st.subheader("Column Mapping")

                # Get first selected sheet for column detection
                first_file, first_sheet = selected_sheets_list[0]
                sample_df = get_sheet_data(first_file, first_sheet)

                if sample_df is not None:
                    col_info = st.session_state.loader.get_column_info(sample_df)
                    data_type = st.session_state.loader.detect_data_type(sample_df)

                    st.info(f"Detected data type: **{data_type.replace('_', ' ').title()}** (from first sheet)")

                    all_cols = col_info['all_columns']

                    gene_col = st.selectbox(
                        "Gene/Protein column",
                        options=all_cols,
                        index=all_cols.index(col_info['gene_col']) if col_info['gene_col'] in all_cols else 0
                    )

                    # Gene ID conversion (simplified)
                    convert_ids = st.checkbox(
                        "Convert gene IDs to symbols",
                        value=False,
                        help="Convert UniProt, Ensembl, or Entrez IDs to gene symbols"
                    )
                    species = 'human'
                    if convert_ids:
                        species = st.selectbox("Species", options=['human', 'mouse'])

                    st.markdown("---")

                    if data_type == 'deg_results':
                        # DEG results - need FC and pvalue columns
                        fc_col_detected = col_info.get('log2fc_col') or col_info.get('fc_col')

                        if col_info.get('log2fc_col'):
                            st.success(f"Detected log2FC column: **{col_info['log2fc_col']}**")
                        elif col_info.get('fc_col'):
                            st.warning(f"Detected FC column: **{col_info['fc_col']}** (will convert to log2)")

                        log2fc_col = st.selectbox(
                            "Fold Change column",
                            options=all_cols,
                            index=all_cols.index(fc_col_detected) if fc_col_detected and fc_col_detected in all_cols else 0
                        )

                        # Check if needs log transform
                        selected_col_lower = log2fc_col.lower().replace(' ', '').replace('_', '').replace('-', '').replace('.', '')
                        is_log_col = any(x in selected_col_lower for x in ['log2', 'log', 'lfc'])
                        needs_log_transform = not is_log_col

                        if needs_log_transform:
                            st.info("📊 Will convert to log2FC")

                        pvalue_col = st.selectbox(
                            "P-value column",
                            options=all_cols,
                            index=all_cols.index(col_info['pvalue_col']) if col_info['pvalue_col'] and col_info['pvalue_col'] in all_cols else 0
                        )

                        padj_options = ['(Use p-value)'] + all_cols
                        padj_default = 0
                        if col_info['padj_col'] and col_info['padj_col'] in all_cols:
                            padj_default = padj_options.index(col_info['padj_col'])
                        padj_col = st.selectbox("Adjusted p-value column", options=padj_options, index=padj_default)
                        if padj_col == '(Use p-value)':
                            padj_col = None

                        if st.button(f"🚀 Process {len(selected_sheets_list)} Dataset(s)", type="primary", use_container_width=True):
                            process_all_deg_datasets(
                                selected_sheets_list, gene_col, log2fc_col, pvalue_col, padj_col,
                                convert_ids=convert_ids, species=species,
                                needs_log_transform=needs_log_transform
                            )

                    elif data_type == 'intensity_data':
                        # Intensity data - need to calculate DE
                        st.subheader("🧬 Sample Group Configuration")
                        numeric_cols = sample_df.select_dtypes(include=['number']).columns.tolist()
                        detected_groups = detect_sample_groups(numeric_cols)

                        if len(detected_groups) >= 2:
                            st.success(f"Detected **{len(detected_groups)} groups**: {', '.join(detected_groups.keys())}")

                            with st.expander("View detected groups"):
                                for group_name, cols in detected_groups.items():
                                    st.write(f"**{group_name}**: {', '.join(cols)}")

                            group_names = list(detected_groups.keys())
                            control_group = st.selectbox(
                                "Control group",
                                options=group_names,
                                index=group_names.index('Control') if 'Control' in group_names else 0
                            )

                            other_groups = [g for g in group_names if g != control_group]
                            treatment_groups = st.multiselect("Treatment groups", options=other_groups, default=other_groups)

                            is_log = st.checkbox("Data is log2-transformed", value=True)

                            if treatment_groups:
                                if st.button(f"🧬 Calculate DE for {len(selected_sheets_list)} sheet(s)", type="primary", use_container_width=True):
                                    process_all_intensity_datasets(
                                        selected_sheets_list, gene_col,
                                        detected_groups, control_group, treatment_groups,
                                        is_log_transformed=is_log,
                                        convert_ids=convert_ids, species=species
                                    )
                        else:
                            st.warning("Could not auto-detect groups. Select manually:")
                            group_a = st.multiselect("Control samples", options=numeric_cols)
                            group_b = st.multiselect("Treatment samples", options=numeric_cols)
                            is_log = st.checkbox("Data is log2-transformed", value=True)

                            if group_a and group_b:
                                if st.button("🧬 Calculate DE", type="primary"):
                                    process_all_intensity_manual(
                                        selected_sheets_list, gene_col, group_a, group_b,
                                        is_log_transformed=is_log,
                                        convert_ids=convert_ids, species=species
                                    )

                    else:  # raw_counts
                        st.markdown("**Configure sample groups for DE analysis:**")

                        numeric_cols = sample_df.select_dtypes(include=['number']).columns.tolist()

                        group_a = st.multiselect("Control samples", options=numeric_cols)
                        group_b = st.multiselect("Treatment samples", options=numeric_cols)

                        if group_a and group_b:
                            if st.button("🧬 Run Differential Expression", type="primary"):
                                process_all_count_datasets(
                                    selected_sheets_list, gene_col, group_a, group_b,
                                    convert_ids=convert_ids, species=species
                                )

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


def get_sheet_data(filename: str, sheetname: str) -> Optional[pd.DataFrame]:
    """Get DataFrame for a specific sheet from uploaded files."""
    if filename not in st.session_state.uploaded_files_data:
        return None

    file_data = st.session_state.uploaded_files_data[filename]

    # Check if already loaded
    if sheetname in file_data['dataframes']:
        return file_data['dataframes'][sheetname]

    # Load from Excel file
    if file_data['excel_file'] is not None:
        try:
            df = pd.read_excel(file_data['excel_file'], sheet_name=sheetname)
            file_data['dataframes'][sheetname] = df
            return df
        except Exception as e:
            st.warning(f"Could not load {sheetname} from {filename}: {e}")
            return None

    return None


def process_all_deg_datasets(
    sheets_list: List[Tuple[str, str]],  # [(filename, sheetname), ...]
    gene_col: str, log2fc_col: str, pvalue_col: str, padj_col: Optional[str],
    convert_ids: bool = False, species: str = 'human',
    needs_log_transform: bool = False
):
    """Process DEG datasets from multiple files/sheets."""
    preprocessor = st.session_state.preprocessor
    mapper = st.session_state.gene_mapper

    processed = 0
    with st.spinner(f"Processing {len(sheets_list)} dataset(s)..."):
        for filename, sheetname in sheets_list:
            try:
                df = get_sheet_data(filename, sheetname)
                if df is None:
                    st.warning(f"Could not load {sheetname} from {filename}")
                    continue

                # Auto-detect columns for this sheet if they differ
                col_info = st.session_state.loader.get_column_info(df)

                # Use specified columns, or fall back to auto-detected
                actual_gene_col = gene_col if gene_col in df.columns else col_info.get('gene_col')
                actual_fc_col = log2fc_col if log2fc_col in df.columns else (col_info.get('log2fc_col') or col_info.get('fc_col'))
                actual_pval_col = pvalue_col if pvalue_col in df.columns else col_info.get('pvalue_col')
                actual_padj_col = padj_col if padj_col and padj_col in df.columns else col_info.get('padj_col')

                if not actual_gene_col or not actual_fc_col or not actual_pval_col:
                    st.warning(f"Sheet '{sheetname}' missing required columns. Skipping.")
                    continue

                # Check if FC column needs log transform
                fc_col_lower = actual_fc_col.lower().replace(' ', '').replace('_', '')
                sheet_needs_log = not any(x in fc_col_lower for x in ['log2', 'log', 'lfc'])

                # Convert IDs if requested
                if convert_ids:
                    try:
                        df = mapper.convert_to_symbols(df, actual_gene_col, species)
                        actual_gene_col = 'Gene_Symbol'
                    except:
                        pass

                standardized = preprocessor.standardize_dataset(
                    df, actual_gene_col, actual_fc_col, actual_pval_col, actual_padj_col,
                    needs_log_transform=sheet_needs_log or needs_log_transform
                )

                # Use sheetname as dataset name (add filename prefix if multiple files)
                if len(st.session_state.uploaded_files_data) > 1:
                    dataset_name = f"{Path(filename).stem}_{sheetname}"
                else:
                    dataset_name = sheetname

                st.session_state.datasets[dataset_name] = standardized
                processed += 1

            except Exception as e:
                st.warning(f"Failed to process {sheetname}: {str(e)}")

    st.success(f"Successfully processed **{processed}** dataset(s)")
    st.rerun()


def process_all_intensity_datasets(
    sheets_list: List[Tuple[str, str]],
    gene_col: str,
    detected_groups: dict, control_group: str, treatment_groups: List[str],
    is_log_transformed: bool = True,
    convert_ids: bool = False, species: str = 'human'
):
    """Process intensity data from multiple sheets."""
    mapper = st.session_state.gene_mapper
    control_cols = detected_groups[control_group]

    processed = 0
    with st.spinner(f"Calculating DE for {len(sheets_list)} sheet(s)..."):
        for filename, sheetname in sheets_list:
            try:
                df = get_sheet_data(filename, sheetname)
                if df is None:
                    continue

                # Auto-detect gene column if needed
                if gene_col not in df.columns:
                    col_info = st.session_state.loader.get_column_info(df)
                    actual_gene_col = col_info.get('gene_col')
                    if not actual_gene_col:
                        st.warning(f"Sheet '{sheetname}' missing gene column. Skipping.")
                        continue
                else:
                    actual_gene_col = gene_col

                for treatment in treatment_groups:
                    treatment_cols = detected_groups[treatment]

                    # Check columns exist
                    missing = [c for c in control_cols + treatment_cols if c not in df.columns]
                    if missing:
                        st.warning(f"Sheet '{sheetname}' missing columns: {missing[:3]}...")
                        continue

                    result = run_intensity_de(
                        df, actual_gene_col, control_cols, treatment_cols,
                        is_log_transformed=is_log_transformed
                    )

                    if convert_ids:
                        try:
                            result = mapper.convert_to_symbols(result, 'Gene', species)
                            result['Gene'] = result['Gene_Symbol']
                            result = result.drop(columns=['Gene_Symbol'])
                        except:
                            pass

                    # Dataset naming
                    base_name = sheetname
                    if len(st.session_state.uploaded_files_data) > 1:
                        base_name = f"{Path(filename).stem}_{sheetname}"

                    if len(treatment_groups) > 1:
                        dataset_name = f"{base_name}_{treatment}_vs_{control_group}"
                    else:
                        dataset_name = base_name

                    st.session_state.datasets[dataset_name] = result
                    processed += 1

            except Exception as e:
                st.warning(f"Failed to process {sheetname}: {str(e)}")

    st.success(f"Generated **{processed}** dataset(s)")
    st.rerun()


def process_all_intensity_manual(
    sheets_list: List[Tuple[str, str]],
    gene_col: str,
    group_a_cols: List[str], group_b_cols: List[str],
    is_log_transformed: bool = True,
    convert_ids: bool = False, species: str = 'human'
):
    """Process intensity data with manually selected groups."""
    mapper = st.session_state.gene_mapper

    processed = 0
    with st.spinner(f"Processing {len(sheets_list)} sheet(s)..."):
        for filename, sheetname in sheets_list:
            try:
                df = get_sheet_data(filename, sheetname)
                if df is None:
                    continue

                col_info = st.session_state.loader.get_column_info(df)
                actual_gene_col = gene_col if gene_col in df.columns else col_info.get('gene_col')

                if not actual_gene_col:
                    continue

                result = run_intensity_de(
                    df, actual_gene_col, group_a_cols, group_b_cols,
                    is_log_transformed=is_log_transformed
                )

                if convert_ids:
                    try:
                        result = mapper.convert_to_symbols(result, 'Gene', species)
                        result['Gene'] = result['Gene_Symbol']
                        result = result.drop(columns=['Gene_Symbol'])
                    except:
                        pass

                dataset_name = sheetname
                if len(st.session_state.uploaded_files_data) > 1:
                    dataset_name = f"{Path(filename).stem}_{sheetname}"

                st.session_state.datasets[dataset_name] = result
                processed += 1

            except Exception as e:
                st.warning(f"Failed: {sheetname}: {e}")

    st.success(f"Processed **{processed}** dataset(s)")
    st.rerun()


def process_all_count_datasets(
    sheets_list: List[Tuple[str, str]],
    gene_col: str,
    group_a: List[str], group_b: List[str],
    convert_ids: bool = False, species: str = 'human'
):
    """Run DE analysis on count data from multiple sheets."""
    mapper = st.session_state.gene_mapper

    processed = 0
    with st.spinner(f"Running DE for {len(sheets_list)} sheet(s)..."):
        for filename, sheetname in sheets_list:
            try:
                df = get_sheet_data(filename, sheetname)
                if df is None:
                    continue

                col_info = st.session_state.loader.get_column_info(df)
                actual_gene_col = gene_col if gene_col in df.columns else col_info.get('gene_col')

                if not actual_gene_col:
                    continue

                try:
                    result = run_deseq2(df, actual_gene_col, group_a, group_b)
                except ImportError:
                    result = run_simple_de(df, actual_gene_col, group_a, group_b)

                if convert_ids:
                    try:
                        result = mapper.convert_to_symbols(result, 'Gene', species)
                        result['Gene'] = result['Gene_Symbol']
                        result = result.drop(columns=['Gene_Symbol'])
                    except:
                        pass

                dataset_name = sheetname
                if len(st.session_state.uploaded_files_data) > 1:
                    dataset_name = f"{Path(filename).stem}_{sheetname}"

                st.session_state.datasets[dataset_name] = result
                processed += 1

            except Exception as e:
                st.warning(f"Failed: {sheetname}: {e}")

    st.success(f"Analyzed **{processed}** dataset(s)")
    st.rerun()


# Legacy functions for backwards compatibility
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


def process_intensity_datasets(
    sheets: List[str], gene_col: str,
    detected_groups: dict, control_group: str, treatment_groups: List[str],
    is_log_transformed: bool = True,
    convert_ids: bool = False, species: str = 'human'
):
    """Process intensity/proteomics data to calculate DE for each treatment vs control."""
    mapper = st.session_state.gene_mapper

    with st.spinner("Calculating differential expression..."):
        control_cols = detected_groups[control_group]

        for sheet in sheets:
            try:
                df = st.session_state.loader.get_sheet(sheet)

                # Check if gene column exists
                if gene_col not in df.columns:
                    st.warning(f"Sheet '{sheet}' missing gene column '{gene_col}'. Skipping.")
                    continue

                # Run DE for each treatment group
                for treatment in treatment_groups:
                    treatment_cols = detected_groups[treatment]

                    # Check if all columns exist
                    missing_cols = [c for c in control_cols + treatment_cols if c not in df.columns]
                    if missing_cols:
                        st.warning(f"Sheet '{sheet}' missing columns: {missing_cols[:3]}... Skipping.")
                        continue

                    result = run_intensity_de(
                        df,
                        gene_col,
                        control_cols,
                        treatment_cols,
                        is_log_transformed=is_log_transformed
                    )

                    # Convert IDs if requested
                    if convert_ids:
                        try:
                            result = mapper.convert_to_symbols(result, 'Gene', species)
                            result['Gene'] = result['Gene_Symbol']
                            result = result.drop(columns=['Gene_Symbol'])
                        except Exception as e:
                            st.warning(f"ID conversion failed: {str(e)}")

                    # Create dataset name
                    if len(treatment_groups) > 1:
                        dataset_name = f"{sheet}_{treatment}_vs_{control_group}"
                    else:
                        dataset_name = sheet

                    st.session_state.datasets[dataset_name] = result

            except Exception as e:
                st.warning(f"Failed to analyze {sheet}: {str(e)}")

    n_datasets = len(st.session_state.datasets)
    st.success(f"Generated {n_datasets} dataset(s) with log2FC and p-values")
    st.rerun()


def process_intensity_manual(
    sheets: List[str], gene_col: str,
    group_a_cols: List[str], group_b_cols: List[str],
    is_log_transformed: bool = True,
    convert_ids: bool = False, species: str = 'human'
):
    """Process intensity data with manually selected sample groups."""
    mapper = st.session_state.gene_mapper

    with st.spinner("Calculating differential expression..."):
        for sheet in sheets:
            try:
                df = st.session_state.loader.get_sheet(sheet)

                if gene_col not in df.columns:
                    st.warning(f"Sheet '{sheet}' missing gene column '{gene_col}'. Skipping.")
                    continue

                result = run_intensity_de(
                    df,
                    gene_col,
                    group_a_cols,
                    group_b_cols,
                    is_log_transformed=is_log_transformed
                )

                # Convert IDs if requested
                if convert_ids:
                    try:
                        result = mapper.convert_to_symbols(result, 'Gene', species)
                        result['Gene'] = result['Gene_Symbol']
                        result = result.drop(columns=['Gene_Symbol'])
                    except Exception as e:
                        st.warning(f"ID conversion failed: {str(e)}")

                st.session_state.datasets[sheet] = result

            except Exception as e:
                st.warning(f"Failed to analyze {sheet}: {str(e)}")

    st.success(f"Analyzed {len(st.session_state.datasets)} datasets")
    st.rerun()


def dashboard_tab(log2fc_threshold: float, pvalue_threshold: float, use_padj: bool):
    """Dashboard overview of all datasets and key statistics."""
    st.header("📊 Dashboard - Dataset Overview")

    if not st.session_state.datasets:
        st.info("Please load datasets in the Data Import tab first.")
        return

    datasets = st.session_state.datasets
    pval_col = 'padj' if use_padj else 'pvalue'

    # Summary metrics
    st.subheader("Summary Statistics")

    total_datasets = len(datasets)
    col_m1, col_m2, col_m3, col_m4 = st.columns(4)

    # Total genes across all datasets
    all_genes = set()
    for df in datasets.values():
        all_genes.update(df['Gene'].values)

    col_m1.metric("Total Datasets", total_datasets)
    col_m2.metric("Unique Genes", len(all_genes))

    # Genes in multiple datasets
    analyzer = DEGAnalyzer(log2fc_threshold, pvalue_threshold, use_padj)
    deg_sets = analyzer.get_deg_sets(datasets, 'both')

    genes_in_2plus = sum(1 for gene in all_genes if sum(1 for s in deg_sets.values() if gene in s) >= 2)
    genes_in_all = sum(1 for gene in all_genes if all(gene in s for s in deg_sets.values()))

    col_m3.metric("In 2+ Datasets", genes_in_2plus)
    col_m4.metric("In All Datasets", genes_in_all)

    st.markdown("---")

    # Per-dataset stats table
    st.subheader("Per-Dataset Statistics")

    stats_data = []
    for name, df in datasets.items():
        sig_mask = (df['log2FC'].abs() >= log2fc_threshold) & (df[pval_col] <= pvalue_threshold)
        up_mask = (df['log2FC'] >= log2fc_threshold) & (df[pval_col] <= pvalue_threshold)
        down_mask = (df['log2FC'] <= -log2fc_threshold) & (df[pval_col] <= pvalue_threshold)

        stats_data.append({
            'Dataset': name,
            'Total Genes': len(df),
            'Significant': sig_mask.sum(),
            'Upregulated': up_mask.sum(),
            'Downregulated': down_mask.sum(),
            '% Significant': f"{sig_mask.sum() / len(df) * 100:.1f}%"
        })

    stats_df = pd.DataFrame(stats_data)
    st.dataframe(stats_df, use_container_width=True, hide_index=True)

    # Download stats
    csv = stats_df.to_csv(index=False)
    st.download_button(
        "📥 Download Statistics (CSV)",
        csv,
        "dataset_statistics.csv",
        "text/csv"
    )

    st.markdown("---")

    # Intersection overview
    st.subheader("Intersection Overview")

    col_viz1, col_viz2 = st.columns(2)

    with col_viz1:
        # Bar chart: genes appearing in N datasets
        gene_counts = {}
        for gene in all_genes:
            count = sum(1 for s in deg_sets.values() if gene in s)
            if count > 0:
                gene_counts[count] = gene_counts.get(count, 0) + 1

        import plotly.graph_objects as go

        if st.session_state.dark_mode:
            paper_bg = '#1E1E1E'
            plot_bg = '#2D2D2D'
            font_color = '#FFFFFF'
            bar_color = '#4ECDC4'
        else:
            paper_bg = '#FFFFFF'
            plot_bg = '#FAFAFA'
            font_color = '#2C3E50'
            bar_color = '#3498DB'

        fig = go.Figure(data=[
            go.Bar(
                x=list(gene_counts.keys()),
                y=list(gene_counts.values()),
                marker_color=bar_color,
                text=list(gene_counts.values()),
                textposition='outside'
            )
        ])

        fig.update_layout(
            title="Gene Distribution Across Datasets",
            xaxis_title="Number of Datasets",
            yaxis_title="Number of Genes",
            paper_bgcolor=paper_bg,
            plot_bgcolor=plot_bg,
            font=dict(color=font_color),
            xaxis=dict(color=font_color, gridcolor='#404040' if st.session_state.dark_mode else '#E0E0E0'),
            yaxis=dict(color=font_color, gridcolor='#404040' if st.session_state.dark_mode else '#E0E0E0'),
            height=400
        )

        st.plotly_chart(fig, use_container_width=True)

    with col_viz2:
        # Pie chart: up vs down regulation
        total_up = sum(
            ((df['log2FC'] >= log2fc_threshold) & (df[pval_col] <= pvalue_threshold)).sum()
            for df in datasets.values()
        )
        total_down = sum(
            ((df['log2FC'] <= -log2fc_threshold) & (df[pval_col] <= pvalue_threshold)).sum()
            for df in datasets.values()
        )

        fig = go.Figure(data=[go.Pie(
            labels=['Upregulated', 'Downregulated'],
            values=[total_up, total_down],
            marker=dict(colors=['#FF6B6B' if st.session_state.dark_mode else '#E74C3C',
                               '#4ECDC4' if st.session_state.dark_mode else '#3498DB'])
        )])

        fig.update_layout(
            title="Overall Regulation Direction",
            paper_bgcolor=paper_bg,
            font=dict(color=font_color),
            height=400
        )

        st.plotly_chart(fig, use_container_width=True)

    st.markdown("---")

    # Quick gene lists
    st.subheader("Quick Gene Lists")

    if genes_in_all > 0:
        with st.expander(f"✅ Genes in ALL {total_datasets} datasets ({genes_in_all} genes)"):
            genes_list = sorted([gene for gene in all_genes if all(gene in s for s in deg_sets.values())])
            st.code(", ".join(genes_list), language=None)
            st.download_button(
                "📋 Download",
                "\n".join(genes_list),
                "genes_in_all_datasets.txt",
                key="dash_all"
            )

    if genes_in_2plus > 0:
        with st.expander(f"🔗 Genes in 2+ datasets ({genes_in_2plus} genes)"):
            genes_list = sorted([gene for gene in all_genes if sum(1 for s in deg_sets.values() if gene in s) >= 2])
            if len(genes_list) <= 100:
                st.code(", ".join(genes_list), language=None)
            else:
                st.caption(f"Showing first 100 of {len(genes_list)} genes:")
                st.code(", ".join(genes_list[:100]), language=None)
            st.download_button(
                "📋 Download",
                "\n".join(genes_list),
                "genes_in_2plus_datasets.txt",
                key="dash_2plus"
            )


def concordance_tab(log2fc_threshold: float, pvalue_threshold: float, use_padj: bool):
    """Show genes regulated in the same direction across datasets (concordance)."""
    st.header("🔄 Concordance Analysis")

    if len(st.session_state.datasets) < 2:
        st.info("Please load at least 2 datasets to analyze concordance.")
        return

    datasets = st.session_state.datasets
    pval_col = 'padj' if use_padj else 'pvalue'

    st.markdown("""
    **Concordance** shows whether genes are regulated in the **same direction** across datasets.
    - **Concordant**: Same direction (all up or all down)
    - **Discordant**: Mixed directions (up in some, down in others)

    This helps identify **robust biological hits** vs artifacts.
    """)

    st.markdown("---")

    # Build concordance matrix
    analyzer = DEGAnalyzer(log2fc_threshold, pvalue_threshold, use_padj)

    # Get all genes that appear in at least 2 datasets
    all_genes = set()
    for df in datasets.values():
        all_genes.update(df['Gene'].values)

    # Filter to genes in 2+ datasets
    min_datasets_filter = st.slider(
        "Minimum datasets gene must appear in",
        min_value=2,
        max_value=len(datasets),
        value=2
    )

    concordance_data = []
    for gene in all_genes:
        fc_values = []
        sig_in_datasets = []

        for name, df in datasets.items():
            gene_row = df[df['Gene'] == gene]
            if len(gene_row) > 0:
                fc = gene_row.iloc[0]['log2FC']
                pval = gene_row.iloc[0][pval_col]

                if abs(fc) >= log2fc_threshold and pval <= pvalue_threshold:
                    fc_values.append(fc)
                    sig_in_datasets.append(name)

        if len(fc_values) >= min_datasets_filter:
            # Check concordance
            all_positive = all(fc > 0 for fc in fc_values)
            all_negative = all(fc < 0 for fc in fc_values)

            if all_positive:
                direction = 'Concordant Up'
                color_code = 1
            elif all_negative:
                direction = 'Concordant Down'
                color_code = -1
            else:
                direction = 'Discordant'
                color_code = 0

            concordance_data.append({
                'Gene': gene,
                'Datasets': len(fc_values),
                'Direction': direction,
                'Mean_log2FC': np.mean(fc_values),
                'FC_Values': fc_values,
                'Datasets_List': sig_in_datasets,
                'Color': color_code
            })

    if not concordance_data:
        st.warning("No genes found in multiple datasets with current thresholds.")
        return

    concordance_df = pd.DataFrame(concordance_data)

    # Summary metrics
    st.subheader("Concordance Summary")
    col1, col2, col3, col4 = st.columns(4)

    n_concordant_up = len(concordance_df[concordance_df['Direction'] == 'Concordant Up'])
    n_concordant_down = len(concordance_df[concordance_df['Direction'] == 'Concordant Down'])
    n_discordant = len(concordance_df[concordance_df['Direction'] == 'Discordant'])
    total = len(concordance_df)

    col1.metric("Total Genes", total)
    col2.metric("✅ Concordant Up", n_concordant_up, delta=f"{n_concordant_up/total*100:.1f}%")
    col3.metric("✅ Concordant Down", n_concordant_down, delta=f"{n_concordant_down/total*100:.1f}%")
    col4.metric("⚠️ Discordant", n_discordant, delta=f"{n_discordant/total*100:.1f}%")

    st.markdown("---")

    # Filter options
    direction_filter = st.selectbox(
        "Show genes:",
        options=['All', 'Concordant Only', 'Concordant Up', 'Concordant Down', 'Discordant Only']
    )

    if direction_filter == 'Concordant Only':
        filtered_df = concordance_df[concordance_df['Direction'].str.contains('Concordant')]
    elif direction_filter == 'Concordant Up':
        filtered_df = concordance_df[concordance_df['Direction'] == 'Concordant Up']
    elif direction_filter == 'Concordant Down':
        filtered_df = concordance_df[concordance_df['Direction'] == 'Concordant Down']
    elif direction_filter == 'Discordant Only':
        filtered_df = concordance_df[concordance_df['Direction'] == 'Discordant']
    else:
        filtered_df = concordance_df

    # Sort by number of datasets
    filtered_df = filtered_df.sort_values(['Datasets', 'Mean_log2FC'], ascending=[False, False])

    st.subheader(f"{direction_filter} ({len(filtered_df)} genes)")

    # Heatmap visualization
    if len(filtered_df) > 0:
        # Build matrix: genes x datasets
        max_genes_to_show = st.slider("Max genes to show in heatmap", 10, 200, 50)
        genes_to_show = filtered_df.head(max_genes_to_show)['Gene'].tolist()

        heatmap_data = []
        for gene in genes_to_show:
            row_data = {'Gene': gene}
            for ds_name in datasets.keys():
                df = datasets[ds_name]
                gene_row = df[df['Gene'] == gene]
                if len(gene_row) > 0:
                    fc = gene_row.iloc[0]['log2FC']
                    pval = gene_row.iloc[0][pval_col]
                    if abs(fc) >= log2fc_threshold and pval <= pvalue_threshold:
                        row_data[ds_name] = fc
                    else:
                        row_data[ds_name] = 0  # Not significant
                else:
                    row_data[ds_name] = np.nan

            heatmap_data.append(row_data)

        heatmap_df = pd.DataFrame(heatmap_data).set_index('Gene')

        # Create heatmap
        import plotly.graph_objects as go

        if st.session_state.dark_mode:
            colorscale = [[0, '#4ECDC4'], [0.5, '#2D2D2D'], [1, '#FF6B6B']]
            paper_bg = '#1E1E1E'
            plot_bg = '#2D2D2D'
            font_color = '#FFFFFF'
        else:
            colorscale = [[0, '#3498DB'], [0.5, '#FFFFFF'], [1, '#E74C3C']]
            paper_bg = '#FFFFFF'
            plot_bg = '#FAFAFA'
            font_color = '#2C3E50'

        max_val = min(5, heatmap_df.abs().max().max()) if not heatmap_df.empty else 2

        fig = go.Figure(data=go.Heatmap(
            z=heatmap_df.values,
            x=heatmap_df.columns.tolist(),
            y=heatmap_df.index.tolist(),
            colorscale=colorscale,
            zmin=-max_val,
            zmax=max_val,
            colorbar=dict(title=dict(text="log2FC", font=dict(color=font_color))),
            hovertemplate='<b>%{y}</b><br>%{x}<br>log2FC: %{z:.2f}<extra></extra>'
        ))

        fig.update_layout(
            title=f"Concordance Heatmap - {direction_filter}",
            xaxis_title="Dataset",
            yaxis_title="Gene",
            paper_bgcolor=paper_bg,
            plot_bgcolor=plot_bg,
            font=dict(color=font_color),
            xaxis=dict(tickangle=45, color=font_color),
            yaxis=dict(autorange="reversed", color=font_color),
            height=max(400, len(genes_to_show) * 20 + 150)
        )

        st.plotly_chart(fig, use_container_width=True)

        # Gene lists
        st.markdown("---")
        st.subheader("Gene Lists")

        # Concordant up
        conc_up = filtered_df[filtered_df['Direction'] == 'Concordant Up']
        if len(conc_up) > 0:
            with st.expander(f"🔺 Concordant Up ({len(conc_up)} genes)"):
                genes_list = sorted(conc_up['Gene'].tolist())
                st.code(", ".join(genes_list), language=None)
                st.download_button(
                    "📋 Download",
                    "\n".join(genes_list),
                    "concordant_up.txt",
                    key="conc_up"
                )

        # Concordant down
        conc_down = filtered_df[filtered_df['Direction'] == 'Concordant Down']
        if len(conc_down) > 0:
            with st.expander(f"🔻 Concordant Down ({len(conc_down)} genes)"):
                genes_list = sorted(conc_down['Gene'].tolist())
                st.code(", ".join(genes_list), language=None)
                st.download_button(
                    "📋 Download",
                    "\n".join(genes_list),
                    "concordant_down.txt",
                    key="conc_down"
                )

        # Discordant
        disc = filtered_df[filtered_df['Direction'] == 'Discordant']
        if len(disc) > 0:
            with st.expander(f"⚠️ Discordant ({len(disc)} genes)"):
                genes_list = sorted(disc['Gene'].tolist())
                if len(genes_list) <= 50:
                    st.code(", ".join(genes_list), language=None)
                else:
                    st.caption(f"Showing first 50 of {len(genes_list)} genes:")
                    st.code(", ".join(genes_list[:50]), language=None)
                st.download_button(
                    "📋 Download",
                    "\n".join(genes_list),
                    "discordant.txt",
                    key="disc"
                )


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
                plot_bgcolor=paper_bgcolor,
                font=dict(color=font_color),
                height=max(400, len(selected_datasets) * 50),
                xaxis=dict(tickangle=45, color=font_color),
                yaxis=dict(color=font_color),
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
                    textfont=dict(color=font_color),
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
                    annotation_text=f"All {total_datasets} datasets",
                    annotation_font_color=font_color
                )

                fig.update_layout(
                    title=f"Top {len(display_df)} Recurring DEGs",
                    xaxis_title="Gene",
                    yaxis_title="Number of Datasets",
                    paper_bgcolor=paper_bgcolor,
                    plot_bgcolor=plot_bgcolor,
                    font=dict(color=font_color),
                    xaxis=dict(tickangle=45, color=font_color, gridcolor='#404040' if st.session_state.dark_mode else '#E0E0E0'),
                    yaxis=dict(color=font_color, gridcolor='#404040' if st.session_state.dark_mode else '#E0E0E0'),
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
                    xaxis=dict(tickangle=45, color=font_color),
                    yaxis=dict(autorange="reversed", color=font_color),
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
