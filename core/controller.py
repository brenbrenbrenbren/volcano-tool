"""
Controller module for handling high-level data processing and analysis tasks.
This module separates business logic from the Streamlit UI and implements caching.
"""

import streamlit as st
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Optional, Dict

from core.deseq_runner import run_deseq2, run_simple_de, run_intensity_de
from core.preprocessing import DataPreprocessor
from utils.gene_mapping import GeneMapper
from core.data_loader import DataLoader

# Initialize components
# We use st.cache_resource for objects that should persist across reruns
@st.cache_resource
def get_preprocessor():
    return DataPreprocessor()

@st.cache_resource
def get_gene_mapper():
    return GeneMapper()

@st.cache_resource
def get_data_loader():
    return DataLoader()

# Cached processing functions
@st.cache_data(show_spinner=False)
def process_single_deg_sheet(
    df: pd.DataFrame,
    filename: str,
    sheetname: str,
    gene_col: str,
    fc_col: str,
    pval_col: str,
    padj_col: Optional[str],
    convert_ids: bool = False,
    species: str = 'human',
    needs_log_transform: bool = False
) -> Tuple[str, Optional[pd.DataFrame], Optional[str]]:
    """
    Process a single DEG sheet. 
    Returns (dataset_name, processed_df, error_message).
    """
    try:
        preprocessor = get_preprocessor()
        mapper = get_gene_mapper()
        
        # Convert IDs if requested
        actual_gene_col = gene_col
        if convert_ids:
            try:
                df = mapper.convert_to_symbols(df, gene_col, species)
                actual_gene_col = 'Gene_Symbol'
            except Exception as e:
                pass # Fail gracefully on ID conversion

        standardized = preprocessor.standardize_dataset(
            df, actual_gene_col, fc_col, pval_col, padj_col,
            needs_log_transform=needs_log_transform
        )
        
        # Dataset naming
        if filename:
            dataset_name = f"{Path(filename).stem}_{sheetname}"
        else:
            dataset_name = sheetname
            
        return dataset_name, standardized, None
        
    except Exception as e:
        return sheetname, None, str(e)

@st.cache_data(show_spinner=False)
def process_intensity_analysis(
    df: pd.DataFrame,
    filename: str,
    sheetname: str,
    gene_col: str,
    control_cols: List[str],
    treatment_cols: List[str],
    treatment_name: str,
    control_name: str,
    is_log_transformed: bool = True,
    convert_ids: bool = False,
    species: str = 'human'
) -> Tuple[str, Optional[pd.DataFrame], Optional[str]]:
    """
    Run intensity analysis for a single comparison.
    Returns (dataset_name, processed_df, error_message).
    """
    try:
        mapper = get_gene_mapper()
        
        result = run_intensity_de(
            df, gene_col, control_cols, treatment_cols,
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
        if filename:
            base_name = f"{Path(filename).stem}_{sheetname}"
        
        dataset_name = f"{base_name}_{treatment_name}_vs_{control_name}"
        
        return dataset_name, result, None

    except Exception as e:
        return sheetname, None, str(e)

@st.cache_data(show_spinner=False)
def process_count_analysis(
    df: pd.DataFrame,
    filename: str,
    sheetname: str,
    gene_col: str,
    group_a: List[str],
    group_b: List[str],
    convert_ids: bool = False,
    species: str = 'human'
) -> Tuple[str, Optional[pd.DataFrame], Optional[str]]:
    """
    Run count-based DE analysis (DESeq2 or simple).
    Returns (dataset_name, processed_df, error_message).
    """
    try:
        mapper = get_gene_mapper()
        
        try:
            result = run_deseq2(df, gene_col, group_a, group_b)
        except ImportError:
            result = run_simple_de(df, gene_col, group_a, group_b)

        if convert_ids:
            try:
                result = mapper.convert_to_symbols(result, 'Gene', species)
                result['Gene'] = result['Gene_Symbol']
                result = result.drop(columns=['Gene_Symbol'])
            except:
                pass

        dataset_name = sheetname
        if filename:
            dataset_name = f"{Path(filename).stem}_{sheetname}"

        return dataset_name, result, None

    except Exception as e:
        return sheetname, None, str(e)

# Cached volcano plot functions
@st.cache_data(show_spinner=False, ttl=600)
def create_volcano_cached(
    df_bytes: bytes,
    df_hash: str,
    title: str,
    log2fc_threshold: float,
    pvalue_threshold: float,
    use_padj: bool,
    dark_mode: bool,
    show_labels: bool,
    top_n_labels: int
):
    """
    Create and cache a volcano plot from DataFrame bytes.
    The df_bytes and df_hash are used for cache invalidation.

    Args:
        df_bytes: Serialized dataframe for caching
        df_hash: Hash of dataframe for cache key
        title: Plot title
        log2fc_threshold: Fold change threshold
        pvalue_threshold: P-value threshold
        use_padj: Use adjusted p-values
        dark_mode: Dark theme
        show_labels: Show gene labels
        top_n_labels: Number of top labels

    Returns:
        Plotly Figure object
    """
    import pickle
    from plotting.volcano import create_volcano_plot

    # Reconstruct dataframe from bytes
    df = pickle.loads(df_bytes)

    # Create volcano plot
    fig = create_volcano_plot(
        df,
        title=title,
        log2fc_threshold=log2fc_threshold,
        pvalue_threshold=pvalue_threshold,
        use_padj=use_padj,
        dark_mode=dark_mode,
        show_labels=show_labels,
        top_n_labels=top_n_labels
    )

    return fig

@st.cache_data(show_spinner=False, ttl=600)
def create_volcano_with_highlights_cached(
    df_bytes: bytes,
    df_hash: str,
    gene_index_hash: str,
    title: str,
    log2fc_threshold: float,
    pvalue_threshold: float,
    use_padj: bool,
    highlight_genes_tuple: Tuple,
    dark_mode: bool,
    show_labels: bool,
    top_n_labels: int,
    gene_index_json: str
):
    """
    Create and cache a volcano plot with highlighted genes.
    Optimized with gene_index for O(1) highlighting.

    Args:
        df_bytes: Serialized dataframe for caching
        df_hash: Hash of dataframe for cache key
        gene_index_hash: Hash of gene index for cache key
        title: Plot title
        log2fc_threshold: Fold change threshold
        pvalue_threshold: P-value threshold
        use_padj: Use adjusted p-values
        highlight_genes_tuple: Tuple of genes to highlight
        dark_mode: Dark theme
        show_labels: Show gene labels
        top_n_labels: Number of top labels
        gene_index_json: JSON string of gene index for reconstruction

    Returns:
        Plotly Figure object
    """
    import pickle
    import json
    from plotting.volcano import create_volcano_plot

    # Reconstruct dataframe and gene_index
    df = pickle.loads(df_bytes)
    gene_index = json.loads(gene_index_json) if gene_index_json else None

    # Create volcano plot with highlights (uses O(1) lookups if gene_index provided)
    fig = create_volcano_plot(
        df,
        title=title,
        log2fc_threshold=log2fc_threshold,
        pvalue_threshold=pvalue_threshold,
        use_padj=use_padj,
        highlight_genes=list(highlight_genes_tuple) if highlight_genes_tuple else None,
        dark_mode=dark_mode,
        show_labels=show_labels,
        top_n_labels=top_n_labels,
        gene_index=gene_index
    )

    return fig
