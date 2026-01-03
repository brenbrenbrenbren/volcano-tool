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
