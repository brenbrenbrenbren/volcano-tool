import pandas as pd
import pytest
from core.data_loader import DataLoader
import io

def test_detect_data_type_deg():
    loader = DataLoader()
    df = pd.DataFrame({
        'Gene': ['G1', 'G2'],
        'log2FC': [1.5, -0.5],
        'pvalue': [0.01, 0.05],
        'padj': [0.02, 0.1]
    })
    assert loader.detect_data_type(df) == 'deg_results'

def test_detect_data_type_counts():
    loader = DataLoader()
    df = pd.DataFrame({
        'Gene': ['G1', 'G2'],
        'Sample1': [100, 200],
        'Sample2': [150, 250],
        'Sample3': [120, 220]
    })
    assert loader.detect_data_type(df) == 'raw_counts'

def test_get_column_info():
    loader = DataLoader()
    df = pd.DataFrame({
        'GeneSymbol': ['A', 'B'],
        'Log2_Fold_Change': [1, -1],
        'P-Value': [0.05, 0.01],
        'FDR': [0.1, 0.02]
    })
    
    info = loader.get_column_info(df)
    assert info['gene_col'] == 'GeneSymbol'
    assert info['log2fc_col'] == 'Log2_Fold_Change'
    assert info['pvalue_col'] == 'P-Value'
    assert info['padj_col'] == 'FDR'
