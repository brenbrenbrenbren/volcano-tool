import pandas as pd
import numpy as np
import pytest
from core.deseq_runner import detect_sample_groups, run_simple_de

def test_detect_sample_groups():
    cols = ['Ctrl_1', 'Ctrl_2', 'Treat_1', 'Treat_2', 'Treat_3']
    groups = detect_sample_groups(cols)
    
    assert 'Control' in groups
    assert 'Treat' in groups
    assert len(groups['Control']) == 2
    assert len(groups['Treat']) == 3
    # Check that it normalizes "Ctrl" to "Control"
    assert 'Ctrl' not in groups

def test_run_simple_de():
    # create dummy count data
    df = pd.DataFrame({
        'Gene': ['G1', 'G2'],
        'Ctrl1': [10, 10],
        'Ctrl2': [12, 12],
        'Treat1': [100, 100],
        'Treat2': [110, 110]
    })
    
    # Simple DE should detect upregulation
    # log2(105/11) approx 3.25
    
    res = run_simple_de(df, 'Gene', ['Ctrl1', 'Ctrl2'], ['Treat1', 'Treat2'])
    
    assert len(res) == 2
    assert 'log2FC' in res.columns
    assert 'pvalue' in res.columns
    assert res.iloc[0]['log2FC'] > 2  # G1 should be upregulated
