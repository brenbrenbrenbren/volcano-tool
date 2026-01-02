"""
DESeq2-style differential expression analysis using PyDESeq2.
"""

import pandas as pd
import numpy as np
from typing import List, Optional, Tuple
import warnings


def run_deseq2(
    count_df: pd.DataFrame,
    gene_col: str,
    group_a_cols: List[str],
    group_b_cols: List[str],
    group_a_name: str = "Control",
    group_b_name: str = "Treatment"
) -> pd.DataFrame:
    """
    Run differential expression analysis using PyDESeq2.

    Args:
        count_df: DataFrame with gene counts
        gene_col: Name of gene identifier column
        group_a_cols: Column names for group A samples
        group_b_cols: Column names for group B samples
        group_a_name: Name for group A
        group_b_name: Name for group B

    Returns:
        DataFrame with DEG results (Gene, log2FC, pvalue, padj)
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        raise ImportError(
            "PyDESeq2 is required for count data analysis. "
            "Install with: pip install pydeseq2"
        )

    # Prepare count matrix
    genes = count_df[gene_col].values
    count_cols = group_a_cols + group_b_cols
    counts = count_df[count_cols].values.T  # Samples as rows

    # Create counts DataFrame (samples x genes)
    counts_df = pd.DataFrame(
        counts,
        index=count_cols,
        columns=genes
    )

    # Ensure integer counts
    counts_df = counts_df.round().astype(int)

    # Remove genes with zero counts across all samples
    counts_df = counts_df.loc[:, (counts_df.sum(axis=0) > 0)]

    # Create metadata
    metadata = pd.DataFrame({
        'sample': count_cols,
        'condition': [group_a_name] * len(group_a_cols) + [group_b_name] * len(group_b_cols)
    })
    metadata = metadata.set_index('sample')

    # Run DESeq2
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        dds = DeseqDataSet(
            counts=counts_df,
            metadata=metadata,
            design_factors="condition",
            ref_level=["condition", group_a_name]
        )

        # Run differential expression
        dds.deseq2()

        # Get results
        stat_res = DeseqStats(dds, contrast=["condition", group_b_name, group_a_name])
        stat_res.summary()

        results_df = stat_res.results_df

    # Format output
    output = pd.DataFrame({
        'Gene': results_df.index,
        'log2FC': results_df['log2FoldChange'],
        'pvalue': results_df['pvalue'],
        'padj': results_df['padj']
    })

    # Handle NaN values
    output = output.dropna(subset=['log2FC'])
    output['pvalue'] = output['pvalue'].fillna(1.0)
    output['padj'] = output['padj'].fillna(1.0)

    return output.reset_index(drop=True)


def run_simple_de(
    count_df: pd.DataFrame,
    gene_col: str,
    group_a_cols: List[str],
    group_b_cols: List[str]
) -> pd.DataFrame:
    """
    Run simple fold-change analysis when PyDESeq2 is not available.
    Uses log2 ratio of mean counts with pseudo-count.

    Args:
        count_df: DataFrame with gene counts
        gene_col: Name of gene identifier column
        group_a_cols: Column names for group A samples
        group_b_cols: Column names for group B samples

    Returns:
        DataFrame with basic DE results
    """
    from scipy import stats

    genes = count_df[gene_col].values

    # Calculate means with pseudo-count
    pseudo = 1
    mean_a = count_df[group_a_cols].mean(axis=1) + pseudo
    mean_b = count_df[group_b_cols].mean(axis=1) + pseudo

    # Log2 fold change
    log2fc = np.log2(mean_b / mean_a)

    # T-test for p-values
    pvalues = []
    for i in range(len(count_df)):
        a_vals = count_df[group_a_cols].iloc[i].values
        b_vals = count_df[group_b_cols].iloc[i].values

        # Check for variance
        if np.std(a_vals) == 0 and np.std(b_vals) == 0:
            pvalues.append(1.0)
        else:
            _, p = stats.ttest_ind(a_vals, b_vals)
            pvalues.append(p if not np.isnan(p) else 1.0)

    # Multiple testing correction (Benjamini-Hochberg)
    pvalues_array = np.array(pvalues)
    n = len(pvalues_array)
    sorted_indices = np.argsort(pvalues_array)
    sorted_pvals = pvalues_array[sorted_indices]

    # BH correction
    padj = np.zeros(n)
    cummin = 1.0
    for i in range(n - 1, -1, -1):
        corrected = sorted_pvals[i] * n / (i + 1)
        cummin = min(cummin, corrected)
        padj[sorted_indices[i]] = min(cummin, 1.0)

    output = pd.DataFrame({
        'Gene': genes,
        'log2FC': log2fc,
        'pvalue': pvalues,
        'padj': padj
    })

    return output.dropna(subset=['log2FC']).reset_index(drop=True)
