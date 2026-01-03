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
    padj = _benjamini_hochberg(pvalues)

    output = pd.DataFrame({
        'Gene': genes,
        'log2FC': log2fc,
        'pvalue': pvalues,
        'padj': padj
    })

    return output.dropna(subset=['log2FC']).reset_index(drop=True)


def run_intensity_de(
    intensity_df: pd.DataFrame,
    gene_col: str,
    group_a_cols: List[str],
    group_b_cols: List[str],
    is_log_transformed: bool = True
) -> pd.DataFrame:
    """
    Run differential expression analysis on intensity/abundance data.

    This is for proteomics or other quantitative data where values are
    intensities (or log2-transformed intensities), not raw counts.

    Args:
        intensity_df: DataFrame with intensity values
        gene_col: Name of gene/protein identifier column
        group_a_cols: Column names for group A (control) samples
        group_b_cols: Column names for group B (treatment) samples
        is_log_transformed: If True, data is already log2-transformed

    Returns:
        DataFrame with DE results (Gene, log2FC, pvalue, padj)
    """
    from scipy import stats

    genes = intensity_df[gene_col].values

    # Get intensity values
    group_a_data = intensity_df[group_a_cols].values
    group_b_data = intensity_df[group_b_cols].values

    # Calculate log2 fold change
    if is_log_transformed:
        # For log-transformed data: log2FC = mean(B) - mean(A)
        mean_a = np.nanmean(group_a_data, axis=1)
        mean_b = np.nanmean(group_b_data, axis=1)
        log2fc = mean_b - mean_a
    else:
        # For linear data: log2FC = log2(mean(B) / mean(A))
        mean_a = np.nanmean(group_a_data, axis=1)
        mean_b = np.nanmean(group_b_data, axis=1)
        # Add small value to avoid division by zero
        log2fc = np.log2((mean_b + 1e-10) / (mean_a + 1e-10))

    # T-test for p-values
    pvalues = []
    for i in range(len(intensity_df)):
        a_vals = group_a_data[i]
        b_vals = group_b_data[i]

        # Remove NaN values
        a_vals = a_vals[~np.isnan(a_vals)]
        b_vals = b_vals[~np.isnan(b_vals)]

        # Need at least 2 values in each group for t-test
        if len(a_vals) < 2 or len(b_vals) < 2:
            pvalues.append(1.0)
            continue

        # Check for zero variance
        if np.std(a_vals) == 0 and np.std(b_vals) == 0:
            pvalues.append(1.0)
        else:
            try:
                _, p = stats.ttest_ind(a_vals, b_vals, equal_var=False)  # Welch's t-test
                pvalues.append(p if not np.isnan(p) else 1.0)
            except:
                pvalues.append(1.0)

    # Multiple testing correction (Benjamini-Hochberg)
    padj = _benjamini_hochberg(pvalues)

    output = pd.DataFrame({
        'Gene': genes,
        'log2FC': log2fc,
        'pvalue': pvalues,
        'padj': padj
    })

    # Remove rows with NaN log2FC
    output = output.dropna(subset=['log2FC'])

    return output.reset_index(drop=True)


def _benjamini_hochberg(pvalues: List[float]) -> np.ndarray:
    """Apply Benjamini-Hochberg FDR correction."""
    pvalues_array = np.array(pvalues)
    n = len(pvalues_array)

    if n == 0:
        return np.array([])

    sorted_indices = np.argsort(pvalues_array)
    sorted_pvals = pvalues_array[sorted_indices]

    padj = np.zeros(n)
    cummin = 1.0
    for i in range(n - 1, -1, -1):
        corrected = sorted_pvals[i] * n / (i + 1)
        cummin = min(cummin, corrected)
        padj[sorted_indices[i]] = min(cummin, 1.0)

    return padj


def detect_sample_groups(columns: List[str]) -> dict:
    """
    Auto-detect sample groups from column names.

    Looks for patterns like:
    - ctrl, ctrl.1, ctrl.2 or control, control_1
    - KD1, KD1.1, KD1.2 or treatment, treatment_1
    - WT, WT_1, WT_2 or wt_rep1, wt_rep2

    Args:
        columns: List of column names

    Returns:
        Dictionary mapping group names to list of columns
    """
    import re

    groups = {}

    for col in columns:
        # Normalize the column name
        col_clean = col.strip()

        # Try to extract base group name
        # Only remove patterns that have a separator (., _, -, or space) followed by a number
        # This preserves "KD1" as "KD1" but converts "KD1.1" to "KD1"
        # Pattern: remove [separator][optional spaces][number(s)] at the end
        base_name = re.sub(r'[\s._-]+\d+$', '', col_clean)
        base_name = re.sub(r'\s+$', '', base_name)  # Remove trailing spaces

        if not base_name:
            base_name = col_clean

        # Normalize common control names
        base_lower = base_name.lower()
        if base_lower in ['ctrl', 'control', 'con', 'ctl', 'wt', 'wildtype', 'wild-type', 'untreated', 'vehicle']:
            base_name = 'Control'

        if base_name not in groups:
            groups[base_name] = []
        groups[base_name].append(col)

    return groups


def auto_run_de(
    df: pd.DataFrame,
    gene_col: str,
    control_group: str = None,
    treatment_groups: List[str] = None
) -> List[Tuple[str, pd.DataFrame]]:
    """
    Automatically run DE analysis on a dataframe with sample columns.

    Args:
        df: DataFrame with gene column and sample intensity columns
        gene_col: Name of gene identifier column
        control_group: Name of control group (auto-detected if None)
        treatment_groups: List of treatment group names (auto-detected if None)

    Returns:
        List of (comparison_name, results_df) tuples
    """
    # Get numeric columns (excluding gene column)
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

    # Detect groups
    groups = detect_sample_groups(numeric_cols)

    if len(groups) < 2:
        raise ValueError(f"Need at least 2 groups for comparison. Found: {list(groups.keys())}")

    # Auto-detect control group
    if control_group is None:
        control_candidates = ['Control', 'ctrl', 'Ctrl', 'CTRL', 'WT', 'wt', 'Wt']
        for candidate in control_candidates:
            if candidate in groups:
                control_group = candidate
                break

        if control_group is None:
            # Use first group as control
            control_group = list(groups.keys())[0]

    # Get treatment groups
    if treatment_groups is None:
        treatment_groups = [g for g in groups.keys() if g != control_group]

    # Detect if data is log-transformed (values typically between -10 and 10)
    sample_values = df[numeric_cols].values.flatten()
    sample_values = sample_values[~np.isnan(sample_values)]
    is_log_transformed = np.percentile(np.abs(sample_values), 95) < 20

    # Run DE for each treatment vs control
    results = []
    control_cols = groups[control_group]

    for treatment in treatment_groups:
        treatment_cols = groups[treatment]

        de_result = run_intensity_de(
            df,
            gene_col,
            control_cols,
            treatment_cols,
            is_log_transformed=is_log_transformed
        )

        comparison_name = f"{treatment}_vs_{control_group}"
        results.append((comparison_name, de_result))

    return results
