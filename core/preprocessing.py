"""
Data preprocessing module for cleaning and standardizing datasets.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple


class DataPreprocessor:
    """Handles data cleaning, standardization, and column mapping."""

    def __init__(self):
        self.gene_column_name = 'Gene'
        self.log2fc_column_name = 'log2FC'
        self.pvalue_column_name = 'pvalue'
        self.padj_column_name = 'padj'

    def standardize_dataset(
        self,
        df: pd.DataFrame,
        gene_col: str,
        log2fc_col: str,
        pvalue_col: str,
        padj_col: Optional[str] = None,
        needs_log_transform: bool = False
    ) -> pd.DataFrame:
        """
        Standardize a dataset to common column names.

        Args:
            df: Input DataFrame
            gene_col: Name of gene/protein column
            log2fc_col: Name of log2 fold change column (or FC column if needs_log_transform)
            pvalue_col: Name of p-value column
            padj_col: Name of adjusted p-value column (optional)
            needs_log_transform: If True, convert FC to log2FC

        Returns:
            Standardized DataFrame
        """
        # Create standardized DataFrame
        result = pd.DataFrame()
        result[self.gene_column_name] = df[gene_col].astype(str)

        # Handle fold change - convert to log2 if needed
        fc_values = pd.to_numeric(df[log2fc_col], errors='coerce')

        if needs_log_transform:
            # Convert FC to log2FC
            # Handle negative values (some tools report FC as negative for downregulation)
            # Handle values <= 0 which can't be log-transformed directly
            fc_values = fc_values.replace(0, np.nan)

            # Check if values look like ratios (centered around 1) or already log-scale (centered around 0)
            median_val = fc_values.abs().median()
            if median_val > 0.5 and median_val < 100:  # Likely linear FC ratios
                # Standard FC: >1 = up, <1 = down
                # log2(FC) gives: >0 = up, <0 = down
                result[self.log2fc_column_name] = np.log2(fc_values.abs()) * np.sign(fc_values)
                # For standard ratios where down is <1 (e.g., 0.5 for 2-fold down)
                # We need: log2(0.5) = -1, which is correct
                # But if FC is already signed (e.g., -2 for 2-fold down), handle that
                mask_negative_fc = fc_values < 0
                if mask_negative_fc.any():
                    # Negative values mean it's signed FC, take log2 of absolute
                    result.loc[mask_negative_fc, self.log2fc_column_name] = -np.log2(fc_values[mask_negative_fc].abs())
                mask_ratio_down = (fc_values > 0) & (fc_values < 1)
                if mask_ratio_down.any():
                    result.loc[mask_ratio_down, self.log2fc_column_name] = np.log2(fc_values[mask_ratio_down])
                mask_ratio_up = fc_values >= 1
                if mask_ratio_up.any():
                    result.loc[mask_ratio_up, self.log2fc_column_name] = np.log2(fc_values[mask_ratio_up])
            else:
                # Values might already be log-scale or unusual, use as-is
                result[self.log2fc_column_name] = fc_values
        else:
            result[self.log2fc_column_name] = fc_values

        result[self.pvalue_column_name] = pd.to_numeric(df[pvalue_col], errors='coerce')

        if padj_col and padj_col in df.columns:
            result[self.padj_column_name] = pd.to_numeric(df[padj_col], errors='coerce')
        else:
            # Use p-value as padj if not available
            result[self.padj_column_name] = result[self.pvalue_column_name]

        # Clean gene names
        result[self.gene_column_name] = result[self.gene_column_name].str.strip()

        # Remove rows with missing essential values
        result = result.dropna(subset=[self.gene_column_name, self.log2fc_column_name])

        # Remove duplicate genes (keep first occurrence)
        result = result.drop_duplicates(subset=[self.gene_column_name], keep='first')

        # Reset index
        result = result.reset_index(drop=True)

        return result

    def add_significance_column(
        self,
        df: pd.DataFrame,
        log2fc_threshold: float = 1.0,
        pvalue_threshold: float = 0.05,
        use_padj: bool = True
    ) -> pd.DataFrame:
        """
        Add significance classification column.

        Args:
            df: Standardized DataFrame
            log2fc_threshold: Absolute log2FC threshold
            pvalue_threshold: P-value threshold
            use_padj: Use adjusted p-value if True

        Returns:
            DataFrame with 'Significance' column
        """
        df = df.copy()

        pval_col = self.padj_column_name if use_padj else self.pvalue_column_name

        conditions = [
            (df[self.log2fc_column_name] >= log2fc_threshold) & (df[pval_col] <= pvalue_threshold),
            (df[self.log2fc_column_name] <= -log2fc_threshold) & (df[pval_col] <= pvalue_threshold),
        ]
        choices = ['Up', 'Down']

        df['Significance'] = np.select(conditions, choices, default='NS')

        return df

    def filter_significant(
        self,
        df: pd.DataFrame,
        log2fc_threshold: float = 1.0,
        pvalue_threshold: float = 0.05,
        use_padj: bool = True,
        direction: str = 'both'
    ) -> pd.DataFrame:
        """
        Filter to significant genes only.

        Args:
            df: Standardized DataFrame
            log2fc_threshold: Absolute log2FC threshold
            pvalue_threshold: P-value threshold
            use_padj: Use adjusted p-value if True
            direction: 'up', 'down', or 'both'

        Returns:
            Filtered DataFrame
        """
        pval_col = self.padj_column_name if use_padj else self.pvalue_column_name

        mask = df[pval_col] <= pvalue_threshold

        if direction == 'up':
            mask &= df[self.log2fc_column_name] >= log2fc_threshold
        elif direction == 'down':
            mask &= df[self.log2fc_column_name] <= -log2fc_threshold
        else:  # both
            mask &= df[self.log2fc_column_name].abs() >= log2fc_threshold

        return df[mask].copy()

    def merge_datasets(
        self,
        datasets: Dict[str, pd.DataFrame],
        on_gene: bool = True
    ) -> pd.DataFrame:
        """
        Merge multiple datasets for cross-comparison.

        Args:
            datasets: Dictionary of dataset name to DataFrame
            on_gene: Merge on gene column

        Returns:
            Merged DataFrame with suffixed columns
        """
        if not datasets:
            return pd.DataFrame()

        dataset_list = list(datasets.items())
        merged = dataset_list[0][1].copy()
        merged = merged.add_suffix(f'_{dataset_list[0][0]}')
        merged = merged.rename(columns={f'{self.gene_column_name}_{dataset_list[0][0]}': self.gene_column_name})

        for name, df in dataset_list[1:]:
            df_renamed = df.copy()
            df_renamed = df_renamed.add_suffix(f'_{name}')
            df_renamed = df_renamed.rename(columns={f'{self.gene_column_name}_{name}': self.gene_column_name})

            merged = merged.merge(df_renamed, on=self.gene_column_name, how='outer')

        return merged

    def get_common_genes(
        self,
        datasets: Dict[str, pd.DataFrame],
        significant_only: bool = False,
        log2fc_threshold: float = 1.0,
        pvalue_threshold: float = 0.05
    ) -> set:
        """
        Get genes common to all datasets.

        Args:
            datasets: Dictionary of dataset name to DataFrame
            significant_only: Only consider significant genes
            log2fc_threshold: Threshold for significance
            pvalue_threshold: Threshold for significance

        Returns:
            Set of common gene names
        """
        if not datasets:
            return set()

        gene_sets = []
        for name, df in datasets.items():
            if significant_only:
                df = self.filter_significant(df, log2fc_threshold, pvalue_threshold)
            gene_sets.append(set(df[self.gene_column_name].unique()))

        return set.intersection(*gene_sets)

    def get_overlap_matrix(
        self,
        datasets: Dict[str, pd.DataFrame],
        significant_only: bool = True,
        log2fc_threshold: float = 1.0,
        pvalue_threshold: float = 0.05
    ) -> Tuple[pd.DataFrame, Dict[str, set]]:
        """
        Create overlap matrix showing gene presence across datasets.

        Args:
            datasets: Dictionary of dataset name to DataFrame
            significant_only: Only consider significant genes
            log2fc_threshold: Threshold for significance
            pvalue_threshold: Threshold for significance

        Returns:
            Tuple of (overlap matrix DataFrame, gene sets dictionary)
        """
        gene_sets = {}

        for name, df in datasets.items():
            if significant_only:
                df = self.filter_significant(df, log2fc_threshold, pvalue_threshold)
            gene_sets[name] = set(df[self.gene_column_name].unique())

        # Get all genes
        all_genes = sorted(set.union(*gene_sets.values()) if gene_sets else set())

        # Create binary matrix
        matrix_data = {self.gene_column_name: all_genes}
        for name, genes in gene_sets.items():
            matrix_data[name] = [1 if gene in genes else 0 for gene in all_genes]

        return pd.DataFrame(matrix_data), gene_sets
