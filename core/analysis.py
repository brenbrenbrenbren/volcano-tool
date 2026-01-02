"""
DEG analysis and cross-dataset comparison functions.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Set, Tuple, Optional
from itertools import combinations


class DEGAnalyzer:
    """Analyze differentially expressed genes across datasets."""

    def __init__(
        self,
        log2fc_threshold: float = 1.0,
        pvalue_threshold: float = 0.05,
        use_padj: bool = True
    ):
        self.log2fc_threshold = log2fc_threshold
        self.pvalue_threshold = pvalue_threshold
        self.use_padj = use_padj

    def get_deg_sets(
        self,
        datasets: Dict[str, pd.DataFrame],
        direction: str = 'both'
    ) -> Dict[str, Set[str]]:
        """
        Extract DEG sets from multiple datasets.

        Args:
            datasets: Dictionary of dataset name to standardized DataFrame
            direction: 'up', 'down', or 'both'

        Returns:
            Dictionary mapping dataset names to sets of DEG names
        """
        pval_col = 'padj' if self.use_padj else 'pvalue'
        deg_sets = {}

        for name, df in datasets.items():
            mask = df[pval_col] <= self.pvalue_threshold

            if direction == 'up':
                mask &= df['log2FC'] >= self.log2fc_threshold
            elif direction == 'down':
                mask &= df['log2FC'] <= -self.log2fc_threshold
            else:
                mask &= df['log2FC'].abs() >= self.log2fc_threshold

            deg_sets[name] = set(df[mask]['Gene'].unique())

        return deg_sets

    def compute_overlaps(
        self,
        deg_sets: Dict[str, Set[str]]
    ) -> Dict[Tuple[str, ...], Set[str]]:
        """
        Compute all possible overlaps between DEG sets.

        Args:
            deg_sets: Dictionary of dataset names to gene sets

        Returns:
            Dictionary mapping set combinations to overlapping genes
        """
        dataset_names = list(deg_sets.keys())
        overlaps = {}

        # All possible combinations
        for r in range(1, len(dataset_names) + 1):
            for combo in combinations(dataset_names, r):
                # Genes in ALL sets of this combination
                intersection = set.intersection(*[deg_sets[name] for name in combo])
                # Genes NOT in other sets
                other_sets = [deg_sets[name] for name in dataset_names if name not in combo]
                if other_sets:
                    exclusive = intersection - set.union(*other_sets)
                else:
                    exclusive = intersection
                overlaps[combo] = exclusive

        return overlaps

    def get_upset_data(
        self,
        deg_sets: Dict[str, Set[str]]
    ) -> pd.DataFrame:
        """
        Prepare data for UpSet plot.

        Args:
            deg_sets: Dictionary of dataset names to gene sets

        Returns:
            DataFrame with binary membership matrix and gene names
        """
        all_genes = sorted(set.union(*deg_sets.values()) if deg_sets else set())

        data = {'Gene': all_genes}
        for name, genes in deg_sets.items():
            data[name] = [gene in genes for gene in all_genes]

        return pd.DataFrame(data)

    def get_gene_presence(
        self,
        datasets: Dict[str, pd.DataFrame],
        gene_name: str
    ) -> pd.DataFrame:
        """
        Get presence and values of a gene across all datasets.

        Args:
            datasets: Dictionary of dataset name to DataFrame
            gene_name: Gene to search for

        Returns:
            DataFrame with gene's values across datasets
        """
        results = []

        for name, df in datasets.items():
            gene_data = df[df['Gene'].str.upper() == gene_name.upper()]

            if len(gene_data) > 0:
                row = gene_data.iloc[0]
                results.append({
                    'Dataset': name,
                    'Gene': row['Gene'],
                    'log2FC': row['log2FC'],
                    'pvalue': row.get('pvalue', np.nan),
                    'padj': row.get('padj', np.nan),
                    'Significant': self._is_significant(row)
                })
            else:
                results.append({
                    'Dataset': name,
                    'Gene': gene_name,
                    'log2FC': np.nan,
                    'pvalue': np.nan,
                    'padj': np.nan,
                    'Significant': False
                })

        return pd.DataFrame(results)

    def _is_significant(self, row: pd.Series) -> bool:
        """Check if a gene row is significant."""
        pval = row.get('padj' if self.use_padj else 'pvalue', 1.0)
        log2fc = row.get('log2FC', 0)

        return (abs(log2fc) >= self.log2fc_threshold and
                pval <= self.pvalue_threshold)

    def get_concordance(
        self,
        datasets: Dict[str, pd.DataFrame],
        dataset_a: str,
        dataset_b: str
    ) -> pd.DataFrame:
        """
        Analyze concordance between two datasets.

        Args:
            datasets: Dictionary of dataset name to DataFrame
            dataset_a: First dataset name
            dataset_b: Second dataset name

        Returns:
            Merged DataFrame with concordance classification
        """
        df_a = datasets[dataset_a][['Gene', 'log2FC', 'pvalue', 'padj']].copy()
        df_b = datasets[dataset_b][['Gene', 'log2FC', 'pvalue', 'padj']].copy()

        df_a.columns = ['Gene', f'log2FC_{dataset_a}', f'pvalue_{dataset_a}', f'padj_{dataset_a}']
        df_b.columns = ['Gene', f'log2FC_{dataset_b}', f'pvalue_{dataset_b}', f'padj_{dataset_b}']

        merged = df_a.merge(df_b, on='Gene', how='inner')

        # Add significance columns
        pval_col_a = f'padj_{dataset_a}' if self.use_padj else f'pvalue_{dataset_a}'
        pval_col_b = f'padj_{dataset_b}' if self.use_padj else f'pvalue_{dataset_b}'

        merged[f'Sig_{dataset_a}'] = (
            (merged[f'log2FC_{dataset_a}'].abs() >= self.log2fc_threshold) &
            (merged[pval_col_a] <= self.pvalue_threshold)
        )
        merged[f'Sig_{dataset_b}'] = (
            (merged[f'log2FC_{dataset_b}'].abs() >= self.log2fc_threshold) &
            (merged[pval_col_b] <= self.pvalue_threshold)
        )

        # Classify concordance
        conditions = [
            merged[f'Sig_{dataset_a}'] & merged[f'Sig_{dataset_b}'] &
            (merged[f'log2FC_{dataset_a}'] > 0) & (merged[f'log2FC_{dataset_b}'] > 0),

            merged[f'Sig_{dataset_a}'] & merged[f'Sig_{dataset_b}'] &
            (merged[f'log2FC_{dataset_a}'] < 0) & (merged[f'log2FC_{dataset_b}'] < 0),

            merged[f'Sig_{dataset_a}'] & merged[f'Sig_{dataset_b}'] &
            (merged[f'log2FC_{dataset_a}'] * merged[f'log2FC_{dataset_b}'] < 0),

            merged[f'Sig_{dataset_a}'] & ~merged[f'Sig_{dataset_b}'],

            ~merged[f'Sig_{dataset_a}'] & merged[f'Sig_{dataset_b}'],
        ]

        choices = [
            'Up in both',
            'Down in both',
            'Discordant',
            f'Sig in {dataset_a} only',
            f'Sig in {dataset_b} only'
        ]

        merged['Concordance'] = np.select(conditions, choices, default='Not significant')

        return merged

    def export_deg_summary(
        self,
        datasets: Dict[str, pd.DataFrame]
    ) -> pd.DataFrame:
        """
        Create summary table of DEGs across all datasets.

        Args:
            datasets: Dictionary of dataset name to DataFrame

        Returns:
            Summary DataFrame
        """
        all_genes = set()
        for df in datasets.values():
            all_genes.update(df['Gene'].unique())

        summary_data = []

        for gene in sorted(all_genes):
            row = {'Gene': gene}
            sig_count = 0

            for name, df in datasets.items():
                gene_data = df[df['Gene'] == gene]
                if len(gene_data) > 0:
                    g = gene_data.iloc[0]
                    row[f'{name}_log2FC'] = g['log2FC']
                    row[f'{name}_pvalue'] = g.get('pvalue', np.nan)
                    row[f'{name}_padj'] = g.get('padj', np.nan)
                    row[f'{name}_Sig'] = self._is_significant(g)
                    if row[f'{name}_Sig']:
                        sig_count += 1
                else:
                    row[f'{name}_log2FC'] = np.nan
                    row[f'{name}_pvalue'] = np.nan
                    row[f'{name}_padj'] = np.nan
                    row[f'{name}_Sig'] = False

            row['Significant_in_N_datasets'] = sig_count
            summary_data.append(row)

        return pd.DataFrame(summary_data)
