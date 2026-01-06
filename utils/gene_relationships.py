"""
Gene relationship detection and correlation analysis.

Auto-detects co-regulated genes and patterns across datasets.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple
from scipy.stats import pearsonr, spearmanr


class GeneRelationshipAnalyzer:
    """Analyze relationships between genes across datasets."""

    def __init__(self, datasets: Dict[str, pd.DataFrame]):
        """
        Initialize with datasets.

        Args:
            datasets: Dictionary of dataset name to DataFrame
        """
        self.datasets = datasets

    def find_correlated_genes(
        self,
        gene: str,
        min_datasets: int = 3,
        min_correlation: float = 0.7,
        use_spearman: bool = False
    ) -> List[Dict]:
        """
        Find genes that correlate with the given gene across datasets.

        Args:
            gene: Gene to find correlations for
            min_datasets: Minimum datasets where both genes must be present
            min_correlation: Minimum correlation coefficient (absolute value)
            use_spearman: Use Spearman instead of Pearson correlation

        Returns:
            List of dictionaries with correlated genes and stats
        """
        gene_upper = gene.upper()
        results = []

        # Get log2FC values for the query gene across datasets
        query_values = {}
        for ds_name, df in self.datasets.items():
            gene_data = df[df['Gene'].str.upper() == gene_upper]
            if len(gene_data) > 0:
                query_values[ds_name] = gene_data.iloc[0]['log2FC']

        if len(query_values) < min_datasets:
            return []

        # Get all unique genes across datasets
        all_genes = set()
        for df in self.datasets.values():
            all_genes.update(df['Gene'].str.upper().unique())

        all_genes.discard(gene_upper)

        # Calculate correlations
        for other_gene in all_genes:
            other_values = {}
            for ds_name, df in self.datasets.items():
                if ds_name in query_values:  # Only use datasets where query gene exists
                    gene_data = df[df['Gene'].str.upper() == other_gene]
                    if len(gene_data) > 0:
                        other_values[ds_name] = gene_data.iloc[0]['log2FC']

            # Need at least min_datasets overlapping data points
            common_datasets = set(query_values.keys()) & set(other_values.keys())
            if len(common_datasets) < min_datasets:
                continue

            # Extract values in same order
            x = [query_values[ds] for ds in common_datasets]
            y = [other_values[ds] for ds in common_datasets]

            # Calculate correlation
            try:
                if use_spearman:
                    corr, pval = spearmanr(x, y)
                else:
                    corr, pval = pearsonr(x, y)

                if abs(corr) >= min_correlation and pval < 0.05:
                    results.append({
                        'gene': other_gene,
                        'correlation': corr,
                        'pvalue': pval,
                        'n_datasets': len(common_datasets),
                        'direction': 'same' if corr > 0 else 'opposite'
                    })
            except:
                continue

        # Sort by absolute correlation
        results.sort(key=lambda x: abs(x['correlation']), reverse=True)
        return results

    def find_co_significant(
        self,
        min_datasets: int = 3,
        min_genes: int = 2
    ) -> List[Dict]:
        """
        Find groups of genes that are significant together across datasets.

        Args:
            min_datasets: Minimum datasets where genes are co-significant
            min_genes: Minimum genes in a group

        Returns:
            List of gene groups
        """
        # Build dataset -> significant genes mapping
        sig_by_dataset = {}
        for ds_name, df in self.datasets.items():
            # Assuming standard significance criteria
            sig_genes = set(df[
                (df['log2FC'].abs() >= 1.0) &
                (df.get('padj', df.get('pvalue', pd.Series([1.0]))) <= 0.05)
            ]['Gene'].str.upper().tolist())
            sig_by_dataset[ds_name] = sig_genes

        # Find gene pairs that co-occur
        from collections import defaultdict
        co_occurrence = defaultdict(int)

        for genes in sig_by_dataset.values():
            genes_list = list(genes)
            for i, g1 in enumerate(genes_list):
                for g2 in genes_list[i+1:]:
                    pair = tuple(sorted([g1, g2]))
                    co_occurrence[pair] += 1

        # Filter by min_datasets
        results = []
        for (g1, g2), count in co_occurrence.items():
            if count >= min_datasets:
                results.append({
                    'genes': [g1, g2],
                    'n_datasets': count,
                    'datasets': [ds for ds, genes in sig_by_dataset.items()
                                if g1 in genes and g2 in genes]
                })

        results.sort(key=lambda x: x['n_datasets'], reverse=True)
        return results

    def get_gene_signature_similarity(
        self,
        gene: str,
        top_n: int = 10
    ) -> List[Dict]:
        """
        Find genes with similar expression patterns (signature) across datasets.

        Uses direction (+/-) rather than magnitude.

        Args:
            gene: Gene to compare
            top_n: Top N similar genes to return

        Returns:
            List of similar genes with similarity scores
        """
        gene_upper = gene.upper()

        # Get signature (direction) for query gene
        query_signature = {}
        for ds_name, df in self.datasets.items():
            gene_data = df[df['Gene'].str.upper() == gene_upper]
            if len(gene_data) > 0:
                fc = gene_data.iloc[0]['log2FC']
                query_signature[ds_name] = 1 if fc > 0 else -1 if fc < 0 else 0

        if len(query_signature) < 2:
            return []

        # Get all genes
        all_genes = set()
        for df in self.datasets.values():
            all_genes.update(df['Gene'].str.upper().unique())
        all_genes.discard(gene_upper)

        # Calculate signature similarity
        results = []
        for other_gene in all_genes:
            other_signature = {}
            for ds_name, df in self.datasets.items():
                if ds_name in query_signature:
                    gene_data = df[df['Gene'].str.upper() == other_gene]
                    if len(gene_data) > 0:
                        fc = gene_data.iloc[0]['log2FC']
                        other_signature[ds_name] = 1 if fc > 0 else -1 if fc < 0 else 0

            # Calculate agreement
            common_ds = set(query_signature.keys()) & set(other_signature.keys())
            if len(common_ds) < 2:
                continue

            agreements = sum(1 for ds in common_ds
                           if query_signature[ds] == other_signature[ds])
            similarity = agreements / len(common_ds)

            if similarity >= 0.7:  # 70% agreement
                results.append({
                    'gene': other_gene,
                    'similarity': similarity,
                    'n_datasets': len(common_ds),
                    'agreement_ratio': f"{agreements}/{len(common_ds)}"
                })

        results.sort(key=lambda x: (x['similarity'], x['n_datasets']), reverse=True)
        return results[:top_n]
