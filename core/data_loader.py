"""
Data loading module for Excel and CSV files with multiple sheets/datasets.
"""

import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class DataLoader:
    """Handles loading of Excel and CSV files with multiple datasets."""

    def __init__(self):
        self.file_path: Optional[Path] = None
        self.sheet_names: List[str] = []
        self.datasets: Dict[str, pd.DataFrame] = {}
        self._excel_file = None

    def load_file(self, file_path: str) -> List[str]:
        """
        Load a file and return available sheet names.

        Args:
            file_path: Path to Excel or CSV file

        Returns:
            List of sheet/dataset names
        """
        self.file_path = Path(file_path)
        suffix = self.file_path.suffix.lower()

        if suffix in ['.xlsx', '.xls']:
            return self._load_excel()
        elif suffix == '.csv':
            return self._load_csv()
        else:
            raise ValueError(f"Unsupported file format: {suffix}")

    def load_file_from_buffer(self, uploaded_file) -> List[str]:
        """
        Load from Streamlit UploadedFile buffer.

        Args:
            uploaded_file: Streamlit UploadedFile object

        Returns:
            List of sheet/dataset names
        """
        file_name = uploaded_file.name
        suffix = Path(file_name).suffix.lower()

        if suffix in ['.xlsx', '.xls']:
            excel_file = pd.ExcelFile(uploaded_file)
            self.sheet_names = excel_file.sheet_names
            self._excel_file = excel_file
            return self.sheet_names
        elif suffix == '.csv':
            df = pd.read_csv(uploaded_file)
            dataset_name = Path(file_name).stem
            self.sheet_names = [dataset_name]
            self.datasets[dataset_name] = df
            return self.sheet_names
        else:
            raise ValueError(f"Unsupported file format: {suffix}")

    def _load_excel(self) -> List[str]:
        """Load Excel file and get sheet names."""
        excel_file = pd.ExcelFile(self.file_path)
        self.sheet_names = excel_file.sheet_names
        self._excel_file = excel_file
        return self.sheet_names

    def _load_csv(self) -> List[str]:
        """Load CSV file as single dataset."""
        df = pd.read_csv(self.file_path)
        dataset_name = self.file_path.stem
        self.sheet_names = [dataset_name]
        self.datasets[dataset_name] = df
        return self.sheet_names

    def get_sheet(self, sheet_name: str) -> pd.DataFrame:
        """
        Get a specific sheet/dataset by name.

        Args:
            sheet_name: Name of the sheet to retrieve

        Returns:
            DataFrame containing the sheet data
        """
        if sheet_name in self.datasets:
            return self.datasets[sheet_name]

        if self._excel_file is not None:
            df = pd.read_excel(self._excel_file, sheet_name=sheet_name)
            self.datasets[sheet_name] = df
            return df

        raise ValueError(f"Sheet '{sheet_name}' not found")

    def get_selected_sheets(self, sheet_names: List[str]) -> Dict[str, pd.DataFrame]:
        """
        Get multiple sheets by name.

        Args:
            sheet_names: List of sheet names to retrieve

        Returns:
            Dictionary mapping sheet names to DataFrames
        """
        return {name: self.get_sheet(name) for name in sheet_names}

    def filter_sheet_names(self, search_term: str) -> List[str]:
        """
        Filter sheet names by search term.

        Args:
            search_term: Term to search for in sheet names

        Returns:
            Filtered list of sheet names
        """
        if not search_term:
            return self.sheet_names

        search_lower = search_term.lower()
        return [name for name in self.sheet_names if search_lower in name.lower()]

    def detect_data_type(self, df: pd.DataFrame) -> str:
        """
        Detect if DataFrame contains DEG results or raw counts.

        Args:
            df: DataFrame to analyze

        Returns:
            'deg_results' or 'raw_counts'
        """
        columns_lower = [col.lower().replace(' ', '').replace('_', '').replace('-', '').replace('.', '')
                        for col in df.columns]

        # Check for DEG result indicators (comprehensive list)
        deg_indicators = [
            # Fold change indicators
            'log2foldchange', 'log2fc', 'logfc', 'lfc', 'foldchange', 'fc', 'ratio',
            # P-value indicators
            'pvalue', 'pval', 'p', 'rawp',
            # Adjusted p-value indicators
            'padj', 'fdr', 'qvalue', 'qval', 'adjp', 'adjustedp', 'bh'
        ]

        has_deg_cols = sum(1 for ind in deg_indicators
                          if any(ind in col or col == ind for col in columns_lower))

        if has_deg_cols >= 2:
            return 'deg_results'

        # Check for count data (mostly integers)
        numeric_cols = df.select_dtypes(include=['number']).columns
        if len(numeric_cols) > 0:
            sample_data = df[numeric_cols].iloc[:100]
            is_integer_like = (sample_data == sample_data.round()).all().all()
            if is_integer_like and sample_data.min().min() >= 0:
                return 'raw_counts'

        return 'deg_results'  # Default assumption

    def get_column_info(self, df: pd.DataFrame) -> Dict:
        """
        Analyze DataFrame columns and suggest mappings.
        Uses comprehensive pattern matching to handle various naming conventions.

        Args:
            df: DataFrame to analyze

        Returns:
            Dictionary with column type suggestions
        """
        columns = df.columns.tolist()
        columns_lower = [col.lower().strip() for col in columns]
        # Also create versions without special characters for matching
        columns_normalized = [
            col.lower().replace(' ', '').replace('_', '').replace('-', '').replace('.', '')
            for col in columns
        ]

        info = {
            'gene_col': None,
            'log2fc_col': None,
            'fc_col': None,  # Track non-log FC separately
            'pvalue_col': None,
            'padj_col': None,
            'count_cols': [],
            'all_columns': columns,
            'needs_log_transform': False
        }

        # =====================================================================
        # Gene/protein identifier patterns (order matters - more specific first)
        # =====================================================================
        gene_exact_matches = [
            'gene', 'gene_symbol', 'genesymbol', 'gene symbol', 'gene.symbol',
            'symbol', 'gene_name', 'genename', 'gene name', 'gene.name',
            'protein', 'protein_name', 'proteinname', 'protein name',
            'uniprot', 'uniprot_id', 'uniprotid', 'accession',
            'ensembl', 'ensembl_id', 'ensemblid', 'ensg',
            'entrez', 'entrez_id', 'entrezid', 'gene_id', 'geneid',
            'id', 'name', 'identifier', 'target'
        ]

        # First try exact matches (normalized)
        for pattern in gene_exact_matches:
            pattern_norm = pattern.replace(' ', '').replace('_', '').replace('-', '').replace('.', '')
            for i, col_norm in enumerate(columns_normalized):
                if col_norm == pattern_norm:
                    info['gene_col'] = columns[i]
                    break
            if info['gene_col']:
                break

        # Then try partial matches if no exact match
        if not info['gene_col']:
            gene_partial_patterns = ['gene', 'symbol', 'protein', 'uniprot', 'ensembl', 'entrez', 'accession']
            for pattern in gene_partial_patterns:
                for i, col in enumerate(columns_lower):
                    if pattern in col:
                        info['gene_col'] = columns[i]
                        break
                if info['gene_col']:
                    break

        # Last resort: first non-numeric column
        if not info['gene_col']:
            for i, col in enumerate(columns):
                if df[col].dtype == 'object' or str(df[col].dtype) == 'string':
                    info['gene_col'] = col
                    break

        # =====================================================================
        # Log2 fold change patterns (order matters - log2 first, then plain FC)
        # =====================================================================
        log2fc_exact_matches = [
            'log2foldchange', 'log2_foldchange', 'log2_fold_change', 'log2 fold change',
            'log2fc', 'log2_fc', 'log2 fc', 'log2(fc)', 'log2(foldchange)',
            'logfc', 'log_fc', 'log fc', 'lfc',
            'log2ratio', 'log2_ratio', 'log2 ratio',
            'logratio', 'log_ratio', 'log ratio'
        ]

        # First look for log2 FC columns
        for pattern in log2fc_exact_matches:
            pattern_norm = pattern.replace(' ', '').replace('_', '').replace('-', '').replace('.', '').replace('(', '').replace(')', '')
            for i, col_norm in enumerate(columns_normalized):
                if pattern_norm in col_norm or col_norm == pattern_norm:
                    info['log2fc_col'] = columns[i]
                    break
            if info['log2fc_col']:
                break

        # If no log2FC found, look for plain fold change
        if not info['log2fc_col']:
            fc_patterns = [
                'foldchange', 'fold_change', 'fold change', 'fold.change',
                'fc', 'ratio', 'expression_ratio'
            ]
            for pattern in fc_patterns:
                pattern_norm = pattern.replace(' ', '').replace('_', '').replace('-', '').replace('.', '')
                for i, col_norm in enumerate(columns_normalized):
                    # Avoid matching 'logfc' patterns
                    if 'log' in col_norm:
                        continue
                    if col_norm == pattern_norm or (pattern_norm in col_norm and len(pattern_norm) > 1):
                        info['fc_col'] = columns[i]
                        info['needs_log_transform'] = True
                        break
                if info['fc_col']:
                    break

        # =====================================================================
        # P-value patterns (exclude adjusted p-values)
        # =====================================================================
        pval_exact_matches = [
            'pvalue', 'p_value', 'p-value', 'p.value', 'p value',
            'pval', 'p_val', 'p-val', 'p.val',
            'rawp', 'raw_p', 'raw.p', 'raw p',
            'raw_pvalue', 'raw_p_value', 'rawpvalue'
        ]

        for pattern in pval_exact_matches:
            pattern_norm = pattern.replace(' ', '').replace('_', '').replace('-', '').replace('.', '')
            for i, col_norm in enumerate(columns_normalized):
                # Skip if it looks like adjusted p-value
                if any(adj in col_norm for adj in ['adj', 'fdr', 'bonf', 'bh', 'qval', 'qvalue']):
                    continue
                if col_norm == pattern_norm or pattern_norm == col_norm:
                    info['pvalue_col'] = columns[i]
                    break
            if info['pvalue_col']:
                break

        # Partial match for p-value
        if not info['pvalue_col']:
            for i, col in enumerate(columns_lower):
                col_norm = columns_normalized[i]
                if any(adj in col_norm for adj in ['adj', 'fdr', 'bonf', 'bh', 'qval', 'qvalue']):
                    continue
                if 'pval' in col_norm or col_norm == 'p':
                    info['pvalue_col'] = columns[i]
                    break

        # =====================================================================
        # Adjusted p-value / FDR patterns
        # =====================================================================
        padj_exact_matches = [
            'padj', 'p_adj', 'p.adj', 'p-adj', 'p adj',
            'adjp', 'adj_p', 'adj.p', 'adj-p', 'adj p',
            'adjusted_pvalue', 'adjusted_p_value', 'adjustedpvalue', 'adjusted pvalue',
            'adj_pvalue', 'adj_p_value', 'adjpvalue', 'adj pvalue',
            'fdr', 'false_discovery_rate', 'falsediscoveryrate',
            'qvalue', 'q_value', 'q-value', 'q.value', 'q value',
            'qval', 'q_val', 'q-val', 'q.val',
            'bh', 'bh_pvalue', 'bhpvalue', 'bh_adjusted',
            'bonferroni', 'bonf', 'bonf_p'
        ]

        for pattern in padj_exact_matches:
            pattern_norm = pattern.replace(' ', '').replace('_', '').replace('-', '').replace('.', '')
            for i, col_norm in enumerate(columns_normalized):
                if col_norm == pattern_norm or pattern_norm in col_norm:
                    info['padj_col'] = columns[i]
                    break
            if info['padj_col']:
                break

        # =====================================================================
        # Count columns (numeric columns that aren't FC or p-values)
        # =====================================================================
        identified_cols = {info['gene_col'], info['log2fc_col'], info['fc_col'],
                          info['pvalue_col'], info['padj_col']}
        identified_cols.discard(None)
        numeric_cols = df.select_dtypes(include=['number']).columns
        info['count_cols'] = [col for col in numeric_cols if col not in identified_cols]

        return info
