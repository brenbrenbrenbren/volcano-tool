"""
Gene ID mapping and conversion utilities.
"""

import pandas as pd
import re
from typing import Dict, List, Optional


class GeneMapper:
    """Map various gene/protein identifiers to gene symbols."""

    def __init__(self):
        self._cache: Dict[str, str] = {}
        self._mygene_client = None

    def _get_mygene_client(self):
        """Lazy load mygene client."""
        if self._mygene_client is None:
            try:
                import mygene
                self._mygene_client = mygene.MyGeneInfo()
            except ImportError:
                raise ImportError(
                    "mygene is required for ID conversion. "
                    "Install with: pip install mygene"
                )
        return self._mygene_client

    def detect_id_type(self, ids: List[str]) -> str:
        """
        Detect the type of identifier.

        Args:
            ids: List of identifiers

        Returns:
            Identifier type: 'symbol', 'uniprot', 'ensembl', 'entrez', or 'unknown'
        """
        sample = [str(i).strip() for i in ids[:100] if pd.notna(i)]

        if not sample:
            return 'unknown'

        # UniProt pattern (e.g., P04637, Q9Y6K9)
        uniprot_pattern = re.compile(
            r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$'
        )
        uniprot_matches = sum(1 for i in sample if uniprot_pattern.match(i))

        # Ensembl pattern (e.g., ENSG00000141510)
        ensembl_pattern = re.compile(r'^ENS[A-Z]*G[0-9]{11}$')
        ensembl_matches = sum(1 for i in sample if ensembl_pattern.match(i))

        # Entrez (numeric)
        entrez_matches = sum(1 for i in sample if i.isdigit())

        # Gene symbol (typically 2-10 chars, alphanumeric, may have hyphens)
        symbol_pattern = re.compile(r'^[A-Z][A-Z0-9\-]{1,15}$', re.IGNORECASE)
        symbol_matches = sum(1 for i in sample if symbol_pattern.match(i))

        total = len(sample)
        if total == 0:
            return 'unknown'

        scores = {
            'uniprot': uniprot_matches / total,
            'ensembl': ensembl_matches / total,
            'entrez': entrez_matches / total,
            'symbol': symbol_matches / total
        }

        best = max(scores, key=scores.get)
        if scores[best] > 0.5:
            return best
        return 'symbol'  # Default to symbol

    def convert_to_symbols(
        self,
        df: pd.DataFrame,
        id_column: str,
        species: str = 'human'
    ) -> pd.DataFrame:
        """
        Convert identifiers to gene symbols.

        Args:
            df: DataFrame with gene/protein identifiers
            id_column: Name of the identifier column
            species: Species ('human' or 'mouse')

        Returns:
            DataFrame with converted symbols
        """
        ids = df[id_column].astype(str).tolist()
        id_type = self.detect_id_type(ids)

        if id_type == 'symbol':
            # Already symbols, just clean them
            df = df.copy()
            df[id_column] = df[id_column].str.upper().str.strip()
            return df

        # Convert using mygene
        mg = self._get_mygene_client()

        # Filter out already cached
        ids_to_query = [i for i in set(ids) if i not in self._cache]

        if ids_to_query:
            scope_map = {
                'uniprot': 'uniprot',
                'ensembl': 'ensembl.gene',
                'entrez': 'entrezgene'
            }
            scope = scope_map.get(id_type, 'symbol,alias')

            species_map = {'human': 'human', 'mouse': 'mouse'}
            species_name = species_map.get(species, 'human')

            results = mg.querymany(
                ids_to_query,
                scopes=scope,
                fields='symbol',
                species=species_name,
                returnall=True
            )

            for item in results.get('out', []):
                query = item.get('query', '')
                symbol = item.get('symbol', query)
                self._cache[query] = symbol

        # Apply mapping
        df = df.copy()
        df['Gene_Symbol'] = df[id_column].map(
            lambda x: self._cache.get(str(x), str(x))
        )

        return df

    def standardize_symbols(self, df: pd.DataFrame, gene_column: str) -> pd.DataFrame:
        """
        Standardize gene symbols (uppercase, strip whitespace).

        Args:
            df: DataFrame with gene column
            gene_column: Name of gene column

        Returns:
            DataFrame with standardized symbols
        """
        df = df.copy()
        df[gene_column] = (
            df[gene_column]
            .astype(str)
            .str.strip()
            .str.upper()
        )
        return df

    def extract_symbol_from_name(self, name: str) -> str:
        """
        Extract gene symbol from a full gene/protein name.

        Args:
            name: Full gene or protein name

        Returns:
            Extracted symbol or original name
        """
        # Common patterns
        # "TP53 tumor protein p53" -> "TP53"
        # "BRCA1 DNA repair associated" -> "BRCA1"

        name = str(name).strip()

        # If starts with a standard symbol pattern, extract it
        match = re.match(r'^([A-Z][A-Z0-9\-]{1,15})\s+', name, re.IGNORECASE)
        if match:
            return match.group(1).upper()

        # Return first word if short enough
        first_word = name.split()[0] if name.split() else name
        if len(first_word) <= 15:
            return first_word.upper()

        return name
