"""
Pathway database access and gene set management.
"""

import pandas as pd
import json
import requests
from pathlib import Path
from typing import Dict, List, Optional, Set
import os


class PathwayDatabase:
    """Access pathway gene sets from MSigDB, KEGG, and Reactome."""

    def __init__(self, cache_dir: Optional[str] = None):
        if cache_dir:
            self.cache_dir = Path(cache_dir)
        else:
            self.cache_dir = Path.home() / '.volcano_tool' / 'pathways'
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        self._pathways: Dict[str, Dict[str, Set[str]]] = {}
        self._load_builtin_pathways()

    def _load_builtin_pathways(self):
        """Load commonly used pathway definitions."""
        # Built-in curated pathways
        self._pathways['curated'] = {
            'WNT_SIGNALING': {
                'CTNNB1', 'WNT1', 'WNT2', 'WNT3', 'WNT3A', 'WNT4', 'WNT5A',
                'WNT5B', 'WNT6', 'WNT7A', 'WNT7B', 'WNT8A', 'WNT8B', 'WNT9A',
                'WNT9B', 'WNT10A', 'WNT10B', 'WNT11', 'WNT16', 'FZD1', 'FZD2',
                'FZD3', 'FZD4', 'FZD5', 'FZD6', 'FZD7', 'FZD8', 'FZD9', 'FZD10',
                'LRP5', 'LRP6', 'DVL1', 'DVL2', 'DVL3', 'AXIN1', 'AXIN2', 'APC',
                'GSK3B', 'CSNK1A1', 'TCF7', 'TCF7L1', 'TCF7L2', 'LEF1', 'MYC',
                'CCND1', 'DKK1', 'DKK2', 'DKK3', 'DKK4', 'SFRP1',
                'SFRP2', 'SFRP4', 'SFRP5', 'WIF1', 'RSPO1', 'RSPO2', 'RSPO3', 'RSPO4'
            },
            'IL17_SIGNALING': {
                'IL17A', 'IL17B', 'IL17C', 'IL17D', 'IL17E', 'IL17F', 'IL17RA',
                'IL17RB', 'IL17RC', 'IL17RD', 'IL17RE', 'TRAF6', 'TRAF3',
                'TAB1', 'TAB2', 'TAB3', 'MAP3K7', 'NFKB1', 'NFKB2', 'RELA',
                'RELB', 'REL', 'CXCL1', 'CXCL2', 'CXCL5', 'CXCL8', 'CCL2',
                'CCL7', 'CCL20', 'IL6', 'CSF2', 'CSF3', 'MMP1', 'MMP3', 'MMP9',
                'MMP13', 'DEFB4A', 'S100A7', 'S100A8', 'S100A9', 'LCN2'
            },
            'GLYCOLYSIS': {
                'HK1', 'HK2', 'HK3', 'GCK', 'GPI', 'PFKM', 'PFKL', 'PFKP',
                'ALDOA', 'ALDOB', 'ALDOC', 'TPI1', 'GAPDH', 'PGK1', 'PGK2',
                'PGAM1', 'PGAM2', 'ENO1', 'ENO2', 'ENO3', 'PKM', 'PKLR',
                'LDHA', 'LDHB', 'LDHC', 'PDK1', 'PDK2', 'PDK3', 'PDK4'
            },
            'TCA_CYCLE': {
                'CS', 'ACO1', 'ACO2', 'IDH1', 'IDH2', 'IDH3A', 'IDH3B', 'IDH3G',
                'OGDH', 'DLST', 'DLD', 'SUCLG1', 'SUCLG2', 'SUCLA2', 'SDHA',
                'SDHB', 'SDHC', 'SDHD', 'FH', 'MDH1', 'MDH2', 'PC', 'PCK1', 'PCK2'
            },
            'OXIDATIVE_PHOSPHORYLATION': {
                'NDUFA1', 'NDUFA2', 'NDUFA3', 'NDUFA4', 'NDUFA5', 'NDUFA6',
                'NDUFB1', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6',
                'NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS5', 'NDUFS6',
                'SDHA', 'SDHB', 'SDHC', 'SDHD', 'UQCRB', 'UQCRC1', 'UQCRC2',
                'CYC1', 'CYCS', 'COX4I1', 'COX5A', 'COX5B', 'COX6A1', 'COX6B1',
                'COX7A1', 'COX7B', 'COX8A', 'ATP5F1A', 'ATP5F1B', 'ATP5F1C',
                'ATP5PB', 'ATP5PD', 'ATP5PF', 'ATP5PO'
            },
            'P_BODIES': {
                'DCP1A', 'DCP1B', 'DCP2', 'EDC3', 'EDC4', 'LSM1', 'LSM2',
                'LSM3', 'LSM4', 'LSM5', 'LSM6', 'LSM7', 'LSM14A', 'LSM14B',
                'DDX6', 'XRN1', 'PATL1', 'EIF4E', 'EIF4ENIF1', 'CPEB1',
                'TNRC6A', 'TNRC6B', 'TNRC6C', 'AGO1', 'AGO2', 'AGO3', 'AGO4',
                'MOV10', 'YTHDF2', 'ZFP36', 'ZFP36L1', 'ZFP36L2'
            },
            'APOPTOSIS': {
                'BCL2', 'BCL2L1', 'BCL2L11', 'MCL1', 'BAX', 'BAK1', 'BID',
                'BAD', 'PMAIP1', 'BBC3', 'CASP3', 'CASP7', 'CASP8', 'CASP9',
                'APAF1', 'CYCS', 'DIABLO', 'XIAP', 'BIRC2', 'BIRC3', 'BIRC5',
                'TP53', 'FAS', 'FASLG', 'TNFRSF10A', 'TNFRSF10B', 'TRADD'
            },
            'CELL_CYCLE': {
                'CDK1', 'CDK2', 'CDK4', 'CDK6', 'CCNA1', 'CCNA2', 'CCNB1',
                'CCNB2', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CCNE2',
                'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'RB1', 'E2F1', 'E2F2',
                'E2F3', 'TP53', 'MDM2', 'ATM', 'ATR', 'CHEK1', 'CHEK2'
            },
            'INFLAMMATORY_RESPONSE': {
                'TNF', 'IL1A', 'IL1B', 'IL6', 'IL8', 'CXCL8', 'CCL2', 'CCL3',
                'CCL4', 'CCL5', 'CXCL1', 'CXCL2', 'CXCL10', 'ICAM1', 'VCAM1',
                'SELE', 'SELP', 'PTGS2', 'NOS2', 'NFKB1', 'RELA', 'IKBKB',
                'IKBKG', 'NFKBIA', 'MYD88', 'TLR2', 'TLR4'
            },
            'EMT_MARKERS': {
                'CDH1', 'CDH2', 'VIM', 'FN1', 'SNAI1', 'SNAI2', 'TWIST1',
                'TWIST2', 'ZEB1', 'ZEB2', 'CTNNB1', 'TGFB1', 'TGFB2', 'TGFB3',
                'SMAD2', 'SMAD3', 'SMAD4', 'MMP2', 'MMP9', 'OCLN', 'TJP1'
            }
        }

    def search_pathways(self, query: str) -> List[Dict]:
        """
        Search for pathways matching a query.

        Args:
            query: Search term

        Returns:
            List of matching pathway info dicts
        """
        results = []
        query_lower = query.lower()

        # Search built-in pathways
        for source, pathways in self._pathways.items():
            for name, genes in pathways.items():
                if query_lower in name.lower():
                    results.append({
                        'name': name,
                        'source': source,
                        'size': len(genes),
                        'genes': genes
                    })

        return results

    def get_pathway_genes(self, pathway_name: str) -> Set[str]:
        """
        Get genes in a pathway.

        Args:
            pathway_name: Name of the pathway

        Returns:
            Set of gene symbols
        """
        for source, pathways in self._pathways.items():
            if pathway_name in pathways:
                return pathways[pathway_name]
        return set()

    def list_all_pathways(self) -> List[str]:
        """Get list of all available pathway names."""
        all_pathways = []
        for source, pathways in self._pathways.items():
            all_pathways.extend(pathways.keys())
        return sorted(all_pathways)

    def fetch_kegg_pathway(self, pathway_id: str) -> Optional[Set[str]]:
        """
        Fetch genes from KEGG pathway.

        Args:
            pathway_id: KEGG pathway ID (e.g., 'hsa04310' for Wnt signaling)

        Returns:
            Set of gene symbols or None if failed
        """
        cache_file = self.cache_dir / f'kegg_{pathway_id}.json'

        # Check cache
        if cache_file.exists():
            with open(cache_file) as f:
                return set(json.load(f))

        try:
            url = f'https://rest.kegg.jp/get/{pathway_id}'
            response = requests.get(url, timeout=10)

            if response.status_code == 200:
                genes = set()
                in_gene_section = False

                for line in response.text.split('\n'):
                    if line.startswith('GENE'):
                        in_gene_section = True
                        # Extract gene from first GENE line
                        parts = line.split()
                        if len(parts) >= 3:
                            genes.add(parts[2].rstrip(';'))
                    elif in_gene_section and line.startswith(' '):
                        parts = line.split()
                        if len(parts) >= 2:
                            genes.add(parts[1].rstrip(';'))
                    elif not line.startswith(' ') and in_gene_section:
                        break

                # Cache results
                with open(cache_file, 'w') as f:
                    json.dump(list(genes), f)

                return genes
        except Exception:
            pass

        return None

    def get_pathway_overlap(
        self,
        pathway_genes: Set[str],
        deg_genes: Set[str]
    ) -> Dict:
        """
        Calculate overlap between pathway and DEG set.

        Args:
            pathway_genes: Genes in pathway
            deg_genes: Differentially expressed genes

        Returns:
            Dictionary with overlap statistics
        """
        overlap = pathway_genes & deg_genes

        return {
            'pathway_size': len(pathway_genes),
            'deg_count': len(deg_genes),
            'overlap_count': len(overlap),
            'overlap_genes': overlap,
            'overlap_fraction': len(overlap) / len(pathway_genes) if pathway_genes else 0
        }

    def add_custom_pathway(self, name: str, genes: Set[str]):
        """Add a custom pathway to the database."""
        if 'custom' not in self._pathways:
            self._pathways['custom'] = {}
        self._pathways['custom'][name] = genes
