"""
Gene information fetching and caching for tooltips.

Fetches gene metadata from public databases for context tooltips.
"""

import requests
from typing import Dict, Optional
import streamlit as st


class GeneInfoFetcher:
    """Fetch and cache gene information from public APIs."""

    def __init__(self):
        """Initialize with caching."""
        self.cache = {}

    @st.cache_data(ttl=86400, show_spinner=False)  # Cache for 24 hours
    def get_gene_info(_self, gene_symbol: str) -> Optional[Dict]:
        """
        Fetch gene information from MyGene.info API.

        Args:
            gene_symbol: Gene symbol (e.g., "ANXA2")

        Returns:
            Dictionary with gene info or None if not found
        """
        # Check memory cache first
        if gene_symbol in _self.cache:
            return _self.cache[gene_symbol]

        try:
            # Use MyGene.info API (free, no key required)
            url = f"https://mygene.info/v3/query"
            params = {
                'q': gene_symbol,
                'species': 'human',
                'fields': 'symbol,name,summary,pathway,go,interpro,uniprot',
                'size': 1
            }

            response = requests.get(url, params=params, timeout=3)
            response.raise_for_status()
            data = response.json()

            if 'hits' in data and len(data['hits']) > 0:
                gene_data = data['hits'][0]

                # Extract relevant info
                info = {
                    'symbol': gene_data.get('symbol', gene_symbol),
                    'name': gene_data.get('name', 'Unknown'),
                    'summary': gene_data.get('summary', 'No summary available'),
                    'pathways': [],
                    'go_terms': [],
                    'uniprot': None
                }

                # Extract pathways
                if 'pathway' in gene_data:
                    pathways = gene_data['pathway']
                    if isinstance(pathways, dict):
                        if 'kegg' in pathways:
                            kegg = pathways['kegg']
                            if isinstance(kegg, list):
                                info['pathways'] = [p.get('name', '') for p in kegg[:3]]
                            elif isinstance(kegg, dict):
                                info['pathways'] = [kegg.get('name', '')]

                # Extract GO terms (molecular function)
                if 'go' in gene_data:
                    go_data = gene_data['go']
                    if isinstance(go_data, dict) and 'MF' in go_data:
                        mf = go_data['MF']
                        if isinstance(mf, list):
                            info['go_terms'] = [term.get('term', '') for term in mf[:3]]
                        elif isinstance(mf, dict):
                            info['go_terms'] = [mf.get('term', '')]

                # Extract UniProt ID
                if 'uniprot' in gene_data:
                    uniprot = gene_data['uniprot']
                    if isinstance(uniprot, dict):
                        info['uniprot'] = uniprot.get('Swiss-Prot', '')
                    elif isinstance(uniprot, str):
                        info['uniprot'] = uniprot

                # Cache result
                _self.cache[gene_symbol] = info
                return info

        except Exception as e:
            # Silently fail - don't disrupt user experience
            return None

        return None

    def get_quick_context(self, gene_symbol: str) -> str:
        """
        Get a one-line context string for a gene.

        Args:
            gene_symbol: Gene symbol

        Returns:
            Brief context string
        """
        info = self.get_gene_info(gene_symbol)

        if not info:
            return f"{gene_symbol}"

        # Build quick context
        parts = [info['symbol']]

        # Add function if available
        if info['go_terms']:
            parts.append(info['go_terms'][0][:40])  # Truncate long terms

        return " - ".join(parts)

    def format_tooltip(self, gene_symbol: str) -> str:
        """
        Format gene info as HTML tooltip content.

        Args:
            gene_symbol: Gene symbol

        Returns:
            HTML string for tooltip
        """
        info = self.get_gene_info(gene_symbol)

        if not info:
            return f"""
            <div style="padding: 8px; min-width: 200px;">
                <strong>{gene_symbol}</strong><br>
                <span style="color: #888;">No information available</span>
            </div>
            """

        # Build HTML
        html = f"""
        <div style="padding: 12px; min-width: 300px; max-width: 400px; font-size: 13px;">
            <strong style="font-size: 15px;">{info['symbol']}</strong><br>
            <span style="color: #666; font-style: italic;">{info['name']}</span>
            <hr style="margin: 8px 0; border: none; border-top: 1px solid #ddd;">
        """

        # Add GO terms
        if info['go_terms']:
            html += f"""
            <div style="margin-top: 8px;">
                <strong style="color: #2196F3;">Function:</strong><br>
                <span style="font-size: 12px;">{info['go_terms'][0][:100]}</span>
            </div>
            """

        # Add pathways
        if info['pathways']:
            html += f"""
            <div style="margin-top: 8px;">
                <strong style="color: #4CAF50;">Pathways:</strong><br>
                <span style="font-size: 12px;">{"<br>".join(['• ' + p[:60] for p in info['pathways'][:2]])}</span>
            </div>
            """

        # Add links
        links = []
        if info['uniprot']:
            links.append(f'<a href="https://www.uniprot.org/uniprot/{info["uniprot"]}" target="_blank">UniProt</a>')
        links.append(f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene_symbol}" target="_blank">GeneCards</a>')
        links.append(f'<a href="https://pubmed.ncbi.nlm.nih.gov/?term={gene_symbol}" target="_blank">PubMed</a>')

        if links:
            html += f"""
            <div style="margin-top: 12px; padding-top: 8px; border-top: 1px solid #ddd;">
                <span style="font-size: 11px;">{' | '.join(links)}</span>
            </div>
            """

        html += "</div>"
        return html
