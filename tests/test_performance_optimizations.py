"""
Test suite for performance optimizations in gene search.

Tests for:
- Gene indexing
- Fuzzy matching
- Caching functions
- Search history
"""

import pytest
import pandas as pd
import numpy as np
from unittest.mock import MagicMock
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


def test_build_gene_index():
    """Test that gene index is built correctly for O(1) lookups."""
    from app import build_gene_index

    # Create mock datasets
    datasets = {
        'Dataset1': pd.DataFrame({
            'Gene': ['ANXA2', 'LARP4', 'IL17RC'],
            'log2FC': [1.5, -2.0, 0.5],
            'pvalue': [0.01, 0.001, 0.05]
        }),
        'Dataset2': pd.DataFrame({
            'Gene': ['ANXA2', 'CMTM4', 'AIFM2'],
            'log2FC': [1.2, -1.8, 0.3],
            'pvalue': [0.02, 0.005, 0.1]
        })
    }

    # Build index (returns tuple: gene_index, sorted_genes)
    gene_index, sorted_genes = build_gene_index(datasets)

    # Verify sorted_genes list
    assert 'ANXA2' in sorted_genes
    assert sorted_genes == sorted(sorted_genes)  # Should be pre-sorted

    # Verify index structure
    assert 'ANXA2' in gene_index
    assert 'Dataset1' in gene_index['ANXA2']
    assert 'Dataset2' in gene_index['ANXA2']
    assert gene_index['ANXA2']['Dataset1'] == 0  # First row
    assert gene_index['ANXA2']['Dataset2'] == 0

    # Verify all genes are indexed
    assert 'LARP4' in gene_index
    assert 'IL17RC' in gene_index
    assert 'CMTM4' in gene_index
    assert 'AIFM2' in gene_index

    # Verify row indices are correct
    assert gene_index['LARP4']['Dataset1'] == 1
    assert gene_index['IL17RC']['Dataset1'] == 2

    print("✓ Gene index built correctly")


def test_fuzzy_match_genes():
    """Test fuzzy matching functionality for typo tolerance."""
    from app import fuzzy_match_genes

    gene_list = ['ANXA2', 'ANXA1', 'LARP4', 'LARP4B', 'IL17RC', 'IL17RA']

    # Test exact match
    matches = fuzzy_match_genes('ANXA2', gene_list, threshold=80)
    assert len(matches) > 0
    assert matches[0][0] == 'ANXA2'
    assert matches[0][1] == 100  # Perfect match

    # Test typo tolerance (ANAX2 -> ANXA2)
    matches = fuzzy_match_genes('ANAX2', gene_list, threshold=60)
    assert len(matches) > 0
    assert 'ANXA2' in [m[0] for m in matches]

    # Test starts-with matching
    matches = fuzzy_match_genes('LARP', gene_list, threshold=60)
    assert len(matches) >= 2
    larp_genes = [m[0] for m in matches]
    assert 'LARP4' in larp_genes

    # Test case insensitivity
    matches = fuzzy_match_genes('anxa2', gene_list, threshold=80)
    assert matches[0][0] == 'ANXA2'

    print("✓ Fuzzy matching works correctly")


def test_search_history():
    """Test search history tracking."""
    from app import add_to_search_history

    # Mock session state
    class MockSessionState:
        def __init__(self):
            self.search_history = []

    import app
    app.st = MagicMock()
    app.st.session_state = MockSessionState()

    # Add genes to history
    add_to_search_history('ANXA2', max_history=5)
    add_to_search_history('LARP4', max_history=5)
    add_to_search_history('IL17RC', max_history=5)

    # Check order (most recent first)
    assert app.st.session_state.search_history[0] == 'IL17RC'
    assert app.st.session_state.search_history[1] == 'LARP4'
    assert app.st.session_state.search_history[2] == 'ANXA2'

    # Test deduplication (adding same gene moves it to front)
    add_to_search_history('ANXA2', max_history=5)
    assert app.st.session_state.search_history[0] == 'ANXA2'
    assert len(app.st.session_state.search_history) == 3  # No duplicates

    # Test max history limit
    for i in range(10):
        add_to_search_history(f'GENE{i}', max_history=5)
    assert len(app.st.session_state.search_history) == 5

    print("✓ Search history tracking works correctly")


def test_datasets_hash():
    """Test dataset hash generation for cache invalidation."""
    from app import get_datasets_hash

    # Mock session state with datasets
    class MockSessionState:
        def __init__(self):
            self.datasets = {
                'Dataset1': pd.DataFrame({
                    'Gene': ['A', 'B'],
                    'log2FC': [1, 2],
                    'pvalue': [0.01, 0.02]
                })
            }

    import app
    app.st = MagicMock()
    app.st.session_state = MockSessionState()

    # Get initial hash
    hash1 = get_datasets_hash()

    # Same data should produce same hash
    hash2 = get_datasets_hash()
    assert hash1 == hash2

    # Different data should produce different hash
    app.st.session_state.datasets['Dataset2'] = pd.DataFrame({
        'Gene': ['C', 'D'],
        'log2FC': [3, 4],
        'pvalue': [0.03, 0.04]
    })
    hash3 = get_datasets_hash()
    assert hash1 != hash3

    print("✓ Dataset hashing works correctly")


def test_gene_index_performance():
    """Test that gene index provides significant speedup."""
    import time
    from app import build_gene_index

    # Create large mock dataset
    n_genes = 20000
    genes = [f'GENE{i}' for i in range(n_genes)]

    datasets = {
        'Dataset1': pd.DataFrame({
            'Gene': genes,
            'log2FC': np.random.randn(n_genes),
            'pvalue': np.random.rand(n_genes)
        })
    }

    # Build index (returns tuple: gene_index, sorted_genes)
    start = time.time()
    gene_index, sorted_genes = build_gene_index(datasets)
    index_time = time.time() - start

    # Verify sorted_genes is pre-sorted
    assert len(sorted_genes) == n_genes
    assert sorted_genes == sorted(sorted_genes)

    # Test lookup using index (O(1))
    start = time.time()
    for _ in range(1000):
        gene_to_find = f'GENE{np.random.randint(0, n_genes)}'
        if gene_to_find in gene_index:
            idx = gene_index[gene_to_find].get('Dataset1', -1)
    indexed_time = time.time() - start

    # Test lookup using scan (O(n))
    df = datasets['Dataset1']
    start = time.time()
    for _ in range(1000):
        gene_to_find = f'GENE{np.random.randint(0, n_genes)}'
        result = df[df['Gene'] == gene_to_find]
    scan_time = time.time() - start

    # Index lookup should be significantly faster
    speedup = scan_time / indexed_time
    print(f"  Index build time: {index_time:.3f}s")
    print(f"  Indexed lookup time: {indexed_time:.3f}s")
    print(f"  Scan lookup time: {scan_time:.3f}s")
    print(f"  Speedup: {speedup:.1f}x")

    assert speedup > 5, f"Expected >5x speedup, got {speedup:.1f}x"
    print("✓ Gene index provides significant speedup")


if __name__ == '__main__':
    print("Running performance optimization tests...\n")

    test_build_gene_index()
    test_fuzzy_match_genes()
    test_search_history()
    test_datasets_hash()
    test_gene_index_performance()

    print("\n✅ All tests passed!")
