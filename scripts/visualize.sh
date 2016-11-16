#!/bin/bash

python scripts/visualize_genomic_elements.py --TAD-Boundary 'hESC' 
python scripts/visualize_genomic_elements.py --TAD-Boundary 'IMR90'
python scripts/visualize_genomic_elements.py --TAD-Boundary 'mESC'
python scripts/visualize_genomic_elements.py --TAD-Boundary 'cortex'

python scripts/visualize_gc_and_divergence.py --TAD-Boundary 'hESC'
python scripts/visualize_gc_and_divergence.py --TAD-Boundary 'IMR90'
python scripts/visualize_gc_and_divergence.py --TAD-Boundary 'mESC'
python scripts/visualize_gc_and_divergence.py --TAD-Boundary 'cortex'

python scripts/visualize_gwas_distribution.py
