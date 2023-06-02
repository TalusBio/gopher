# gopher

Gene ontology enrichment analysis using protein expression.

Gopher uses a Mann-Whitney U Test to look for enriched gene ontology terms that are present when proteins are ranked by a quantitative value, such as the difference between two conditions or the abundance of the protein.

See the documentation at https://TalusBio.github.io/gopher

## Installation

Gopher can be pip installed using:

``` sh
pip install gopher-enrich
```

## TLDR

```python
import gopher                                   # import the package

gopher.test_enrichment(proteins=protein_quant)  # run the enrichment 
```