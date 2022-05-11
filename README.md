# gopher

Gene ontology enrichment analysis using protein expression.

Gopher uses a Mann-Whitney U Test to look for enriched gene ontology terms that are present when proteins are ranked by a quantitative value, such as the difference between two conditions or the abundance of the protein.

## Installation

Gopher can be pip installed directly from GitHub:

``` sh
pip install git+https://github.com/TalusBio/gopher.git
```

## Usage

Gopher can be used within Python. 
Currently, Gopher is pretty slow---we recommend making use of the `go_subset` argument to limit the number of GO terms considered by Gopher. 
Here is an example:

``` python
# Import packages
import gopher
import vizta
import matplotlib.pyplot as plt
import seaborn as sns

# Use Talus Bio's plotting theme:
vizta.mpl.set_theme()

# Only test a few cellular compartment GO terms:
terms = [
    "nucleus",
    "nuclear chromosome",
    "nucleoplasm",
    "euchromatin",
    "heterochromatin",
    "protein-DNA complex",
    "transcription regulator complex",
    "inner mitochondrial membrane protein complex",
    "mitochondrial nucleoid",
    "cell surface",
    "ER to Golgi transport vesicle membrane",
    "organelle membrane",
    "lysosome",
    "cytoplasm",
]

# Read data from an EncyclopeDIA output:
proteins = gopher.read_encyclopedia("RESULTS-quant.elib.proteins.txt")

# Peform the GO enrichment analysis:
results = gopher.test_enrichment(
    proteins=df,
    aspect="cc",
    go_subset=terms,
    progress=True,
)

# Create a plot:
plt.figure(figsize=(8, 3))
sns.barplot(data=results_df, x="GO Name", y="pvalue (-log10)", hue="Run")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()
```
