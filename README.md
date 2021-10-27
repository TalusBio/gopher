# gopher

Gene ontology enrichment analysis using protein expression.

Try it. It's super slow right now:

``` python
import gopher
df = gopher.read_encyclopedia("RESULTS-quant.elib.proteins.txt")
results = gopher.test_enrichment(df)
```
