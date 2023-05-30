# API

## Enrichment
 
### ***gopher.test_enrichment***
    test_enrichment(proteins, desc, aspect, species, release, go_subset, contaminants_filter, fetch, progress, annotations, mapping, aggregate_terms, alternative)

Test for the enrichment of Gene Ontology terms from protein abundance.

The Mann-Whitney U Test is applied to each column of proteins dataframe and
for each Gene Ontology (GO) term. The p-values are then corrected for
multiple hypothesis testing across all of the columns using the
Benjamini-Hochberg procedure.

***Parameters***

> proteins : pandas.DataFrame

A dataframe where the indices are UniProtKB accessions and each column is
an experiment to test. The values in this dataframe should be some
measure of protein abundance: these could be the raw measurement if
originating from a single sample or a fold-change/p-value if looking at
the difference between two conditions.

> desc : bool, *default = True*

Rank proteins in descending order?

> aspect : str, {"cc", "mf", "bp", "all"}, *default = "all"*

The Gene Ontology aspect to use. Use "cc" for "Cellular Compartment",
"mf" for "Molecular Function", "bp" for "Biological Process", or "all"
for all three.

> species : str, {"human", "yeast", ...}, *default = "human"*

The species for which to retrieve GO annotations. If not "human" or
"yeast", see
[here](http://current.geneontology.org/products/pages/downloads.html).

> release : str, *default = "current"*

The Gene Ontology release version. Using "current" will look up the
most current version.

> go_subset: List[str], *default = None*

The go terms of interest. Should consists of the go term names such
as 'nucleus' or 'cytoplasm'.

> contaminants_filter: List[str], *default = None*

A list of UniProtKB accessions for common contaminants such as
Keratin to filter out.

> fetch : bool, *default = False*

Download the GO annotations even if they have been downloaded before?

> progress : bool, *default = False*

Show a progress bar during enrichment tests?

> annotations: pandas.DataFrame, *default = None*

A custom annotations dataframe.

> mapping: defaultdict, *default = None*

A custom mapping of the GO term relationships.

> aggregate_terms : bool, *default = True*

Aggregate the terms and do the graph search.

> alternative : str, {"greater", "less", "two-sided"}, *default = "greater"*

Type of test that should be run.
Could be "greater", "less", or "two-sided".

***Returns***
> pandas.DataFrame

The adjusted p-value for each tested GO term in each sample.
            


## GO Annotations

### ***gopher.generate_annotations***

    generate_annotations(proteins, aspect, go_name, go_id)

Generate an annotation file for a list of proteins that are correlated to a single term and aspect.

The term can be in the GO database or a new term.

***Parameters***
> proteins : List[str]

List of proteins (UniProtKB accessions) that will be annotated to a term.

> aspect: str

String specifying the aspect the term is in ("C", "F", "P").

> go_name : str

String of the GO name for the proteins

> go_id : str, *default = None*

String of the GO ID. If in the GO database, the go id and go name should match the database.

***Returns***
> pandas.DataFrame

An annotations dataframe with a single go term.

### ***gopher.download_annotations***

    download_annotations(stem, release, fetch)

Download the annotation file.

See http://current.geneontology.org/annotations/index.html for details.

***Parameters***
> stem : str

The stem of the annotation file name to retrieve.

> release: str, *default = "current"*

The Gene Ontology release version. Using "current" will look up the
most current version.

> fetch : bool, *default = False*

Download the file even if it already exists?

***Returns***
> pathlib.Path

The path to downloaded file.

## Normalization

### ***gopher.normalize_values***

    normalize_values(proteins, fasta)

Normalize intensity values.

Normalize using the proteomic ruler approach outlined by Wiśniewski et al.
(doi: https://doi.org/10.1074/mcp.M113.037309)

***Parameters***

> proteins : pandas.DataFrame

A dataframe where the indices are UniProtKB accessions and each column is
an experiment to test. The values in this dataframe raw protein abundance

> fasta : str or pathlib.Path

Use the FASTA file to generate molecular weights for normalization

***Returns***

> pandas.DataFrame

The normalized intensities for every protein in each sample.


## Parsers

### ***gopher.read_encyclopedia***

    read_encyclopedia(proteins_txt)

Read results from EncyclopeDIA.

***Parameters***

> proteins_txt : str

The EncyclopeDIA protein output.

***Returns***

> pandas.DataFrame

The EncyclopeDIA results in a format for gopher.

### ***gopher.read_metamorpheus***

    read_metamorpheus(proteins_txt)

Read results from Metamorpheus.

***Parameters***

> proteins_txt : str

The Metamorpheus protein output file.

***Returns***

> pandas.DataFrame

The Metamorpheus results in a format for gopher.


## Config

### ***gopher.get_data_dir***

    get_data_dir()

Retrieve the current data directory for ppx

### ***gopher.set_data_dir***

    set_data_dir(path):

Set the ppx data directory.

***Parameters***

> path : str or pathlib.Path object, *default = None*

The path for ppx to use as its data directory.


## Other

### ***gopher.map_proteins***

    map_proteins(protein_list, aspect, species, release, fetch)

Map the proteins to the GO terms

***Parameters***

> protein_list : List[str]

A list of UniProtKB accessions.

> aspect : str, {"cc", "mf", "bp", "all"}, *default = "all"*

The Gene Ontology aspect to use. Use "cc" for "Cellular Compartment",
"mf" for "Molecular Function", "bp" for "Biological Process", or "all"
for all three.

> species : str, {"human", "yeast", ...}, *default = "human"*

The species for which to retrieve GO annotations. If not "human" or
"yeast", see
[here](http://current.geneontology.org/products/pages/downloads.html).

> release : str, *default = "current"*

The Gene Ontology release version. Using "current" will look up the
most current version.

> fetch : bool, *default = False*

Download the GO annotations even if they have been downloaded before?

***Returns***

> pandas.DataFrame

Dataframe with protein accessions and GO terms


### ***gopher.get_rankings***

    get_rankings(proteins, go_term, aspect, species, release, fetch)

Rank the proteins and show whether proteins are in a specified term

***Parameters***

> proteins : pd.DataFrame

Dataframe of protein quant data

> go_term : str

String of specified GO term name

> aspect : str, {"cc", "mf", "bp", "all"}, *default = "all"*

The Gene Ontology aspect to use. Use "cc" for "Cellular Compartment",
"mf" for "Molecular Function", "bp" for "Biological Process", or "all"
for all three.

> species : str, {"human", "yeast", ...}, *default = "human"*

The species for which to retrieve GO annotations. If not "human" or
"yeast", see
[here](http://current.geneontology.org/products/pages/downloads.html).

> release : str, *default = "current"*

The Gene Ontology release version. Using "current" will look up the
most current version.

> fetch : bool, *default = False*

Download the GO annotations even if they have been downloaded before?

***Returns***

> pandas.DataFrame

Dataframe with protein rankings and whether or not the protein is in the specified term


### ***gopher.get_annotations***

    get_annotations(proteins, go_term, aspect, species, release, fetch, go_subset)

Gets the annotations for proteins in a dataset

***Parameters***

> proteins : pd.DataFrame

Dataframe of protein quant data

> go_term : str

String of specified GO term name

> aspect : str, {"cc", "mf", "bp", "all"}, *default = "all"*

The Gene Ontology aspect to use. Use "cc" for "Cellular Compartment",
"mf" for "Molecular Function", "bp" for "Biological Process", or "all"
for all three.

> species : str, {"human", "yeast", ...}, *default = "human"*

The species for which to retrieve GO annotations. If not "human" or
"yeast", see
[here](http://current.geneontology.org/products/pages/downloads.html).

> release : str, *default = "current"*

The Gene Ontology release version. Using "current" will look up the
most current version.

> fetch : bool, *default = False*

Download the GO annotations even if they have been downloaded before?

> go_subset: List[str], *default = None*

The go terms of interest. Should consists of the go term names such
as 'nucleus' or 'cytoplasm'.

***Returns***

> pandas.DataFrame

Dataframe with protein annotations


### ***gopher.in_term***

    in_term(proteins, go_term, annot)

See if proteins are associated with a specific term

***Parameters***

> proteins : pd.DataFrame

Dataframe of proteins and quantifications

> go_term : str

String of specified GO term name

> annot : pd.DataFrame

Annotation file for the dataset

***Returns***

> pandas.DataFrame

Dataframe with protein quant and if protein is in the given term


### ***gopher.roc***

    roc(proteins, go_term, aspect, species, release, fetch)

Plot the ROC curve for a go term in each sample

***Parameters***

> proteins : pd.DataFrame

Dataframe of proteins and quantifications

> go_term : str

String of specified GO term name
> aspect : str, {"cc", "mf", "bp", "all"}, *default = "all"*

The Gene Ontology aspect to use. Use "cc" for "Cellular Compartment",
"mf" for "Molecular Function", "bp" for "Biological Process", or "all"
for all three.

> species : str, {"human", "yeast", ...}, *default = "human"*

The species for which to retrieve GO annotations. If not "human" or
"yeast", see
[here](http://current.geneontology.org/products/pages/downloads.html).

> release : str, *default = "current"*

The Gene Ontology release version. Using "current" will look up the
most current version.

> fetch : bool, *default = False*

Download the GO annotations even if they have been downloaded before?

***Returns***

> matplotlib.pyplot

Plot of ROC curve for a GO term