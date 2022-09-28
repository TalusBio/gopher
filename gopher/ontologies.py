"""Download the GO ontologies"""
from . import utils, config


def download_ontology():
    """download the Gene Ontology terms.

    Returns
    -------
    Path
        The downloaded OBO file.
    """
    url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
    out_file = config.get_data_dir() / "ontologies" / "go-basic.obo"
    if out_file.exists():
        return out_file

    out_file.parent.mkdir(exist_ok=True, parents=True)
    utils.http_download(url, out_file)
    return out_file


def load_ontology():
    """Load the Gene Ontology terms.

    We use the basic version of GO.

    Returns
    -------
    dict of str: str
        The GO accession mapped to the name.
    """
    obo_file = download_ontology()
    with obo_file.open("r") as obo_ref:
        data = obo_ref.read().split("\n\n")[1:]

    terms = {}
    mapping = {}
    for term in data:
        term_data = term.splitlines()
        term_id, term_name = None, None
        for line in term_data:
            try:
                key, val = line.split(": ", 1)
            except ValueError:
                continue

            if key == "id":
                term_id = val
            elif key == "name":
                term_name = val
            elif key == "is_a":
                val = val.split(" ", 1)[0]
                if val in mapping.keys():
                    mapping[val].append(term_id)
                else:
                    mapping[val] = [term_id]

            if term_id is not None and term_name is not None:
                terms[term_id] = term_name

    return terms, mapping
