from collections import defaultdict
from datetime import datetime
import os
import pickle
import sqlite3

import pandas as pd

DB_PATH = 'PaIntDB.db'


def make_metabolite_mapping():
    """Generates a JSON dictionary with metabolites mapped to their corresponding gene,"""
    with sqlite3.connect(DB_PATH) as db_connection:
        cursor = db_connection.cursor()
        cursor.execute('SELECT id, kegg, pubchem, cas, chebi, ecocyc FROM metabolite')
        gene_metabolite_map = defaultdict(list)  # Creates a dictionary with empty lists as default values.


def make_go_association_dict(path):
    """Creates and pickles a nested GO association dictionary that maps genes to GO terms which is used by GOAtools
    for GO enrichment, using the Ontology file in data directory."""

    def create_dict(row):
        """Assign gene in DataFrame row to corresponding GO term in dictionary."""
        locus_tag, go_id, namespace = row['Locus Tag'], row['Accession'], row['Namespace']
        if domain == namespace:
            if locus_tag not in go_dict[domain_short]:
                go_dict[domain_short][locus_tag] = set()
                go_dict[domain_short][locus_tag].add(go_id)
            else:
                go_dict[domain_short][locus_tag].add(go_id)

    go_dict = {'BP': dict(), 'CC': dict(), 'MF': dict()}

    go_domains = {'BP': 'biological_process',
                  'CC': 'cellular_component',
                  'MF': 'molecular_function'}

    go_associations = pd.read_csv(path)

    for domain_short, domain in go_domains.items():
        go_associations.apply(create_dict, 'columns')

    with open(os.path.join('data', 'go_association.pickle'), 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(go_dict, f, pickle.HIGHEST_PROTOCOL)

    return go_dict
