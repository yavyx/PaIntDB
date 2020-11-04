from collections import defaultdict
import json
import os
import pickle
import sqlite3

import pandas as pd

DB_PATH = 'PaIntDB.db'
ONTOLOGY_PATH = os.path.join('data', 'PAO1_gene_ontology.csv')


def make_metabolite_mapping():
    """Generates a JSON dictionary with metabolites mapped to their corresponding gene,"""
    with sqlite3.connect(DB_PATH) as db_connection:
        cursor = db_connection.cursor()
        cursor.execute('SELECT id, kegg, pubchem, cas, chebi, ecocyc FROM metabolite')
        gene_metabolite_map = defaultdict(list)  # Creates a dictionary with empty lists as default values.


def make_interactome(strain):
    """Generates interactomes to use with OmicsIntegrator."""
    with sqlite3.connect(DB_PATH) as db_connection:
        cursor = db_connection.cursor()
        # Query database with selected strain
        cursor.execute("""SELECT interactor_id, interaction.id, type, is_experimental
                                  FROM interaction_participants
                                  INNER JOIN interaction_sources
                                  USING (interaction_id)
                                  INNER JOIN interaction_source
                                  ON interaction_sources.data_source = interaction_source.id
                                  INNER JOIN interaction
                                  ON interaction_id = interaction.id
                                  WHERE strain = ?
                                  AND type = 'p-p'""", [strain])
        interactors = cursor.fetchall()
        interactome = dict()
        # Create a dictionary edge list from the list of interactions (two rows per interaction)
        for i in range(0, len(interactors), 2):
            interactome[interactors[i][1]] = (interactors[i][0],  # 1st interactor
                                              interactors[i+1][0],  # 2nd interactor
                                              interactors[i][3])  # Confidence
        # Create data frame from dictionary
        interactome_df = pd.DataFrame.from_dict(interactome, orient='index',
                                                columns=['protein1', 'protein2', 'confidence'])
        # Convert confidence values into edge costs
        interactome_df['cost'] = 1.5 - interactome_df['confidence'] / 2
        del interactome_df['confidence']
        interactome_df.to_csv(os.path.join('data', '{}_interactome.tsv'.format(strain)), sep='\t')
    return interactome_df


def make_go_association_dict():
    """Creates and pickles a nested GO association dictionary that maps genes to GO terms which is used by GOAtools
    for GO enrichment, using the Ontology file in data directory."""

    def create_dict(row):
        """Assign gene in DataFrame row to corresponding GO term in dictionary."""
        locus_tag, go_id, namespace = row['Locus Tag'], row['Accession'], row['Namespace']
        if domain == namespace:
            if locus_tag not in go_dict[domain_short]:
                go_dict[domain_short][locus_tag] = set()  # Use sets to avoid duplicate GO ID's
                go_dict[domain_short][locus_tag].add(go_id)
            else:
                go_dict[domain_short][locus_tag].add(go_id)

    go_dict = {'BP': dict(), 'CC': dict(), 'MF': dict()}
    go_domains = {'BP': 'biological_process',
                  'CC': 'cellular_component',
                  'MF': 'molecular_function'}

    go_associations = pd.read_csv(ONTOLOGY_PATH)

    for domain_short, domain in go_domains.items():
        go_associations.apply(create_dict, 'columns')

    # Save dictionary as pickle
    with open(os.path.join('data', 'go_association.pickle'), 'wb') as f:
        pickle.dump(go_dict, f)

    return go_dict
