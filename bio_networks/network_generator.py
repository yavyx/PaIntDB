import os
import sqlite3

import networkx as nx
import pandas as pd

import bio_networks.helpers as helpers

DB_PATH = '/home/javier/PycharmProjects/PaIntDB/PaIntDB.db'


class BioNetwork:
    """Creates NetworkX networks with additional biological attributes for use with PaintDB."""

    def __init__(self, gene_list_path, strain, order, detection_method, metabolites=False):
        self._strain = strain
        self._order = order
        self._detection_method = detection_method
        self._metabolites = metabolites
        self._genes_of_interest = helpers.get_genes(gene_list_path)
        self._raw_info = dict()
        self._interactions_of_interest = []
        self._network = nx.Graph()

    def query_db(self):
        """Queries PaintDB depending on the selected filters and adds the raw information to the network."""
        detection_methods = {'computational': 0, 'mixed': 1, 'experimental': 2, 'all': 3}

        with sqlite3.connect(DB_PATH) as db_connection:
            cursor = db_connection.cursor()
            if self._metabolites is True:
                cursor.execute('SELECT id, kegg, pubchem, cas, chebi, ecocyc FROM metabolite')
                self._raw_info['metabolite'] = cursor.fetchall()
                interaction_type = "'p-p', 'p-m', 'm-p'"
            else:
                interaction_type = 'p-p'

            if detection_methods[self._detection_method] in [0, 1, 2]:
                # Node info (lists to generate dictionaries later)
                cursor.execute("""SELECT interactor_id, interaction.id, type, is_experimental
                                  FROM interaction_participants
                                  INNER JOIN interaction_sources
                                  USING (interaction_id)
                                  INNER JOIN interaction_source
                                  ON interaction_sources.data_source = interaction_source.id
                                  INNER JOIN interaction
                                  ON interaction_id = interaction.id
                                  WHERE strain = ?
                                  AND is_experimental = ?
                                  AND type IN (?)""",
                               [self._strain, detection_methods[self._detection_method], interaction_type])
                self._raw_info['interaction_participants'] = cursor.fetchall()

                # Edge info (dataFrames to merge with the edge list dataFrame)
                self._raw_info['sources'] = pd.read_sql_query("""SELECT is_experimental, interaction_id
                                                                 FROM interaction_source 
                                                                 INNER JOIN interaction_sources
                                                                 ON interaction_source.id = 
                                                                 interaction_sources.data_source
                                                                 WHERE is_experimental = ?""",
                                                              con=db_connection,
                                                              params=[detection_methods[self._detection_method]])

            elif detection_methods[self._detection_method] == 3:
                # Node info (lists to generate dictionaries later)
                cursor.execute("""SELECT interactor_id, interaction_id, type
                                  FROM interaction_participants 
                                  INNER JOIN interaction
                                  ON interaction_participants.interaction_id = interaction.id 
                                  WHERE strain = ?
                                  AND type = ?""",
                               [self._strain])
                self._raw_info['interaction_participants'] = cursor.fetchall()

                # Edge info (dataFrames to merge with the edge list dataFrame)
                self._raw_info['sources'] = pd.read_sql_query("""SELECT is_experimental, interaction_id
                                                                 FROM interaction_source 
                                                                 INNER JOIN interaction_sources
                                                                 ON interaction_source.id = 
                                                                 interaction_sources.data_source""",
                                                              con=db_connection)

            cursor.execute("""SELECT id, product_name, ncbi_acc, uniprotkb 
                              FROM protein
                              WHERE strain = ?""",
                           [self._strain])
            self._raw_info['proteins'] = cursor.fetchall()

            cursor.execute('SELECT id, name, type FROM interactor')
            self._raw_info['short_names'] = cursor.fetchall()

            cursor.execute("""SELECT protein_id, localization
                              FROM localization
                              INNER JOIN protein_localizations
                              ON id = localization_id""")
            self._raw_info['localization'] = cursor.fetchall()
        # Change sources ID's to numeric
        self._raw_info['sources'].interaction_id = pd.to_numeric(self._raw_info['sources'].interaction_id,
                                                                 downcast='integer')

    def format_attribute_dictionaries(self):
        """Returns nested dictionaries of node attributes that can be added directly to a NetworkX graph."""
        db_tidy_info = dict()

        proteins = self._raw_info['proteins']
        node_protein_info = dict()
        for protein in proteins:
            node_protein_info[protein[0]] = dict(description=protein[1],
                                                 ncbi_acc=protein[2],
                                                 uniprotkb=protein[3])
        for key, value in node_protein_info.items():
            node_protein_info[key] = helpers.remove_nones(value)
        db_tidy_info['proteins'] = node_protein_info

        localizations = self._raw_info['localization']
        node_localizations = dict()
        for localization in localizations:
            node_localizations[localization[0]] = localization[1]
        db_tidy_info['localization'] = node_localizations

        names = self._raw_info['short_names']
        short_names = dict()
        for name in names:
            short_names[name[0]] = dict(short_name=name[1],
                                        type=name[2])
        for key, value in short_names.items():
            short_names[key] = helpers.remove_nones(value)
        db_tidy_info['short_names'] = short_names

        return db_tidy_info

    def make_edge_list(self):
        """Returns a Pandas edge list dataFrame that can be directly used to generate a network with NetworkX, and
        filters the genes of interest."""
        interaction_participants = self._raw_info['interaction_participants']
        interaction_edges = dict()
        # Create a dictionary edge list from the list of interactions (two rows per interaction)
        for i in range(0, len(interaction_participants), 2):
            interaction_edges[interaction_participants[i][1]] = (interaction_participants[i][0],
                                                                 interaction_participants[i+1][0],
                                                                 interaction_participants[i][2])
        if self._order == 0:
            interactions_of_interest = [interactionID for interactionID, genes in interaction_edges.items()
                                        if genes[0] in self._genes_of_interest and genes[1] in self._genes_of_interest]
        elif self._order == 1:
            if self._metabolites is True:
                interactions_of_interest = [interactionID for interactionID, genes in interaction_edges.items()
                                            if genes[0] in self._genes_of_interest
                                            or genes[1] in self._genes_of_interest]
            else:
                interactions_of_interest = [interactionID for interactionID, genes in interaction_edges.items()
                                            if (genes[0] in self._genes_of_interest
                                                or genes[1] in self._genes_of_interest) and genes[2] == 'p-p']
        self._interactions_of_interest = interactions_of_interest

        edge_list_df = (pd.DataFrame.from_dict(interaction_edges, orient='index',
                                               columns=['interactor1', 'interactor2', 'type'])
                        .merge(self._raw_info['sources'], how='left', left_index=True, right_on='interaction_id')
                        .query('interaction_id in @interactions_of_interest')
                        )
        return edge_list_df

    def build_network(self, edge_list_df, db_tidy_info):
        """Creates a network from a edge list DataFrame, adds node attributes."""
        self._network = nx.convert_matrix.from_pandas_edgelist(edge_list_df, source='interactor1', target='interactor2',
                                                               edge_attr=['is_experimental', 'interaction_id'])
        nx.set_node_attributes(self._network, db_tidy_info['proteins'])
        nx.set_node_attributes(self._network, db_tidy_info['short_names'])
        nx.set_node_attributes(self._network, db_tidy_info['localization'], name='localization')

        if self._order == 1:  # Label seed proteins in first-order networks.
            for node in self._network.nodes():
                if node in self._genes_of_interest:
                    self._network.node[node]['seed'] = 1
                else:
                    self._network.node[node]['seed'] = 0
        # Remove orphan nodes and self-loop edges.
        self._network.remove_edges_from(self._network.selfloop_edges())
        self._network.remove_nodes_from(list(nx.isolates(self._network)))

    def add_locus_tags(self):
        """Replace 'NA' strings in the short names attribute with their corresponding locus tag."""
        for node in self._network.nodes:
            if self._network.node[node]['short_name'] == 'NA':
                self._network.node[node]['short_name'] = node

    def make_network(self):
        """Generates a PPI network from a list of genes."""
        BioNetwork.query_db(self)
        tidy_db_info = BioNetwork.format_attribute_dictionaries(self)
        network_data = BioNetwork.make_edge_list(self)
        BioNetwork.build_network(self, network_data, tidy_db_info)
        BioNetwork.add_locus_tags(self)

    def write_gml(self, file_name):
        """Export the network as a GraphML file."""
        nx.write_graphml(self._network, '/'.join([os.getcwd(), file_name]))


class DENetwork(BioNetwork):
    """BioNetwork subclass with additional differential expression (DE)-related methods."""

    def __init__(self, gene_list_path, strain, order, detection_method, metabolites=False):
        super().__init__(gene_list_path, strain, order, detection_method, metabolites)
        self._genes_of_interest = helpers.get_genes(gene_list_path)
        self._de_info = DENetwork.process_de_genes_list(gene_list_path)

    @staticmethod
    def process_de_genes_list(gene_list_path):
        """Reads in a DESeq2 output gene list and returns a dictionary with differential expression information
        to use as node attributes. """
        de_genes = pd.read_csv(gene_list_path)
        de_genes = de_genes.rename(columns={de_genes.columns[0]: 'gene'})  # Rename first column (usually unnamed)
        de_info = dict()
        raw_de_info = de_genes[['gene', 'log2FoldChange', 'padj']].to_dict(orient='index')
        # Format dictionary to use as input for networkX node attributes
        for key, value in raw_de_info.items():
            de_info[value['gene']] = dict(log2FoldChange=value['log2FoldChange'],
                                          padj=value['padj'])
        return de_info

    def make_network(self):
        super().make_network()
        nx.set_node_attributes(self._network, self._de_info)


class CombinedNetwork(DENetwork):
    """DENetwork subclass with combined RNASeq (DE) and TnSeq information."""

    def __init__(self, de_genes_path, tnseq_genes_path, strain, order, detection_method, metabolites=False):
        super().__init__(de_genes_path, strain, order, detection_method, metabolites)
        self._de_genes = helpers.get_genes(de_genes_path)
        self._tnseq_genes = helpers.get_genes(tnseq_genes_path)
        self._genes_of_interest = list(set(self._de_genes).union(set(self._tnseq_genes)))

    def add_significance_source(self):
        """Adds a significance source indicating if it is from RNASeq or TnSeq"""
        for node in self._network.nodes():
            if (node in self._de_genes) and (node in self._tnseq_genes):
                self._network.node[node]['significance_source'] = 'both'
            elif node in self._de_genes:
                self._network.node[node]['significance_source'] = 'transcriptional'
            elif node in self._tnseq_genes:
                self._network.node[node]['significance_source'] = 'phenotypical'
            else:
                self._network.node[node]['significance_source'] = 'none'

    def make_network(self):
        super().make_network()
        CombinedNetwork.add_significance_source(self)

