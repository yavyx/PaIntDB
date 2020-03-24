import pandas as pd
import csv
import os

from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.godag_plot import plot_gos, plot_results

import bio_networks.network_generator as ng


def make_go_association_dict(path):
    """Returns a GO association dictionary that maps genes to GO terms can be used by GOAtools for GO enrichment."""
    pseudomonas_go_associations = dict()
    go_domains = {'BP': 'biological_process',
                  'CC': 'cellular_component',
                  'MF': 'molecular_function'}

    go_associations = pd.read_csv(path)
    #locus_tag, go_id, namespace = [list(go_associations[col]) for col in ['Locus Tag', 'Accession', 'Namespace']]
    #print(locus_tag)

    for domain_short, domain in go_domains.items():
        with open(path) as go_associations:
            go_associations_reader = csv.reader(go_associations)
            pseudomonas_go_associations[domain_short] = dict()
            for row in go_associations_reader:
                locus_tag, go_id, namespace = row[0], row[4], row[6]
                if domain == namespace:
                    if locus_tag not in pseudomonas_go_associations[domain_short]:
                        pseudomonas_go_associations[domain_short][locus_tag] = set()
                        pseudomonas_go_associations[domain_short][locus_tag].add(go_id)
                    else:
                        pseudomonas_go_associations[domain_short][locus_tag].add(go_id)
    return pseudomonas_go_associations


def get_genes(path):
    """Returns a list of genes from a DE results table"""
    with open(path) as gene_list:
        gene_list = csv.reader(gene_list)
        gene_list = [row[0] for row in gene_list if row[0].startswith('P')]
    return gene_list


def get_enrichment_results(enrichment_results):
    """Returns a dataFrame of the enrichment results."""
    fields = ['id', 'name', 'namespace', 'enrichment', 'ratio_in_study', 'ratio_in_pop', 'p_fdr_bh', 'study_count',
              'study_items']
    enrichment_results = [go_term.get_field_values(fldnames=fields) for go_term in enrichment_results]
    enrichment_results_df = pd.DataFrame.from_records(enrichment_results, columns=fields)
    # enrichment_results_df.study_items = enrichment_results_df.study_items.str.split(", ")
    return enrichment_results_df


def get_gene_ratios(enrichment_results):
    """Calculate gene ratios for each GO term and add them to the results DataFrame."""
    genes_per_go_term = dict()

    with open('/home/javier/Documents/Evelyn/pseudomonas_GO_gene_sets.tsv') as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        for row in reader:
            genes_per_go_term[row[1]] = len(row) - 2

    genes_per_go_term = pd.DataFrame.from_dict(genes_per_go_term, orient='index', columns=['n_total_genes'])
    enrichment_results = pd.merge(enrichment_results, genes_per_go_term, how='left', left_on='id', right_index=True)
    return enrichment_results


def go_enrichment_stats(go_results, tigs=None, pigs=None):
    """Returns dataFrame with ratios of genes assigned to significant GO terms. """
    results = go_results
    if tigs:
        n_tigs = [len(set(tigs) & set(results.study_items[row])) for row in range(len(results))]
        results['n_tigs'] = n_tigs
        tigs_ratio = [n_tigs[row] / len(results.study_items[row])
                      for row in range(len(results))]
        results['tigs_ratio'] = tigs_ratio
    if pigs:
        n_pigs = [len(set(pigs) & set(results.study_items[row])) for row in range(len(results))]
        results['n_pigs'] = n_pigs
        pigs_ratio = [n_pigs[row] / len(results.study_items[row])
                      for row in range(len(results))]
        results['pigs_ratio'] = pigs_ratio
    return results


def map_pa14_genes(gene_list):
    """Takes a list of PA14 genes and returns the corresponding PAO1 names."""
    pao1_pa14_mapping = dict()
    mapping_path = os.path.join(os.getcwd(), 'data', 'ortholuge_pa14_to_pao1_20190708.tsv')
    with open(mapping_path) as mapping:
        reader = csv.reader(mapping, delimiter='\t')
        for row in reader:
            pao1_pa14_mapping[row[10]] = row[4]

    pao1_genes = [pao1_pa14_mapping[gene] for gene in gene_list if gene in pao1_pa14_mapping.keys()]
    return pao1_genes


def map_pao1_genes(gene_list):
    """Takes a list of PAO1 genes and returns the corresponding PA14 names."""
    pa14_pao1_mapping = dict()
    mapping_path = os.path.join(os.getcwd(), 'data', 'ortholuge_pa14_to_pao1_20190708.tsv')
    with open(mapping_path) as mapping:
        reader = csv.reader(mapping, delimiter='\t')
        for row in reader:
            pa14_pao1_mapping[row[4]] = row[10]

    pa14_genes = [pa14_pao1_mapping[gene] for gene in gene_list if gene in pa14_pao1_mapping.keys()]
    return pa14_genes


def run_go_enrichment(strain, genes_of_interest, significant=True, cutoff=0.05,
                      use_parent_terms=True):
    go_association = make_go_association_dict(os.path.join('data', 'PAO1_gene_ontology.csv'))
    background_genes = get_genes(os.path.join('data', strain + '_all_genes.csv'))
    obo_go_fname = download_go_basic_obo()
    obo_dag = GODag('go-basic.obo')
    # genes_of_interest = list(network.nodes)

    if strain == 'PA14':
        genes_of_interest = map_pa14_genes(genes_of_interest)
        background_genes = map_pa14_genes(background_genes)

    goea_obj = GOEnrichmentStudyNS(background_genes, go_association, obo_dag,
                                   propagate_counts=use_parent_terms,
                                   alpha=cutoff,
                                   methods=['fdr_bh'])
    goea_results = goea_obj.run_study(genes_of_interest)

    if significant is True:
        goea_results = [result for result in goea_results if result.p_fdr_bh < cutoff]

    enrichment_results = get_enrichment_results(goea_results)
    return [enrichment_results, goea_results]


def network_from_go_term(go_term_of_interest, enrichment_results, combined=False, rnaseq_genes=None, tnseq_genes=None,
                         order=0):
    """Takes GO enrichment results and returns the list of genes of interest of the selected GO term."""
    db_path = '/home/javier/Documents/PaIntDB2/PaIntDB.db'
    genes_in_go_term = enrichment_results.loc[enrichment_results.id == go_term_of_interest, 'study_items'].tolist()
    genes_in_go_term = genes_in_go_term[0]

    if combined:
        network = ng.make_network(genes_in_go_term, db_path, strain='PAO1', detection_method='all', order=order)
        network = ng.add_significance_source(network, rnaseq_genes, tnseq_genes)
    else:
        gene_list = rnaseq_genes
        network = ng.make_network(genes_in_go_term, db_path, strain='PAO1', detection_method='all', rna_seq=True,
                                  order=order)

    return network

