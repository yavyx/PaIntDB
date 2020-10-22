import pandas as pd
import csv
import os
import pickle5 as pickle

from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS


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
    # Load GO term association dictionary
    with open(os.path.join('data', 'go_association.pickle'), 'rb') as handle:
        go_association = pickle.load(handle)

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


