import csv


def remove_nones(dictionary):
    """Replace None values with 'NA' strings, because NetworkX's GraphML does not accept None values as attributes."""
    for key, value in dictionary.items():
        if value is None:
            dictionary[key] = 'NA'
    return dictionary


def get_genes(path):
    """Returns a list of genes from common results tables. Assumes the genes are in the first column."""
    with open(path) as full_gene_list:
        full_gene_list = csv.reader(full_gene_list)
        gene_list = [row[0] for row in full_gene_list]
        del gene_list[0]  # Remove first "gene" (header)
    return gene_list
