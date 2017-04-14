from __future__ import print_function, unicode_literals
import re
import pandas
from indra.databases import hgnc_client
from indra.util import write_unicode_csv
from indra.literature import pubmed_client
import random

ab_file = 'MDACC_RPPA_Standard Ab List_Updated.xlsx'
ab_sheet = 'Current_Standard Ab List_304_'

ab_gene_map = {
        'ABL': ['ABL1', 'ABL2'],
        'AKT1,2,3': ['AKT1', 'AKT2', 'AKT3'],
        'H2AX': ['H2AFX'],
        'OCT4': ['POU5F1'],
        'RAB11A,B': ['RAB11A', 'RAB11B'],
        'RPS6K': ['RPS6KA1', 'RPS6KA2', 'RPS6KA3', 'RPS6KA4', 'RPS6KA5',
                  'RPS6KA6', 'RPS6KB1', 'RPS6KB2'],
        'TAU': ['MAPT'],
        'EMA': ['MUC1'],
        'GLUD': ['GLUD1', 'GLUD2'],
        'INSRB': ['INSR'],
        # https://www.cellsignal.com/products/primary-antibodies/hsp70-antibody/4872
        'HSP70': ['HSPA1L', 'HSPA1A', 'HSPA8'],
        'PKM2': ['PKM'],
        'C12ORF5': ['TIGAR'],
        'PDHK1': ['PDK1'],
        'XPF': ['ERCC4'],
        'CD26': ['DPP4'],
        'CNST43': ['GJA1'],
        'FRA1': ['FOSL1'],
        'HSP27': ['HSPB1'],
        # This appears to be the wrong/invalid gene name, the antibody is
        # listed as binding COX IV
        # https://www.cellsignal.com/products/primary-antibodies/cox-iv-3e11-rabbit-mab/4850
        'PTGS3': ['COX4I1', 'COX4I2'],
        'RIP': ['RIPK1'],
        'LC3AB': ['MAP1LC3A', 'MAP1LC3B', 'MAP1LC3C'],
        'RPA32': ['RPA2'],
        'HISTH3': ['HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E',
                   'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J'],
        }

# H3K9ME2 refers to Histone H3 dimethylated on lysine 9; an indicator of
# transcriptional silencing

# PAR (http://www.amsbio.com/productpage.aspx?code=4336-BPC-100)
# detects ribosylated proteins (not a specific gene/protein)

def read_ab_table(ab_file=ab_file, ab_sheet=ab_sheet):
    with open(ab_file, 'r') as fh:
        df = pandas.read_excel(ab_file, sheetname=ab_sheet, skiprows=5,
                               skipfooter=5, index_col=None)
    return df

def get_ab_genes(data):
    gene_list_raw = list(data['Gene Name'].values)
    gene_list = []
    for gn in gene_list_raw:
        mapped_genes = ab_gene_map.get(gn)
        if mapped_genes:
            gene_list += mapped_genes
        else:
            terms = re.split(',| ', gn)
            gene_list += [term for term in terms if term]
    gene_list = sorted(list(set(gene_list)))
    return gene_list

def check_ab_genes(genes):
    print('Invalid gene names\n------------------')
    valid_genes = []
    for gene in genes:
        hgnc_id = hgnc_client.get_hgnc_id(gene)
        if not hgnc_id:
            print(gene)
        else:
            valid_genes.append(gene)
    return valid_genes

def get_valid_genes():
    ab_data = read_ab_table()
    valid_genes = check_ab_genes(get_ab_genes(ab_data))
    return valid_genes

if __name__ == '__main__':
    valid_genes = get_valid_genes()
    write_unicode_csv('ab_gene_list.csv', [[gene] for gene in valid_genes])
    gene_pmids = {}
    all_pmids = set([])
    for gene in valid_genes:
        pmids = pubmed_client.get_ids_for_gene(gene)
        gene_pmids[gene] = pmids
        all_pmids = all_pmids.union(set(pmids))
    all_pmids = list(all_pmids)
    # The PMIDs went through a set so they are presumably not ordered, but
    # it never hurts to shuffle for good measure--this keeps the
    # full text/abstract distribution more uniform across the corpus
    random.Random(0).shuffle(all_pmids)
    write_unicode_csv('ab_pmid_list.csv', [[pmid] for pmid in all_pmids])

