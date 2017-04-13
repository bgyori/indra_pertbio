from __future__ import print_function, unicode_literals
import re
import pandas
from indra.databases import hgnc_client

ab_file = 'MDACC_RPPA_Standard Ab List_Updated.xlsx'
ab_sheet = 'Current_Standard Ab List_304_'

ab_gene_map = {
        'ABL': ['ABL1', 'ABL2'],
        'AKT1,2,3': ['AKT1', 'AKT2', 'AKT3'],
        'H2AX': ['H2AFX'],
        'OCT4': ['POU5F1'],
        'RAB11A,B': ['RAB11A', 'RAB11B'],
        'RPS6K': ['RPS6KA1', 'RPS6KA2', 'RPS6KA3', 'RPS6KA4', 'RPS6KA5',
                  'RPS6KA6', 'RPS6KB1', 'RPS6KB2']
        '''
        C12ORF5
        CD26
        CNST43
        EMA
        FRA1
        GLUD
        H3K9ME2
        HISTH3
        HSP27
        HSP70
        INSRB
        LC3AB
        PAR
        PDHK1
        PKM2
        PTGS3
        RIP
        RPA32
        TAU
        XPF
        '''
        }


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
    for gene in genes:
        hgnc_id = hgnc_client.get_hgnc_id(gene)
        if not hgnc_id:
            print(gene)

