from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
from os.path import join as pjoin
import pandas as pd
from indra.tools import assemble_corpus as ac
from process_abs import get_valid_genes
from indra.tools.gene_network import GeneNetwork
from indra import biopax


def read_phosphosite_owl(fname):
    bp = biopax.process_owl(fname)
    for stmt in bp.statements:
        for ev in stmt.evidence:
            ev.source_api = 'phosphosite'
            ev.epistemics = {'direct': True}
    return bp.statements


if __name__ == '__main__':

    if len(sys.argv) < 2:
        assemble_models = ['pysb', 'sif', 'cx']
    else:
        model_types = sys.argv[1:]
        if 'all' in model_types:
            assemble_models = ['pysb', 'sif', 'cx']
        else:
            assemble_models = sys.argv[1:]

    print('Assembling the following model types: %s' % \
          ', '.join(assemble_models))
    print('##############')

    outf = 'output/'
    reassemble = True

    def build_prior(genes, out_file):
        gn = GeneNetwork(genes, basename=pjoin(outf, 'pertbio'))
        stmts = gn.get_statements(filter=False)
        ac.dump_statements(stmts, out_file)
        return stmts

    if not reassemble:
        stmts = ac.load_statements(pjoin(outf, 'preassembled.pkl'))
    else:
        """
        data_genes = get_valid_genes()
        #prior_stmts = build_prior(data_genes, pjoin(outf, 'prior.pkl'))
        #prior_stmts = ac.map_grounding(prior_stmts,
        #                               save=pjoin(outf, 'gmapped_prior.pkl'))
        prior_stmts = ac.load_statements(pjoin(outf, 'prior.pkl'))
        reach_stmts = ac.load_statements(pjoin(outf, 'pertbio_stmts.pkl'))
        reach_stmts = ac.filter_no_hypothesis(reach_stmts)
        phosphosite_stmts = read_phosphosite_owl(
                    '../indra/models/phase3_eval/sources/Kinase_substrates.owl')
        reading_stmts = reach_stmts + phosphosite_stmts
        reading_stmts = ac.map_grounding(reading_stmts,
                                         save=pjoin(outf, 'gmapped_reading.pkl'))
        stmts = prior_stmts + reading_stmts

        stmts = ac.filter_grounded_only(stmts)
        stmts = ac.filter_genes_only(stmts, specific_only=False)
        stmts = ac.filter_human_only(stmts)
        stmts = ac.expand_families(stmts)
        #stmts = ac.filter_gene_list(stmts, data_genes, 'one')
        #stmts = ac.map_sequence(stmts, save=pjoin(outf, 'smapped.pkl'))
        """
        stmts = ac.load_statements(pjoin(outf, 'smapped.pkl'))
        stmts = ac.run_preassembly(stmts, return_toplevel=False,
                                   save=pjoin(outf, 'preassembled.pkl'),
                                   poolsize=4)

    ### PySB assembly
    if 'pysb' in assemble_models:
        pysb_model = assemble_pysb(stmts, data_genes,
                                   pjoin(outf, 'korkut_model_pysb.py'))
    ### SIF assembly
    if 'sif' in assemble_models:
        sif_str = assemble_sif(stmts, data, pjoin(outf, 'PKN-korkut_all_ab.sif'))
    ### CX assembly
    if 'cx' in assemble_models:
        cxa = assemble_cx(stmts, pjoin(outf, 'korkut_full_high_belief.cx'))


    # Plan - get ndex_neighborhood for each gene, plus indirect stmts
    # Get paths between all genes from PC
    
    # reach_stmts = ac.load_statements('pertbio_stmts.pkl')
    # Get statements for 
    # Load 
    # ac.load_statements(
