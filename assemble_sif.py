from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.assemblers import SifAssembler
import indra.tools.assemble_corpus as ac
from indra.statements import *


def assemble_sif(stmts, out_file):
    """Return an assembled SIF."""
    # Filter for high-belief statements
    stmts = ac.filter_belief(stmts, 0.95)
    stmts = ac.filter_top_level(stmts)
    stmts = ac.filter_direct(stmts)
    print(len(stmts))
    # Make the SIF model
    sa = SifAssembler(stmts)
    sa.make_model(use_name_as_key=True, include_mods=True,
                  include_complexes=True)
    sif_str = sa.print_model()
    # assemble and dump a cx of the sif
    with open(out_file, 'wb') as fh:
        fh.write(sif_str.encode('utf-8'))
