import os.path

import pandas as pd


def ERCC():
    ercc = pd.read_table(os.path.dirname(__file__) + '/ERCC.tsv', index_col=1)

    return ercc

def reference_templates():
    ''' Get XML templates to query Biomart with.
    '''
    with open(os.path.dirname(__file__) + '/template_transcriptome.xml') as fh:
        tx_template = fh.read()

    with open(os.path.dirname(__file__) + '/template_gene_annotation.xml') as fh:
        ga_template = fh.read()

    with open(os.path.dirname(__file__) + '/template_genemap.xml') as fh:
        gm_template = fh.read()

    return tx_template, ga_template, gm_template
