import os.path

import pandas as pd


def ERCC():
    ercc = pd.read_table(os.path.dirname(__file__) + '/ERCC.tsv', index_col=1)

    return ercc
