import os
from glob import glob

from tqdm import tqdm
import pandas as pd
import click
from sklearn.linear_model import LogisticRegression
import numpy as np

import readquant
from readquant.utils import BioMartQuery


def get_ERCC():
    return readquant.data.ERCC()['concentration in Mix 1 (attomoles/ul)']


def get_MT():
    ENSEMBL = BioMartQuery("hsapiens_gene_ensembl")
    ENSEMBL.add_filters(chromosome_name="MT")
    ENSEMBL.add_attributes("ensembl_gene_id")
    idx = pd.Index(ENSEMBL.stream()['ensembl_gene_id'])

    return idx


def get_rRNA():
    ENSEMBL = BioMartQuery("hsapiens_gene_ensembl")
    ENSEMBL.add_filters(biotype="rRNA")
    ENSEMBL.add_attributes("ensembl_gene_id")
    idx = pd.Index(ENSEMBL.stream()['ensembl_gene_id'])

    return idx


def get_detection_limit(ercc, quant, det_threshold=0.1):
    X = ercc[:, None]
    y = quant[ercc.index] >= det_threshold

    if y.sum() < 8:
        return np.inf

    lr = LogisticRegression(solver='liblinear', fit_intercept=True)
    lr.fit(X, y)
    midpoint = -lr.intercept_ / lr.coef_[0]

    return np.exp(midpoint[0])


def get_accuracy(ercc, quant, det_threshold=0.1):
    y = np.log(quant[ercc.index]) \
            .replace([np.inf, -np.inf], np.nan) \
            .dropna()

    if (y >= np.log(det_threshold)).sum() < 8:
        return -np.inf

    correlation = y.corr(ercc, method='pearson')

    return correlation


@click.command()
@click.argument('pattern', default='salmon/*_salmon_out')
@click.argument('output', default='sample_bio_qc.csv')
@click.option('--version', default='0.7.2')
def main(pattern='salmon/*_salmon_out', output='sample_bio_qc.csv', version=None):
    ercc = np.log(get_ERCC())
    MT = get_MT()
    rRNA = get_rRNA()

    print('Collected QC values')

    QCs = pd.DataFrame()
    for path in tqdm(glob(pattern)):
        quant = readquant.parse.read_salmon(path, version=version)

        qc_data = pd.Series()
        qc_data['detection_limit'] = get_detection_limit(ercc, quant)
        qc_data['accuracy'] = get_accuracy(ercc, quant)

        try:
            qc_data['ERCC_content'] = quant[ercc.index].sum()
            quant = quant.drop(ercc.index)
            quant = quant / quant.sum() * 1e6

        except ValueError:
            # ERCCs not present
            pass

        qc_data['num_genes'] = (quant > 1.).sum()
        qc_data['MT_content'] = quant[MT].sum()
        qc_data['rRNA_content'] = quant[rRNA].sum()

        QCs[path] = qc_data

    QCs.T.to_csv(output)


if __name__ == '__main__':
    main()
