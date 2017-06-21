from __future__ import print_function
from glob import iglob

import os
import numpy as np
import pandas as pd
from tqdm import tqdm

def read_kallisto(sample_path):
    ''' Function for reading a Kallisto quantification result.

    Returns
    -------
    A pandas.Series with the expression values in the sample.
    '''
    quant_file = sample_path + '/abundance.tsv'
    df = pd.read_table(quant_file, engine='c',
                                   usecols=['target_id', 'tpm'],
                                   index_col=0,
                                   dtype={'target_id': np.str, 'tpm': np.float64})

    df = df.rename(columns={'tpm': 'TPM'})
    return df['TPM']


def read_salmon(sample_path, isoforms=False, version='0.7.2', unit='TPM'):
    ''' Function for reading a Salmon quantification result.

    Parameters
    ----------
    isoforms : bool, default False
        Whether to parse isoform level expression or gene level expression.

    version : str, default '0.7.2'
        The version of Salmon which generated the directory. Currently
        supports '0.7.2', '0.6.0' and '0.4.0'. (Other versions might be compatible
        with these.)

    Returns
    -------
    A pandas.Series with the expression values in the sample.
    '''
    if isoforms:
        quant_file = sample_path + '/quant.sf'
    else:
        quant_file = sample_path + '/quant.genes.sf'

    read_kwargs = {
        '0.7.2': {
            'engine': 'c',
            'usecols': ['Name', unit],
            'index_col': 0,
            'dtype': {'Name': np.str, unit: np.float64}
        },
        '0.6.0': {
            'engine': 'c',
            'usecols': ['Name', unit],
            'index_col': 0,
            'dtype': {'Name': np.str, unit: np.float64}
        },
        '0.4.0': {
            'engine': 'c',
            'comment': '#',
            'header': None,
            'names': ['Name', 'length', 'TPM', 'NumReads'],
            'usecols': ['Name', 'TPM'],
            'index_col': 0,
            'dtype': {'Name': np.str, 'TPM': np.float64}
        }
    }

    if not os.path.isfile(quant_file):
        print("WARNING: Could not find file: %s" % quant_file)
        return
    else:
        df = pd.read_table(quant_file, **read_kwargs[version])

        df = df.rename(columns={'Name': 'target_id'})
        return df[unit]


def read_cufflinks(sample_path, isoforms=False):
    ''' Function for reading a Cufflinks quantification result.

    Returns
    -------
    A pandas.Series with the expression values in the sample.
    '''
    if isoforms:
        quant_file = sample_path + '/isoforms.fpkm_tracking'
    else:
        quant_file = sample_path + '/genes.fpkm_tracking'
    df = pd.read_table(quant_file, engine='c',
                                   usecols=['tracking_id', 'FPKM'],
                                   index_col=0,
                                   dtype={'tracking_id': np.str, 'FPKM': np.float64})

    df['tracking_id'] = df.index
    df = df.groupby('tracking_id').sum()
    df['TPM'] = df['FPKM'] / df['FPKM'].sum() * 1e6

    df = df.rename(columns={'tracking_id': 'target_id'})
    return df['TPM']


def read_quants(pattern='salmon/*_salmon_out', tool='salmon', **kwargs):
    ''' Read quantification results from every directory matching the glob
    in pattern.

    Parameters
    ----------
    tool, str, default 'salmon'
        The quantification tool used to generate the results. Currently
        supports 'salmon', 'sailfish', 'kallisto', and 'cufflinks'.

    **kwargs,
        kwargs are passed on to the tool specific sample parser. See documentation
        for individual parsers for details.

    Returns
    -------
    A pandas.DataFrame where columns are samples, rows are genes, and cells
    contain the expression value.
    '''
    sample_readers = {
        'salmon': read_salmon,
        'sailfish': read_salmon,
        'kallisto': read_kallisto,
        'cufflinks': read_cufflinks
    }

    quant_reader = sample_readers[tool]

    quants = pd.DataFrame()
    for sample_path in tqdm(iglob(pattern)):
        sample_quant = quant_reader(sample_path, **kwargs)
        if sample_quant is not None:
            quants[sample_path] = sample_quant

    return quants


def read_salmon_3p_bias(pattern='salmon/*_salmon_out/'):
    ''' Read a smoothed representation of 3p bias for each sample.
    '''
    from .utils import PosModel
    p = PosModel()

    from sklearn import gaussian_process
    kernel = gaussian_process.kernels.RBF(50., (50., 50.))
    gpr = gaussian_process.GaussianProcessRegressor(kernel=kernel, normalize_y=False)
    
    xx = np.linspace(0, 100, 20)[:, None]
    xxx = np.linspace(0, 100, 100)[:, None]

    pattern += 'aux_info/obs3_pos.gz'

    sample_3p_bias = pd.DataFrame()
    for f in tqdm(iglob(pattern)):
        p.from_file(f)
        yy = pd.DataFrame(p.models).iloc[:, 2]
        yyy = gpr.fit(xx, yy).predict(xxx)
        sample_3p_bias[f.replace('aux_info/obs3_pos.gz', '')] = yyy

    return sample_3p_bias


def read_salmon_qc(sample_path, flen_lim=(100, 100), version='0.7.2'):
    ''' Parse technical quality control data from a Salmon quantification
    result.

    Parameters
    ----------
    flen_lim, tuple (int start, int end), default (100, 100)
        How many bases to remove from start and end of fragment length
        distribution when calculating the robust mode. This is too see if
        things roughly worked out even if the max FLD Salmon parameter was
        set too small.

    version, str, default '0.6.0'
        The version of Salmon which generated the directory. Currently
        supports '0.6.0' and '0.4.0'. (Other versions might be compatible
        with these.)

    Returns
    -------
    A pandas.Series with technical information from the Salmon results for
    the sample.
    '''
    try:
        flen_dist = np.fromfile(sample_path + '/libParams/flenDist.txt', sep='\t')
        global_fl_mode = flen_dist.argmax()
        robust_fl_mode = flen_dist[flen_lim[0]:-flen_lim[1]].argmax() + flen_lim[0]
    except FileNotFoundError:
        global_fl_mode = 0
        robust_fl_mode = 0

    if version == '0.7.2':
        qc_data = pd.read_json(sample_path + '/aux_info/meta_info.json', typ='series')
        qc_data = qc_data[['num_processed', 'num_mapped', 'percent_mapped']]
        qc_data['global_fl_mode'] = global_fl_mode
        qc_data['robust_fl_mode'] = robust_fl_mode

    if version == '0.6.0':
        qc_data = pd.read_json(sample_path + '/aux/meta_info.json', typ='series')
        qc_data = qc_data[['num_processed', 'num_mapped', 'percent_mapped']]
        qc_data['global_fl_mode'] = global_fl_mode
        qc_data['robust_fl_mode'] = robust_fl_mode

    if version == '0.4.0':
        qc_data = pd.Series()
        log_file = sample_path + '/logs/salmon_quant.log'
        with open(log_file) as fh:
            for l in fh:
                if 'Observed ' in l:
                    frags = int(l.split('Observed ')[-1].split(' total')[0])
                    qc_data['num_processed'] = frags

                if 'mapping rate' in l:
                    rate = float(l.split(' = ')[1].split('%')[0])
                    qc_data['percent_mapped'] = rate


        qc_data['global_fl_mode'] = global_fl_mode
        qc_data['robust_fl_mode'] = robust_fl_mode

    return qc_data


def read_tophat_qc(sample_path):
    ''' Parse technical quality control data from TopHat alignment results.

    Parameters
    ----------
    sample_path, str
        The path to the resulting TopHat directory.

    Returns
    -------
    A pandas.Series with technical information from the TopHat alignment results
    for the sample.
    '''
    with open(sample_path + '/align_summary.txt') as fh:
        for l in fh:
            if 'Input' in l:
                n_reads = int(l.split()[-1])
                break

        for l in fh:
            if 'overall' in l:
                pct_mapped = float(l.split('%')[0])
                break

    qc_data = pd.Series({'input_reads': n_reads,
                         'pct_mapped': pct_mapped})

    return qc_data


def read_qcs(pattern='salmon/*_salmon_out', tool='salmon', **kwargs):
    ''' Read technical quality control data results from every directory
    matching the glob in pattern.

    Parameters
    ----------
    tool, str, default 'salmon'
        The quantification tool used to generate the results. Currently
        supports 'salmon' and 'sailfish'.

    **kwargs,
        kwargs are passed on to the tool specific sample parser. See documentation
        for individual parsers for details.

    Returns
    -------
    A pandas.DataFrame where rows are samples, and columns are technical
    features extrated from the tool results.
    '''
    sample_readers = {
        'salmon': read_salmon_qc,
        'sailfish': read_salmon_qc,
        'tophat': read_tophat_qc
    }

    qc_reader = sample_readers[tool]

    QCs = pd.DataFrame()
    for sample_path in tqdm(iglob(pattern)):
        try:
            sample_qc = qc_reader(sample_path, **kwargs)

        except ValueError:
            print('Error parsing {}'.format(sample_path))
            raise

        QCs[sample_path] = sample_qc

    return QCs.T
