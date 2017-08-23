import sys
import collections
from io import BufferedReader
import gzip

import click

import pandas as pd


def stream_fastq(file_handler):
    ''' Generator which gives all four lines if a fastq read as one string
    '''
    next_element = ''
    for i, line in enumerate(file_handler):
        next_element += line
        if i % 4 == 3:
            yield next_element
            next_element = ''


def read_fastq(filename):
    """
    return a stream of FASTQ entries, handling gzipped and empty files
    """
    if not filename:
        return itertools.cycle((None,))
    if filename.endswith('gz'):
        filename_fh = BufferedReader(gzip.open(filename, mode='rt'))
    else:
        filename_fh = open(filename)
    return stream_fastq(filename_fh)


def write_fastq(filename):
    """
    return a handle for FASTQ writing, handling gzipped files
    """
    if filename:
        if filename.endswith('gz'):
            filename_fh = gzip.open(filename, mode='wb')
        else:
            filename_fh = open(filename, mode='w')
    else:
        filename_fh = None
    return filename_fh

@click.command()
@click.argument('fq1')
@click.argument('fq2')
@click.argument('adapter')
def concata_filter(fq1, fq2, adapter):
    copy_hist = collections.defaultdict(int)

    r1 = read_fastq(fq1)
    r2 = read_fastq(fq2)

    for read in zip(r1, r2):

        copy_hist[read[0].count(adapter) + read[1].count(adapter)] += 1

    df = pd.DataFrame({'Fragments': copy_hist})
    df.index.name = 'Copies'
    df.to_csv(sys.stdout)

if __name__ == '__main__':
    concata_filter()
