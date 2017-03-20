from __future__ import print_function

import sys
import struct
import gzip


import pandas as pd
import click
from xml.etree.ElementTree import Element, SubElement, tostring
import requests
import six
import pycurl
from urllib.parse import urlencode
from urllib.parse import quote_plus
from io import BytesIO


class PosModel:
    ''' Helper class for parsing positional bias parameters from Salmon.
    '''
    def __init__(self):
        self.models = {}

    def from_file(self, fn):
        f = gzip.open(fn)
        b = f.read()

        offset = 0
        uint_size = struct.calcsize('I')
        double_size = struct.calcsize('d')

        num_models = struct.unpack_from('I', b[offset:])[0]
        offset += uint_size

        length_bins = []
        for i in range(0, num_models):
            length_bins.append(struct.unpack_from('I', b[offset:])[0])
            offset += uint_size

        models = []
        for i in range(0, num_models):
            model_bins = struct.unpack_from('I', b[offset:])[0]
            offset += uint_size
            model = list(struct.unpack_from('d'*model_bins, b[offset:]))
            offset += model_bins * double_size
            models.append(model)
        f.close()
        self.models = dict(zip(length_bins, models))



class BioMartQuery(object):

    """ Wrapper for making BioMart requests

    Will build and execute an xml request to BioMart using either pycurl or
    requests.


    Example: All human mitochondrial genes

    >>> q = BioMartQuery("hsapiens_gene_ensembl")
    >>> q.add_filters(chromosome_name="MT")
    >>> q.add_attributes("ensembl_gene_id")
    >>> df = q.stream()

    """

    biomart_base = 'http://www.ensembl.org/biomart/martservice?query='

    def __init__(self, dataset, method='pycurl'):
        """

        Parameters
        ----------
        dataset, str
            The BioMart organism dataset to query. E.g. 'hsapiens_gene_ensembl'
        method, str, default 'pycurl'
            Python library to use for the request.

        """
        self._query = Element('Query', attrib={
            'virtualSchemaName': 'default',
            'formatter': 'TSV',
            'header': '0',
            'uniqueRows': '0',
            'count': '',
            'datasetConfigVersion': '0.6'
        })
        self._dataset_query = SubElement(self._query, 'Dataset', attrib={
            'name': dataset,
            'interface': 'default'
        })
        self.headings = []
        if method == 'pycurl':
            self.getter = self._pycurl_piecewise
        else:
            self.getter = self._requests_piecewise

    def add_attributes(self, *attributes):
        """

        Parameters
        ----------
        *attributes, str
            Attributes to be queried. E.g. 'ensembl_gene_id'

        """
        for attr in attributes:
            self.headings.append(attr)
            self._dataset_query.append(
                Element('Attribute', attrib={'name': attr}))

    def add_filters(self, **filters):
        """

        Parameters
        ----------
        **filters, str
            Key-value pairs of filters. E.g. chromosome_name='MT'

        """

        for key, value in six.iteritems(filters):
            self._dataset_query.append(Element(
                'Filter', attrib={'name': key, 'value': value}))

    def _requests_piecewise(self, handle):

        request_url = self.biomart_base + self.encoded_request()
        print("Making request to biomart...")
        response = requests.get(self.biomart_base + str(self), stream=True)
        if not response.ok:
            print("Request failed: %d" % response.status_code)
            return
        print("Response okay. Downloading...")
        total_length = response.headers.get('content-length')
        total_downloaded = 0
        if total_length is None:
            for block in response.iter_content(1024):
                handle.write(block)
                total_downloaded += 1024
                # stdscr.addstr('%d MB Downloaded' %
                #                  (total_downloaded / 1e6))
                # stdscr.refresh()
        else:
            with click.progressbar(length=int(total_length),
                                   label="Downloading reference") as bar:
                for block in response.iter_content(1024):
                    handle.write(block)
                    bar.update(1024)

    def _pycurl_piecewise(self, handle):

        request_url = self.biomart_base + self.encoded_request()
        c = pycurl.Curl()
        c.setopt(c.URL, request_url)
        c.setopt(c.WRITEDATA, handle)
        c.perform()
        c.close()

    def download(self, outfile):
        """ Download request to file

        Parameters
        ----------
        outfile, str
            Location of output file

        """

        with open(outfile, 'wb') as outbuffer:
            self.getter(outbuffer)

    def stream(self):
        """ Stream request to pandas DataFrame

        Returns
        -------
        df, pd.DataFrame
            DataFrame of returned values with headers

        """

        outbuffer = BytesIO()
        self.getter(outbuffer)
        outbuffer.seek(0)

        df = pd.read_csv(outbuffer, sep='\t', encoding='utf-8')
        df.set_axis(1, self.headings)

        return df

    def encoded_request(self):
        result = '<?xml version="1.0" encoding="UTF-8" ?><!DOCTYPE Query>'
        result += tostring(self._query, 'utf-8').decode('utf-8')

        return quote_plus(result)

    def __str__(self):
        result = '<?xml version="1.0" encoding="UTF-8" ?><!DOCTYPE Query>'
        result += tostring(self._query, 'utf-8').decode('utf-8')

        return result
