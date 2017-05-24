import os

import click

import readquant

@click.command()
@click.argument('name', default='hsapiens_gene_ensembl')
def main(name='hsapiens_gene_ensembl'):
    ''' Get files from Ensembl needed to make gene expression references.

    Make a new folder with the data you ran the script for provenence, and run it there.

    By defualt downloads the human reference files (hsapiens_gene_ensembl).
    The name of the Ensembl database
    is given as a an optional argument. For mouse, use mmusculus_gene_ensembl.
    '''
    biomart_url = 'http://www.ensembl.org/biomart/martservice'
    tx_temp, ga_temp, gm_temp = readquant.data.reference_templates()

    tx_temp = tx_temp.replace('hsapiens_gene_ensembl', name)
    ga_temp = ga_temp.replace('hsapiens_gene_ensembl', name)
    gm_temp = gm_temp.replace('hsapiens_gene_ensembl', name)

    with open('cDNA.query.xml', 'w') as fh:
        fh.write(tx_temp)

    with open('gene_annotation.query.xml', 'w') as fh:
        fh.write(ga_temp)

    with open('genemap.query.xml', 'w') as fh:
        fh.write(gm_temp)

    os.system('curl --data-urlencode query@cDNA.query.xml {}?query= -o cDNA.fasta'.format(biomart_url))
    os.system('curl --data-urlencode query@gene_annotation.query.xml {}?query= -o gene_annotation.csv'.format(biomart_url))
    os.system('curl --data-urlencode query@genemap.query.xml {}?query= -o genemap.tsv'.format(biomart_url))

if __name__ == '__main__':
    main()
