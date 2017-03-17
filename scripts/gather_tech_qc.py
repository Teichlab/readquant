import click

import readquant


@click.command()
@click.argument('pattern', default='salmon/*_salmon_out')
@click.argument('output', default='sample_qc.csv')
@click.option('--version', default='0.7.2')
def main(pattern='salmon/*_salmon_out', output='sample_qc.csv', version=None):
    QCs = readquant.read_qcs(pattern=pattern, tool='salmon', version=version)
    QCs.to_csv(output)


if __name__ == '__main__':
    main()
