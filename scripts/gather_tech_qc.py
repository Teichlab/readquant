import click

import readquant


@click.command()
@click.argument('pattern')
@click.argument('output')
@click.option('--version')
def main(pattern='salmon/*_salmon_out', output='sample_qc.csv', version=None):
    QCs = readquant.read_qcs(pattern=pattern, tool='salmon', version=version)
    QCs.to_csv(output)


if __name__ == '__main__':
    main()
