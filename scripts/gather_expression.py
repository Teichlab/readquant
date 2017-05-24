import click

import readquant


@click.command()
@click.argument('pattern', default='salmon/*_salmon_out')
@click.argument('output', default='expression.csv')
@click.option('--version', default='0.7.2')
@click.option('--unit', default='NumReads')
@click.option('--isoforms', default=0)
def main(pattern='salmon/*_salmon_out', output='expression.csv', unit='NumReads', version=None, isoforms=0):
    expr = readquant.read_quants(pattern=pattern, version=version, unit=unit, isoforms=isoforms)
    expr.to_csv(output)


if __name__ == '__main__':
    main()
