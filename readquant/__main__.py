import click
import readquant


@click.group()
def cli():
    pass


_common_args = [
   click.argument('pattern'),
   click.argument('output'),
   click.option('--version', default="0.7.2")
]


def pass_common_args(func):
    for command in _common_args:
        func = command(func)
    return func


@cli.command()
@pass_common_args
def qc(pattern='salmon/*_salmon_out', output='sample_qc.csv', version=None):
    QCs = readquant.read_qcs(pattern=pattern, tool='salmon', version=version)
    QCs.to_csv(output)


@cli.command()
@click.option('--tool', default='salmon')
@pass_common_args
def quant(pattern, output, version="0.7.2", tool='salmon'):
    quants = readquant.read_quants(pattern=pattern, tool=tool, version=version)
    quants.to_csv(output)


if __name__ == "__main__":
    cli()
