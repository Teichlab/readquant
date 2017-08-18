import click

@click.command()
@click.argument("gencode", type=click.File('r'))
@click.argument("output", type=click.File('w'), default='-', required=False)

def main(gencode, output):
    """ in gencode, from release 25,
    gene and trancript ids on the chrY PAR regions have "_PAR_Y" appended
    this script simply remove those transcripts from the fasta file
    there are 74 transcripts like that in gencode.v26.pc_transcripts.fa
    usage: python filter_gencode.py gencode.v26.pc_transcripts.fa [filtered.fa]
    if [filtered.fa] is absent, will output to stdout.
    """
    for line in gencode:
        if line.startswith('>'):
            if '_PAR_Y' in line:
                to_write = False
            else:
                output.write(line)
                to_write = True
        else:
            if to_write:
                output.write(line)

if __name__ == "__main__":
    main()
