#!/usr/bin/env python3


def main():
    import argparse
    import os
    from extra._version import __version__
    import subprocess
    import sys

    parser = argparse.ArgumentParser(description='Divide read sets into subsets',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -i File -r File -s File [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', type=argparse.FileType(), metavar='File',
                         help='Csv file with reads merged by label')
    io_args.add_argument('-r', '--reads', type=argparse.FileType('w'), metavar='File',
                         help='Output csv file with information about reads')
    io_args.add_argument('-s', '--summary', type=argparse.FileType('w'), metavar='File',
                         help='Output csv file with subsets summary')

    div_args = parser.add_argument_group('Division arguments')
    div_args.add_argument('-S', '--significance', type=float, metavar='Float', default=.001,
                          help='Division significance (default: %(default)s)')
    div_args.add_argument('-p', '--pairs', type=int, metavar='Int', default=5,
                          help='Maximum number of evaluated pairs of errors (default: %(default)s)')
    div_args.add_argument('-c', '--coverage', type=int, metavar='Int', default=50,
                          help='Minimum subset size (default: %(default)s)')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    args.reads.write('# %s\n' % ' '.join(sys.argv))
    args.reads.write('read\tsegment\terrors\tlabel\tsubset\n')
    args.reads.flush()
    args.summary.write('# %s\n' % ' '.join(sys.argv))
    args.summary.write('subset\tcoverage\tlabels\tconsensus\n')
    args.summary.flush()
    current_dir = os.path.dirname(__file__)

    returncode = subprocess.run(['%s/bin/divide_reads' % current_dir,
                                 str(args.significance),
                                 str(args.pairs),
                                 str(args.coverage)],
                                stdin=args.input,
                                stdout=args.reads,
                                stderr=args.summary).returncode
    if returncode != 0:
        sys.stderr.write('Error while dividing reads\n')
        exit(1)


if __name__ == '__main__':
    main()
