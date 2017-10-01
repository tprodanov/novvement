#!/usr/bin/env python3


def main():
    import argparse
    import os
    from extra._version import __version__
    import subprocess

    parser = argparse.ArgumentParser(description='Construct Hamming graph on CDRs and label reads by components',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -i File -o File [--tau Int]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', help='Input CDRs (fasta file)',
                         required=True, metavar='File')
    io_args.add_argument('-o', '--output', help='Output csv file',
                         required=True, metavar='File')

    hamming_args = parser.add_argument_group('Hamming graph arguments')
    hamming_args.add_argument('-t', '--tau', metavar='Int', type=int, default=3,
                         help='Edges constructed by Hamming distance --tau (default: %(default)s)')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    current_dir = os.path.dirname(__file__)
    returncode = subprocess.run(['%s/bin/hamming' % current_dir,
                                 args.input,
                                 args.output,
                                 str(args.tau)]).returncode
    if returncode != 0:
        import sys
        sys.stderr.write('Error in hamming graph construction\n')
        exit(1)


if __name__ == '__main__':
    main()
