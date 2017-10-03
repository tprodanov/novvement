#!/usr/bin/env python3


def main():
    import argparse
    import os
    from extra._version import __version__
    import subprocess
    import sys

    parser = argparse.ArgumentParser(description='Change read labels by Hamming graph components',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -i File -o File [-t Int] [-d]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', help='Input csv file (<read name> <label>, first line ignored)',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-o', '--output', help='Output csv file',
                         type=argparse.FileType('w'), required=True, metavar='File')

    hamming_args = parser.add_argument_group('Hamming graph arguments')
    hamming_args.add_argument('-t', '--tau', metavar='Int', type=int, default=3,
                              help='Edges constructed by Hamming distance --tau (default: %(default)s)')
    hamming_args.add_argument('-d', '--accept-diff-lengths', action='store_true',
                              help='Possible edges between labels with different lengths')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    args.output.write('# %s\n' % ' '.join(sys.argv))
    args.output.flush()
    current_dir = os.path.dirname(__file__)
    returncode = subprocess.run(['%s/bin/hamming' % current_dir,
                                 str(args.tau),
                                 str(args.accept_diff_lengths)],
                                stdin=args.input,
                                stdout=args.output).returncode
    if returncode != 0:
        sys.stderr.write('Error in hamming graph construction\n')
        exit(1)


if __name__ == '__main__':
    main()
