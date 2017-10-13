#!/usr/bin/env python3


def main():
    import argparse
    import os
    from extra._version import __version__
    import subprocess
    import sys

    parser = argparse.ArgumentParser(description='Clip sequences and filter by distance to target',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -s File -S File -o File [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-s', '--sequences', help='Input file, every line - single sequence',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-S', '--segments', help='Input fasta with segments (cannot be stdin "-")',
                         required=True, metavar='File')
    io_args.add_argument('-o', '--output', help='Output csv file with filtered and clipped sequences',
                         type=argparse.FileType('w'), required=True, metavar='File')

    cl_args = parser.add_argument_group('Clipping arguments')
    cl_args.add_argument('--no-clipping', help='Do not clip segments', action='store_true')
    cl_args.add_argument('-r', '--range', help='Positions range (default: %(default)s)',
                         metavar=('Int', 'Int'), type=int, nargs=2, default=[40, 290])

    filt_args = parser.add_argument_group('Filtering arguments')
    filt_args.add_argument('-d', '--min-dist', type=int, metavar='Int', default=2,
                           help='Minimal Hamming distance to clipped segments (default: %(default)s)')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    args.output.write('# %s\n' % ' '.join(sys.argv))
    args.output.flush()
    current_dir = os.path.dirname(__file__)
    returncode = subprocess.run(['%s/bin/clip_and_filter' % current_dir,
                                 args.segments,
                                 str(args.min_dist),
                                 str(args.range[0]),
                                 str(args.range[1]),
                                 str(args.no_clipping)],
                                stdin=args.sequences,
                                stdout=args.output).returncode
    if returncode != 0:
        sys.stderr.write('Clipping and filtering: error\n')
        exit(1)


if __name__ == '__main__':
    main()
