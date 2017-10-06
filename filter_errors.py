#!/usr/bin/env python3


def load_segment_coverage(f):
    coverages = {}
    for line in f:
        segment, coverage = line.strip().split()
        coverages[segment] = int(coverage)
    return coverages


def filter_errors(inp, outp, pos_range, coverages, cov_threshold, keep_n):
    import sys

    l, r = pos_range
    outp.write('# %s\n' % ' '.join(sys.argv))

    for line in inp:
        outp.write(line)
        if not line.startswith('#'):
            break
    
    for line in inp:
        read, segment, position, alt = line.strip().split('\t')
        position = int(position)
        if not keep_n and alt == 'N':
            continue
        if l <= position <= r and coverages[segment] >= cov_threshold:
            outp.write(line)


def main():
    import argparse
    from extra._version import __version__

    parser = argparse.ArgumentParser(description='Filter errors based on positions and segment coverage',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -e File -c File -o File [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-e', '--errors', help='Input alignment errors csv (short format)',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-c', '--coverage', help='Csv file with segment coverages',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-o', '--output', help='Output csv file',
                         type=argparse.FileType('w'), required=True, metavar='File')

    filter_args = parser.add_argument_group('Filter arguments')
    filter_args.add_argument('-r', '--range', metavar=('Int', 'Int'), type=int, nargs=2,
                             help='Positions range (default: %(default)s)', 
                             default=[40, 290])
    filter_args.add_argument('-t', '--threshold', metavar='Int', type=int, default=100,
                             help='Minimal segment coverage (default: %(default)s)')
    filter_args.add_argument('--keep-n', action='store_true',
                             help='Do not filter out mismatches with N as a substitute')
    
    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
    coverages = load_segment_coverage(args.coverage)
    filter_errors(args.errors, args.output,
                  pos_range=args.range, coverages=coverages, cov_threshold=args.threshold, keep_n=args.keep_n)


if __name__ == '__main__':
    main()

