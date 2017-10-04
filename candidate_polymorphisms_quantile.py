#!/usr/bin/env python3


def load_errors(f):
    import collections
    errors = collections.defaultdict(collections.Counter)
    
    while next(f).startswith('#'):
        pass
    for line in f:
        read, segment, position, alt = line.strip().split('\t')
        errors[segment][(int(position), alt)] += 1

    return errors


def analyze_segment(counter, quantile, multiplier):
    sorted_errors = counter.most_common()
    q = sorted_errors[(100 - quantile) * len(sorted_errors) // 100][1]
    
    for pos_alt, count in sorted_errors:
        if count >= multiplier * q:
            yield pos_alt
        else:
            break


def analyze_all(inp, outp, quantile, multiplier):
    import sys

    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('segment\tposition\talt\n')
    errors = load_errors(inp)
    for segment, counter in errors.items():
        for pos, alt in analyze_segment(counter, quantile, multiplier):
            outp.write('%s\t%d\t%s\n' % (segment, pos, alt))


def main():
    import argparse
    from extra._version import __version__

    parser = argparse.ArgumentParser(description='Detect candidate polymorphisms',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -i File -o File [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', help='Input alignment errors csv (short format)',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-o', '--output', help='Output csv file',
                         type=argparse.FileType('w'), required=True, metavar='File')

    detection_args = parser.add_argument_group('Detection arguments')
    detection_args.add_argument('-m', '--multiplier', metavar='Float', type=float, default=4,
                                help='Quantile multiplier:\n'
                                'We select positions with number of errors at least\n'
                                '<--multiplier> * <-q quantile of errors per position>\n'
                                '(default: %(default)s)')
    detection_args.add_argument('-q', '--quantile', metavar='Float', type=float, default=95,
                                help='Quantile used for comparison (default: %(default)s)')
    
    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
    analyze_all(args.input, args.output, quantile=args.quantile, multiplier=args.multiplier)


if __name__ == '__main__':
    main()

