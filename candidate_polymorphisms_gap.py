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


def analyze_segment(counter, gap):
    from operator import itemgetter
    sorted_errors = counter.most_common()[::-1]

    for i in range(1, len(sorted_errors)):
        if sorted_errors[i][1] - sorted_errors[i - 1][1] >= gap * sorted_errors[i - 1][1] > 1:
            return list(map(itemgetter(0), sorted_errors[i:]))
    return []


def analyze_all(inp, outp, gap):
    import sys

    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('segment\tposition\talt\n')
    errors = load_errors(inp)
    for segment, counter in errors.items():
        for pos, alt in analyze_segment(counter, gap):
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
    detection_args.add_argument('-g', '--gap', metavar='Float', type=float, default=0.5,
                                help='Relative gap size (default: %(default)s)')
    
    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
    analyze_all(args.input, args.output, gap=args.gap)


if __name__ == '__main__':
    main()

