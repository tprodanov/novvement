#!/usr/bin/python3

import argparse
import sys
import itertools

import combinations_to_segments
from extra._version import __version__

def distance(seq1, seq2):
    return sum(nt1 != nt2 for nt1, nt2 in itertools.zip_longest(seq1, seq2))


def min_distance(seq, source):
    nearest = None
    dist = None

    for entry in source:
        seq2 = entry[0]
        d = distance(seq, seq2)

        if dist is None or d < dist:
            dist = d
            nearest = entry[1:]
    return (dist, *nearest) if dist is not None else None


def keep(name, seq, significance, source, target, args):
    if args.source_dist:
        s_dist, segment = min_distance(seq, source)
        if s_dist < args.source_dist:
            args.log.write('%s is too close to %s (distance %d)\n' % (name, segment, s_dist))
            return 0

    if not args.target_dist:
        return s_dist

    if not target:
        return s_dist

    dist, nearest_significance, oth = min_distance(seq, target)
    if dist == 0:
        args.log.write('%s is equal to %s\n' % (name, oth))
        return 0
    if dist >= args.target_dist:
        return s_dist
    if significance * args.target_mf < nearest_significance:
        args.log.write('%s is too close to %s (distance %d, mf %d / %d = %.2f)\n'
                       % (name, oth, dist, nearest_significance, significance, nearest_significance / significance))
        return 0
    return s_dist


def run(args):
    segments = combinations_to_segments.read_segments(args.v_segments)
    inp = args.combinations
    outp = args.output
    l, r = args.range
    l -= 1

    source = [(segment.seq[l:r], segment.name) for segment in segments.values()]
    target = []

    line = next(inp)
    while line.startswith('#'):
        line = next(inp)
    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write(line)

    for i, line in enumerate(inp):
        split_line = line.strip().split('\t')
        segment = split_line[0]
        significance = int(split_line[1])
        combination = split_line[3]
        seq = combinations_to_segments.generate_novel_segment(segments[segment].seq, combination)[l:r]

        name = '%s-M%d' % (segment, i + 1)

        s_dist = keep(name, seq, significance, source, target, args)
        if s_dist:
            split_line[2] = str(s_dist)
            outp.write('\t'.join(split_line))
            outp.write('\n')
            target.append((seq, significance, name))


def main():
    parser = argparse.ArgumentParser(description='filter out too similar combinations',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -c File -v File -o File [-l File] [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-c', '--combinations', help='Input: target combinations csv file',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-v', '--v-segments', help='Input: source v segments',
                         type=argparse.FileType(), required=True, metavar='File', dest='v_segments')
    io_args.add_argument('-o', '--output', help='Output: filtered combinations csv file',
                         type=argparse.FileType('w'), metavar='File', required=True)
    io_args.add_argument('-l', '--log', help='Output: log file (optional)',
                         type=argparse.FileType('w'), metavar='File', default='/dev/null')

    filter_args = parser.add_argument_group('Filter arguments')
    filter_args.add_argument('--range', help='Positions range (default: [60, 290])',
                             metavar=('Int', 'Int'), type=int, nargs=2, default=[60, 290])
    filter_args.add_argument('--source-dist', help='Minimum required distance to\n'
                                                   'any source segment (default: 2)',
                             metavar='Int', dest='source_dist', type=int, default=2)
    filter_args.add_argument('--target-dist', help='Minimum reliable distance to\n'
                                                   'any target combination (default: 3)',
                             metavar='Int', dest='target_dist', type=int, default=3)
    filter_args.add_argument('--target-mf', help='Dignificance multiplication factor\n'
                                                 'to filter out combination with unreliable\n'
                                                 'target distance (default: 4)',
                             metavar='Float', dest='target_mf', type=float, default=4)

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other.add_argument('--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
