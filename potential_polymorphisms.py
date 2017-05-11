#!/usr/bin/python3

import argparse
import sys
from collections import Counter
from collections import defaultdict
from statistics import median

from extra._version import __version__


class Segment:
    left_margin = None
    right_margin = None
    auto_bad = None
    auto_bad_cov = None
    auto_good = None
    neigh_radius = None
    neigh_mult = None
    neigh_nth = None

    def __init__(self, name, coverage):
        self.name = name
        self.mismatches = Counter()
        self.coverage = coverage

    def add(self, pos, nt):
        if Segment.left_margin <= pos <= Segment.right_margin and nt != 'N':
            self.mismatches[(pos, nt)] += 1

    def get_potential(self):
        highest = [0] * (Segment.right_margin - Segment.left_margin + 1)
        for i in range(Segment.left_margin, Segment.right_margin + 1):
            for nt in 'ACGTN':
                highest[i - Segment.left_margin] = max(highest[i - Segment.left_margin],
                                                       self.mismatches[(i, nt)])

        potential = []
        for (pos, nt), count in self.mismatches.items():
            rate = count / self.coverage
            if rate < Segment.auto_bad or count < Segment.auto_bad_cov:
                continue
            if rate >= Segment.auto_good:
                potential.append((pos, nt))
                continue

            l = pos - Segment.neigh_radius
            r = pos + Segment.neigh_radius
            if l < Segment.left_margin:
                l = Segment.left_margin
                r = l + 2 * Segment.neigh_radius
            if r > Segment.right_margin:
                r = Segment.right_margin
                l = r - 2 * Segment.neigh_radius

            loc = sorted(highest[l - Segment.left_margin : r - Segment.left_margin + 1],
                         reverse=True)
            if count >= loc[Segment.neigh_nth] * Segment.neigh_mult:
                potential.append((pos, nt))
        return sorted(potential)


def read_v_alignment(v_alignment, segments):
    while next(v_alignment).startswith('#'):
            pass

    for line in v_alignment:
        segment, read, error_type, position, nt = line.strip().split('\t')
        position = int(position)
        if error_type != 'mismatch' or segment not in segments:
            continue
        
        segments[segment].add(position, nt)


def write_potential(segments, outp):
    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('segment\tposition\tnt\n')
    for segment in segments.values():
        for pos, nt in segment.get_potential():
            outp.write('%s\t%d\t%s\n' % (segment.name, pos, nt))


def write_full(segments, outp):
    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('segment\tposition\tnt\tpotential\tmut_coverage\tsegment_coverage\n')
    for segment in segments.values():
        potential = set(segment.get_potential())
        for pos_nt, count in sorted(segment.mismatches.items()):
            outp.write('%s\t%d\t%s\t%s\t%d\t%d\n' % (segment.name, *pos_nt,
                    '+' if pos_nt in potential else '-', count, segment.coverage))


def load_segment_coverage(segments_f, threshold):
    segments = {}
    for line in segments_f:
        segment, coverage = line.strip().split()
        coverage = int(coverage)
        if coverage >= threshold:
            segments[segment] = Segment(segment, coverage)
    return segments


def main():
    parser = argparse.ArgumentParser(description='Detecting potential polymorphisms',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -s File -v File -o File [args]')
    input_files = parser.add_argument_group('Input files')
    input_files.add_argument('-s', '--segments', help='File with lines <segment> <coverage>',
                             type=argparse.FileType('r'), required=True, metavar='File')
    input_files.add_argument('-v', '--v-alignment', help='Csv file containig V alignment',
                             type=argparse.FileType('r'), required=True, metavar='File')

    output_files = parser.add_argument_group('Output files')
    output_files.add_argument('-o', '--output', help='Potential polymorphisms in csv format',
                              type=argparse.FileType('w'), required=True, metavar='File')
    output_files.add_argument('--full', help='Output contains more information',
                              action='store_true')

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--range', help='Positions range (default: [45, 290])',
                          metavar=('Int', 'Int'), nargs=2, default=[45, 290], type=int)
    optional.add_argument('--coverage', help='Segment coverage threshold (default: 100)',
                          type=int, default=100, metavar='Int')

    neigh = parser.add_argument_group('Neighborhood arguments')
    neigh.add_argument('--good', help='Polymorphisms that appeared more often are\n'
                                      'automatically treated as potential (default: 0.15)',
                       type=float, default=0.15, metavar='Float')
    neigh.add_argument('--bad', help='Polymorphisms that appeared more rarely are\n'
                                     'automatically filtered out (default: 0.05)',
                       type=float, default=0.05, metavar='Float')
    neigh.add_argument('--bad-cov', default=15, metavar='Int', type=int,
                       help='Polymorphisms that appeared less frequent than\n'
                            '--bad-cov times are automatically filtered out\n'
                            '(default: 15)')
    neigh.add_argument('--mult', help='Potential polymorphism should appear --mult\n'
                                      'times more frequent than the --nth frequent\n'
                                      'in the neighborhood (default: 4)',
                       type=float, default=4, metavar='Float')
    neigh.add_argument('--nth', help='Nth frequent position in the neighborhood for comparison\n'
                                     '0 is the most frequent (useless)\n'
                                     '<--radius> : median\n'
                                     '2 * <--radius> : the least frequent\n'
                                     '(default: 2)',
                       type=int, default=2, metavar='Int')
    neigh.add_argument('--radius', help='Neighborhood radius (default: 5)',
                       type=int, default=5, metavar='Int')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other.add_argument('--version', action='version', help='Show version', version=__version__)
    args = parser.parse_args()

    Segment.left_margin = args.range[0]
    Segment.right_margin = args.range[1]
    Segment.auto_good = args.good
    Segment.auto_bad = args.bad
    Segment.auto_bad_cov = args.bad_cov
    Segment.neigh_nth = args.nth
    Segment.neigh_mult = args.mult
    Segment.neigh_radius = args.radius

    segments = load_segment_coverage(args.segments, args.coverage)
    read_v_alignment(args.v_alignment, segments)
    if args.full:
        write_full(segments, args.output)
    else:
        write_potential(segments, args.output)


if __name__ == '__main__':
    main()
