#!/usr/bin/python3

import argparse
import sys
from collections import Counter
from collections import defaultdict

from _version import __version__


class Parameters:
    left_margin = 60
    right_margin = 290
    coverage_threshold = 200
    mismatch_rate = 0.1


def read_v_alignment(v_alignment, segment_coverages):
    mismatches = defaultdict(Counter)

    while next(v_alignment).startswith('#'):
            pass

    for line in v_alignment:
        segment, read, error_type, position, nt = line.strip().split('\t')
        position = int(position)
        if error_type == 'mismatch' and segment_coverages[segment] >= Parameters.coverage_threshold and \
                Parameters.left_margin <= position <= Parameters.right_margin:
            mismatches[segment][(position, nt)] += 1
    return mismatches


def filter_mismatches(segment, segment_mismatches, segment_coverages, outp):
    i = 1
    segment_important_mismatches = dict()
    segment_mismatch_indexes = dict()
    for pos_nt, coverage in segment_mismatches.items():
        if coverage >= segment_coverages * Parameters.mismatch_rate:
            outp.write('%s\t%d\t%s\t%d\t%d\n' % (segment, pos_nt[0], pos_nt[1], coverage, segment_coverages))
            segment_important_mismatches[pos_nt] = i
            segment_mismatch_indexes[i] = pos_nt
            i *= 2
    return segment_important_mismatches, segment_mismatch_indexes


def extract_mismatches(mismatches, segment_coverages, outp):
    outp.write('segment\tposition\tnt\tmut_coverage\tsegment_coverage\n')
    for segment, segment_mismatches in mismatches.items():
        filter_mismatches(segment, segment_mismatches, segment_coverages[segment], outp)


def calculate_potential_mismatches(segment_coverage, v_alignment, outp):
    mismatches = read_v_alignment(v_alignment, segment_coverage)
    extract_mismatches(mismatches, segment_coverage, outp)


def load_potential_mismatches(f):
    potential_mismatches = defaultdict(Counter)

    while next(f).startswith('#'):
        pass

    for line in f:
        split_line = line.strip().split()
        segment = split_line[0]
        position = int(split_line[1])
        nt = split_line[2]
        mut_coverage = int(split_line[3])
        # segment_coverage = int(split_line[4])

        potential_mismatches[segment][(position, nt)] = mut_coverage
    return potential_mismatches


def load_segment_coverage(segment_coverage_f):
    segment_coverage = dict()
    for line in segment_coverage_f:
        segment, coverage = line.strip().split()
        segment_coverage[segment] = int(coverage)
    return segment_coverage


def main():
    parser = argparse.ArgumentParser(description='Potential mismatches extraction',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -s File -v File -o File [args]')
    input_files = parser.add_argument_group('Input files')
    input_files.add_argument('-s', '--segment-coverage', help='File with lines <segment> <coverage>',
                             type=argparse.FileType('r'), required=True, metavar='File',
                             dest='segment_coverage')
    input_files.add_argument('-v', '--v-alignment', help='Csv file containig V alignment', dest='v_alignment',
                             type=argparse.FileType('r'), required=True, metavar='File')

    output_files = parser.add_argument_group('Output files')
    output_files.add_argument('-o', '--output', help='Potential mismatches in csv format',
                              type=argparse.FileType('w'), required=True, metavar='File')

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--range', help='Positions range (default: [60, 290])',
                          metavar=('Int', 'Int'), nargs=2, default=[60, 290], type=int)
    optional.add_argument('--coverage', help='Segment coverage threshold (default: 200)',
                          type=int, default=200, metavar='Int')
    optional.add_argument('--rate', help='Mismatch rate threshold (default: 0.1)',
                          type=float, default=0.1, metavar='Float')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other.add_argument('--version', action='version', help='Show version', version=__version__)
    args = parser.parse_args()

    Parameters.left_margin = args.range[0]
    Parameters.right_margin = args.range[1]
    Parameters.coverage_threshold = args.coverage
    Parameters.mismatch_rate = args.rate

    segment_coverage = load_segment_coverage(args.segment_coverage)
    args.output.write('# %s\n' % ' '.join(sys.argv))
    calculate_potential_mismatches(segment_coverage, args.v_alignment, args.output)


if __name__ == '__main__':
    main()
