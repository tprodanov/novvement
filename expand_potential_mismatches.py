#!/usr/bin/python3

import argparse
import sys
from collections import Counter
from collections import defaultdict
from operator import itemgetter
import potential_mismatches
import combinations


class Parameters:
    rate = 0.9
    left_margin = 10
    right_margin = 290
    combination_coverage = 50


class Combination:
    def __init__(self, line, mismatches_ind):
        split_line = line.strip().split('\t')
        self.segment = split_line[0]
        self.coverage = int(split_line[1])

        combination = [(int(position), nt) for position, nt in
                       [mismatch.split(':') for mismatch in split_line[3].split(',')]]
        self.hash = combinations.combination_hash(combination, mismatches_ind[self.segment])
        self.initial = set(combination)
        self.new_mismatches = Counter()


def load_combinations(combinations_file, mismatches_ind):
    combinations_dict = defaultdict(dict)
    while next(combinations_file).startswith('#'):
        pass

    for line in combinations_file:
        combination = Combination(line, mismatches_ind)
        if combination.coverage < Parameters.combination_coverage:
            continue
        combinations_dict[combination.segment][combination.hash] = combination
    return combinations_dict


def process_read(read_mismatches, segment_combinations, segment_mismatches, segment_mismatches_ind):
    potential = []
    for pos, nt in read_mismatches:
        if (pos, nt) in segment_mismatches:
            potential.append((pos, nt))
    potential.sort(key=itemgetter(0))
    h = combinations.combination_hash(potential, segment_mismatches_ind)

    if h not in segment_combinations:
        return
    read_combination = segment_combinations[h]
    for pos, nt in read_mismatches:
        if (pos, nt) not in read_combination.initial:
            read_combination.new_mismatches[(pos, nt)] += 1


def process_v_alignment(v_alignment, combinations_dict, mismatches, mismatches_ind):
    mismatches_coverage = defaultdict(Counter)
    segments_coverage = Counter()

    while next(v_alignment).startswith('#'):
        pass

    prev_read = None
    prev_segment = None
    read_mismatches = None

    for line in v_alignment:
        segment, read, mut_type, position, nt = line.strip().split()
        position = int(position)

        if prev_read != read:
            if read_mismatches and prev_segment in combinations_dict:
                process_read(read_mismatches, combinations_dict[prev_segment],
                             mismatches[prev_segment], mismatches_ind[prev_segment])
            prev_read = read
            prev_segment = segment
            segments_coverage[prev_segment] += 1
            read_mismatches = []

        if mut_type == 'mismatch' and Parameters.left_margin <= position <= Parameters.right_margin:
            mismatches_coverage[segment][(position, nt)] += 1
            read_mismatches.append((position, nt))

    if read_mismatches and prev_segment in combinations_dict:
        process_read(read_mismatches, combinations_dict[prev_segment],
                     mismatches[prev_segment], mismatches_ind[prev_segment])
    segments_coverage[prev_segment] += 1
    return mismatches_coverage, segments_coverage


def print_mismatches(outp, expanded_output, mismatches, combinations_dict,
                    mismatches_coverage, segments_coverage):
    outp.write('segment\tposition\tnt\tmut_coverage\tsegment_coverage\n')
    if expanded_output:
        expanded_output.write('segment\tposition\tnt\t'
                              'cov_in_comb\tcomb_coverage\tmut_coverage\tsegment_coverage\n')

    for segment, segment_coverage in segments_coverage.items():
        if segment not in mismatches:
            continue
        segment_mismatches = set(mismatches[segment])

        for combination in combinations_dict[segment].values():
            for new_mut, cov in combination.new_mismatches.items():
                if cov >= Parameters.rate * combination.coverage:
                    segment_mismatches.add(new_mut)
                    if expanded_output:
                        expanded_output.write('%s\t%d\t%s\t%d\t%d\t%d\t%d\n' %
                                              (segment, *new_mut, cov, combination.coverage,
                                               mismatches_coverage[segment][new_mut], segment_coverage))

        for mut in segment_mismatches:
            outp.write('%s\t%d\t%s\t%d\t%d\n' % (segment, *mut, mismatches_coverage[segment][mut],
                                                 segment_coverage))


def main():
    parser = argparse.ArgumentParser(description='expand potential mismatches', add_help=False)
    input_files = parser.add_argument_group('input files')
    input_files.add_argument('-v', '--v-alignment', help='v_alignment.csv', dest='v_alignment',
                             type=argparse.FileType('r'), required=True, metavar='V')
    input_files.add_argument('-m', '--mismatches', help='potential_mismatches.csv',
                             type=argparse.FileType('r'), required=True, metavar='M')
    input_files.add_argument('-c', '--combinations', help='combinations.csv',
                             type=argparse.FileType('r'), required=True, metavar='C')

    output_files = parser.add_argument_group('output files')
    output_files.add_argument('-o', '--output', help='expanded potential mismatches in csv format',
                              type=argparse.FileType('w'), required=True, metavar='O')
    output_files.add_argument('-e', '--expanded', help='only expanded mismatches (optional)',
                              type=argparse.FileType('w'), metavar='E')

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('--range', help='positions range (default: [60, 290])',
                          metavar=('L', 'R'), nargs=2, default=[60, 290])
    optional.add_argument('--coverage', help='combination coverage threshold (default: 50)',
                          type=int, default=50)
    optional.add_argument('--rate', help='mismatch coverage threshold (rate to combination coverage)'
                                          ' (default: 0.9)',
                          type=float, default=0.9)

    other = parser.add_argument_group('other arguments')
    other.add_argument('-h', '--help', action='help', help='show this help message and exit')

    args = parser.parse_args()

    Parameters.left_margin = int(args.range[0])
    Parameters.right_margin = int(args.range[1])
    Parameters.combination_coverage = args.coverage
    Parameters.rate = args.rate

    mismatches = potential_mismatches.load_potential_mismatches(args.mismatches)
    mismatches_ind = combinations.index_mismatches(mismatches)[1]

    combinations_dict = load_combinations(args.combinations, mismatches_ind)
    mismatches_coverage, segments_coverage = process_v_alignment(args.v_alignment, combinations_dict,
                                                             mismatches, mismatches_ind)
    args.output.write('# %s\n' % ' '.join(sys.argv))
    if args.expanded:
        args.expanded.write('# %s\n' % ' '.join(sys.argv))
    print_mismatches(args.output, args.expanded, mismatches, combinations_dict,
                    mismatches_coverage, segments_coverage)


if __name__ == '__main__':
    main()
