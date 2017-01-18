#!/usr/bin/python3

import argparse
from collections import Counter
from collections import defaultdict
from operator import itemgetter
import os
import re
import sys

from extra import utilities
from extra._version import __version__
import potential_mismatches


def index_mismatches(mismatches):
    mismatches_list = dict()
    mismatches_ind = dict()
    for segment, segment_mismatches in mismatches.items():
        mismatches_list[segment] = list(sorted(segment_mismatches, key=itemgetter(0)))
        mismatches_ind[segment] = dict()
        for i in range(len(mismatches_list[segment])):
            mismatches_ind[segment][mismatches_list[segment][i]] = i
    return mismatches_list, mismatches_ind


def combination_hash(combination, mismatches_ind):
    h = 0
    for mut in combination:
        h += 2 ** (mismatches_ind[mut])
    return h


class Combination:
    def __init__(self, combination):
        self.combination = combination
        self.coverage = 0
        self.j_hit = Counter()
        self.datasets = Counter()


class Combinations:
    def __init__(self, v_alignment, mismatches, j_hit):
        # self.mismatches, self.mismatches_ind = index_mismatches(mismatches)
        self.mismatches = mismatches
        self.mismatches_ind = index_mismatches(mismatches)[1]

        self.mut_combinations = defaultdict(dict)
        for read, segment, read_mismatches in utilities.v_alignment_to_mismatches(v_alignment,
                                                                             self.mismatches):
            h = combination_hash(read_mismatches, self.mismatches_ind[segment])
            if h not in self.mut_combinations[segment]:
                self.mut_combinations[segment][h] = Combination(read_mismatches)

            self.mut_combinations[segment][h].coverage += 1
            self.mut_combinations[segment][h].j_hit[j_hit[read]] += 1

    def write(self, outp, length, cov1, covM, validation):
        outp.write('segment\tcoverage\tlength\tcombination\tj_hit')
        if validation:
            outp.write('\tdatasets')
        outp.write('\n')

        for segment, segment_combinations in self.mut_combinations.items():
            for combination in segment_combinations.values():
                if len(combination.combination) < length:
                    continue
                if len(combination.j_hit) == 1 and combination.coverage < cov1:
                    continue
                if len(combination.j_hit) > 1 and combination.coverage < covM:
                    continue
                outp.write('%s\t%d\t%d\t%s\t%s' %
                           (segment, combination.coverage, len(combination.combination),
                            ','.join('%d:%s' % mut for mut in combination.combination),
                            ','.join('%s:%d' % j for j in combination.j_hit.items())))
                if validation:
                    outp.write('\t')
                    outp.write(','.join('%s:%d' % d for d in combination.datasets.items()))
                outp.write('\n')

    def write_human(self, outp, length, cov1, covM):
        for segment, segment_combinations in self.mut_combinations.items():
            outp.write('-------\n%s\n' % segment)
            for combination in segment_combinations.values():
                if len(combination.combination) < length:
                    continue
                if len(combination.j_hit) == 1 and combination.coverage < cov1:
                    continue
                if len(combination.j_hit) > 1 and combination.coverage < covM:
                    continue

                outp.write('\t%s:\n' % ', '.join('%3d %s' % mut for mut in combination.combination))
                outp.write('\t\tCoverage: %d\n' % combination.coverage)
                outp.write('\t\tLength: %d\n' % len(combination.combination))
                if combination.j_hit:
                    outp.write('\t\tJ hits:\n\t\t\t')
                    outp.write('\n\t\t\t'.join('%s: %d' % j for j in combination.j_hit.items()))
                    outp.write('\n')
                if combination.datasets:
                    outp.write('\t\tDatasets:\n\t\t\t')
                    outp.write('\n\t\t\t'.join('%s: %d' % d for d in combination.datasets.items()))
                    outp.write('\n')

    def validate_v_alignment(self, v_alignment, name):
        for _, segment, read_mismatches in utilities.v_alignment_to_mismatches(v_alignment, self.mismatches):
            h = combination_hash(read_mismatches, self.mismatches_ind[segment])
            if h not in self.mut_combinations[segment]:
                continue

            self.mut_combinations[segment][h].datasets[name] += 1


def load_j_segment(cdr_details):
    j_hit = dict()
    for line in cdr_details:
        read_name, j_segment = line.strip().split('\t')
        j_hit[read_name] = j_segment
    return j_hit


def validate(mut_combinations, validation_f, my_v_alignment):
    for line in validation_f:
        line = line.strip()
        if line[0] == '#':
            continue
        name, v_alignment = re.findall(r"[^ \t:,]+", line)
        full_v_alignment = os.path.dirname(validation_f.name)
        if full_v_alignment:
            full_v_alignment += '/'
        full_v_alignment += v_alignment

        if os.path.abspath(full_v_alignment) == os.path.abspath(my_v_alignment):
            continue
        # print(full_v_alignment)

        with open(full_v_alignment) as v_alignment_f:
            mut_combinations.validate_v_alignment(v_alignment_f, name)


def main():
    parser = argparse.ArgumentParser(description='Combine potential mismatches into combinations',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -v File -j File -m File [-d File] -o File [-u File] [args]')
    input_files = parser.add_argument_group('Input files')
    input_files.add_argument('-v', '--v-alignment', help='Csv file containig V alignment', dest='v_alignment',
                             type=argparse.FileType('r'), required=True, metavar='File')
    input_files.add_argument('-j', '--j-hit', help='File with lines <read name> <j segment>',
                             dest='j_hit', type=argparse.FileType('r'), required=True, metavar='File')
    input_files.add_argument('-m', '--mismatches', help='Csv file contating potential mismatches',
                             type=argparse.FileType('r'), required=True, metavar='File')
    input_files.add_argument('-d', '--datasets', type=argparse.FileType('r'),
                             help='File with lines <name> <v_alignment.csv>\n'
                                  'for the validation (optional)', metavar='File')

    output_files = parser.add_argument_group('Output files')
    output_files.add_argument('-o', '--output', help='Csv combinations output',
                              type=argparse.FileType('w'), required=True, metavar='File')
    output_files.add_argument('-u', '--human', help='Human readable combinations output (optional)',
                              type=argparse.FileType('w'), metavar='File')

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--length', help='Min combination length (default: 2)',
                          default=2, type=int, metavar='Int')
    optional.add_argument('--cov-single-j', help='Min coverage with the single J hit (default: 15)',
                          default=15, type=int, metavar='Int', dest='cov1')
    optional.add_argument('--cov-mult-j', help='Min coverage with multiple J hits (default: 5)',
                          default=5, type=int, metavar='Int', dest='cov')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other.add_argument('--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    mismatches = potential_mismatches.load_potential_mismatches(args.mismatches)
    j_hit = load_j_segment(args.j_hit)
    mut_combinations = Combinations(args.v_alignment, mismatches, j_hit)

    if args.datasets:
        validate(mut_combinations, args.datasets, args.v_alignment.name)

    args.output.write('# %s\n' % ' '.join(sys.argv))
    mut_combinations.write(args.output, args.length, args.cov, args.cov, bool(args.datasets))
    if args.human:
        mut_combinations.write_human(args.human, args.length, args.cov1, args.cov)


if __name__ == '__main__':
    main()
