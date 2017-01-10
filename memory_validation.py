#!/usr/bin/python3

import argparse
import re
import math
import os
import sys
from collections import Counter
from operator import itemgetter
from extra import nt_string


class Parameters:
    left_margin = 60
    right_margin = 290
    neighborhood = 3
    neighborhood_ratio = 1.1
    neighborhood_auto_good = .2
    neighborhood_auto_bad = .8


def load_combinations(file):
    while next(file).startswith('#'):
        pass

    combinations = []
    for line in file:
        split_line = line.strip('\n').split('\t')
        segment = split_line[0]
        combination = split_line[3]
        positions = set(int(mut.split(':')[0]) for mut in combination.split(',')) if combination else set()
        combinations.append((segment, positions))
    return combinations


def add_read(fasta, segments_coverage, mismatches, combinations, mut_combinations_count):
    read = next(fasta)
    segment = next(fasta)

    segment_name = re.search('GENE:([^|]*)\|', segment.name).group(1)

    m = re.search('-M([0-9]*)', segment_name)
    if not m:
        return
    segment_number = int(m.group(1)) - 1
    segments_coverage[segment_number] += 1

    position = 0
    combination_mismatches = 0
    for i in range(len(read.seq)):
        if segment.seq[i] != '-':
            position += 1
        if segment.seq[i] == '-' or read.seq[i] == '-':
            continue

        if Parameters.left_margin <= position <= Parameters.right_margin and segment.seq[i] != read.seq[i]:
            mismatches[segment_number][position - Parameters.left_margin][read.seq[i]] += 1
            if position in combinations[segment_number][1]:
                combination_mismatches += 1
    mut_combinations_count[segment_number][combination_mismatches] += 1


def analyze_fasta(fasta_file, combinations):
    segments_coverage = [0] * len(combinations)
    mismatches = [[Counter() for _ in range(Parameters.right_margin - Parameters.left_margin + 1)]
                  for _ in range(len(combinations))]
    mut_combinations_count = [Counter() for _ in range(len(combinations))]

    try:
        fasta = nt_string.read_fasta(fasta_file)
        while True:
            add_read(fasta, segments_coverage, mismatches, combinations, mut_combinations_count)
    except StopIteration:
        pass
    return segments_coverage, mismatches, mut_combinations_count


def information_content(counter, ratio):
    n = sum(counter.values()) - counter['N']
    if not n:
        return 0

    h = 0.
    for nt, v in counter.items():
        if nt != 'N':
            h -= (v / n) * math.log(v / n, 2)
    return (math.log(3, 2) - h) * ratio


def print_mismatches(segments_coverage, mismatches, combinations, outp):
    ratio_information = [[] for _ in range(len(combinations))]
    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('segment\tnovel\tposition\tin_comb\tsegment_cov\tratio\tinformation\n')
    for i in range(len(combinations)):
        if not segments_coverage[i]:
            continue
        segment, positions = combinations[i]
        prefix = '%s-M%d\t%d' % (segment, i + 1, i + 1)

        for pos in range(Parameters.left_margin, Parameters.right_margin + 1):
            outp.write('%s\t%d\t' % (prefix, pos))
            outp.write('+' if pos in positions else '-')
            counter = mismatches[i][pos - Parameters.left_margin]
            ratio = sum(counter.values()) / segments_coverage[i]
            information = information_content(counter, ratio)

            outp.write('\t%d\t%.3f\t%.3f\n' % (segments_coverage[i], ratio, information))
            ratio_information[i].append((ratio, information))
    return ratio_information


def print_mut_count(segments_coverage, mut_combinations_count, combinations, outp):
    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('segment\tnovel\tlength\tmut_count\tcoverage\tsegment_cov\tratio\n')

    for i in range(len(combinations)):
        if not segments_coverage[i]:
            continue
        segment, positions = combinations[i]
        prefix = '%s-M%d\t%d\t%d' % (segment, i + 1, i + 1, len(positions))

        for count, coverage in mut_combinations_count[i].most_common():
            outp.write('%s\t%d\t%d\t%d\t%.3f\n' % (prefix, count, coverage, segments_coverage[i],
                                                   coverage / segments_coverage[i]))


class KeyPositionsOffRange(Exception):
    pass


class SegmentScore:
    def __init__(self, combination_number, coverage, mut_combinations_count,
                 mismatches, ratio_information, combination):
        self.segment, self.positions = combination

        self.combination_number = combination_number
        self.coverage = coverage
        self.highest_mut_count = mut_combinations_count.most_common(1)[0]

        key_mut_sum = 0
        oth_mut_sum = 0
        for pos in range(Parameters.left_margin, Parameters.right_margin + 1):
            if pos in self.positions:
                key_mut_sum += sum(mismatches[pos - Parameters.left_margin].values())
            else:
                oth_mut_sum += sum(mismatches[pos - Parameters.left_margin].values())
        self.key_mut_ratio = key_mut_sum / len(self.positions) / coverage
        self.oth_mut_ratio = oth_mut_sum / (Parameters.right_margin +
                                            1 - Parameters.left_margin - len(self.positions)) / coverage

        self.count_good_positions(ratio_information)
        self.count_score()

    def count_good_positions(self, ratio_information):
        self.good_positions = 0
        for pos in self.positions:
            if pos < Parameters.left_margin or pos > Parameters.right_margin:
                raise KeyPositionsOffRange
            if ratio_information[pos - Parameters.left_margin][0] <= Parameters.neighborhood_auto_good:
                self.good_positions += 1
                continue
            if ratio_information[pos - Parameters.left_margin][0] >= Parameters.neighborhood_auto_bad:
                continue

            for i in range(max(Parameters.left_margin, pos - Parameters.neighborhood),
                           min(Parameters.right_margin, pos + Parameters.neighborhood) + 1):
                if i in self.positions:
                    continue
                if ratio_information[pos - Parameters.left_margin][0] <= ratio_information[i - Parameters.left_margin][0] * \
                        Parameters.neighborhood_ratio and ratio_information[pos - Parameters.left_margin][1] <= \
                        ratio_information[i - Parameters.left_margin][1] * Parameters.neighborhood_ratio:
                    self.good_positions += 1
                    break

    def count_score(self):
        self.score = 2 / (1 + math.exp(-self.coverage / 20)) - 1
        if self.highest_mut_count[0] == 0:
            self.score += .5
            if self.highest_mut_count[1] >= 20:
                self.score += .5
        self.score += 1 - self.key_mut_ratio
        self.score += self.good_positions / len(self.positions)

    def __str__(self):
        return '%s-M%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%d\t%.3f\t%.3f' \
               % (self.segment, self.combination_number + 1, self.combination_number + 1,
                  len(self.positions), self.coverage, self.highest_mut_count[0],
                  self.highest_mut_count[1], self.key_mut_ratio, self.oth_mut_ratio,
                  self.good_positions, self.good_positions / len(self.positions), self.score)


def print_segment_scores(combinations, segment_coverages, mut_combinations_count,
                         mismatches, ratio_information, db_name, outp):
    segment_scores = []
    for i in range(len(combinations)):
        if not segment_coverages[i]:
            continue
        segment_scores.append(SegmentScore(i, segment_coverages[i], mut_combinations_count[i],
                                           mismatches[i], ratio_information[i], combinations[i]))
    segment_scores.sort(key=lambda x: x.score, reverse=True)

    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('# hmc = mismatch count with the highest coverage\n')
    outp.write('# hmcc = it\'s coverage\n')
    outp.write('segment\tnovel\tlen\tcov\thmc\thmcc\t'
               'key_ratio\toth_ratio\tgood_pos\tgood_pos_ratio\tscore\tdb_name\n')
    for segment_score in segment_scores:
        outp.write(str(segment_score))
        outp.write('\t%s\n' % db_name)


def main():
    parser = argparse.ArgumentParser(description='validate novel segments', add_help=False)
    io_args = parser.add_argument_group('input/output arguments')
    io_args.add_argument('-v', '--v-alignment', help='v_alignment fasta file', dest='v_alignment',
                         type=argparse.FileType('r'), required=True, metavar='FILE')
    io_args.add_argument('-c', '--combinations', help='combinations csv file',
                         type=argparse.FileType('r'), required=True, metavar='FILE')
    io_args.add_argument('-n', '--name', help='this dataset name', required=True, metavar='FILE')
    io_args.add_argument('-d', '--directory', help='output directory', required=True, metavar='FILE')

    val_args = parser.add_argument_group('validation arguments')
    val_args.add_argument('--range', help='positions range (default: [60, 290])',
                          metavar=('INT', 'INT'), nargs=2, default=[60, 290], type=int)
    val_args.add_argument('--neigh', help='position neighborhood size (default: 3)',
                          metavar='INT', default=3, type=int)
    val_args.add_argument('--neigh-overhead', help='acceptable mismatch rate overhead '
                                                   'in the position neighborhood (default: 1.1)',
                          metavar='FLOAT', default=1.1, type=float, dest=neigh_overhead)
    val_args.add_argument('--neigh-good', help='position is treated good, if its mismatch rate is less '
                                               'than --neigh-good (default: 0.2)',
                          metavar='FLOAT', default=0.2, type=float, dest=neigh_good)
    val_args.add_argument('--neigh-bad', help='position is treated good, if its mismatch rate is more '
                                               'than --neigh-bad (default: 0.8)',
                          metavar='FLOAT', default=0.8, type=float, dest=neigh_good)

    other = parser.add_argument_group('other arguments')
    other.add_argument('-h', '--help', action='help', help='show this help message and exit')
    args = parser.parse_args()

    Parameters.left_margin, Parameters.right_margin = args.range
    Parameters.neighborhood = args.neigh
    Parameters.neighborhood_ratio = args.neigh_overhead
    Parameters.neighborhood_auto_good = args.neigh_good
    Parameters.neighborhood_auto_bad = args.neigh_bad

    directory = args.directory
    if not os.path.exists(directory):
        os.mkdir(directory)

    combinations = load_combinations(args.combinations)
    segments_coverage, mismatches, mut_combinations_count = analyze_fasta(args.v_alignment, combinations)
    with open(os.path.join(directory, 'positions.csv'), 'w') as positions:
        ratio_information = print_mismatches(segments_coverage, mismatches, combinations, positions)
    with open(os.path.join(directory, 'mut_count.csv'), 'w') as mut_count:
        print_mut_count(segments_coverage, mut_combinations_count, combinations, mut_count)
    with open(os.path.join(directory, 'segments_summary.csv'), 'w') as segments_summary:
        print_segment_scores(combinations, segments_coverage, mut_combinations_count, mismatches,
                             ratio_information, args.name, segments_summary)


if __name__ == '__main__':
    main()
