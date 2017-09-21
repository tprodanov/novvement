#!/usr/bin/python3

import argparse
from collections import Counter
import sys

from extra._version import __version__
from extra import nt_string


def str_items(items):
    return ','.join('%s:%s' % (x, y) for x, y in items)


class Combination:
    def __init__(self, segment, mismatches):
        self.segment = segment
        self.mismatches = mismatches
        self.coverage = 0
        self.j_hits = Counter()
        self.joint = Counter()

    def add_hit(self, j_hit):
        self.coverage += 1
        self.j_hits[j_hit] += 1

    def add_joint_mismatch(self, pos_nt):
        self.joint[pos_nt] += 1

    def __len__(self):
        return len(self.mismatches)

    def __str__(self):
        return '%s\t%d\t%d\t%s\t%s' % (self.segment, self.coverage, len(self),
                                       str_items(self.mismatches),
                                       str_items(self.j_hits.most_common()))

    def joint_polymorphisms(self, rate):
        for pos_nt, coverage in self.joint.items():
            if coverage >= rate * self.coverage:
                yield pos_nt, coverage


class Segment:
    def __init__(self, segment):
        self.mismatches = set()
        self.combinations = {}
        self.segment = segment
    
    def add_mismatch(self, pos, nt):
        self.mismatches.add((pos, nt))

    def add_read(self, read, j_hit):
        read = set(read)
        intersection = self.mismatches & read
        intersection = tuple(sorted(intersection))

        if intersection not in self.combinations:
            self.combinations[intersection] = Combination(self.segment, intersection)

        combination = self.combinations[intersection]
        combination.add_hit(j_hit)

        for pos_nt in read - self.mismatches:
            combination.add_joint_mismatch(pos_nt)
        

    def get_combinations(self):
        from operator import attrgetter
        return sorted(self.combinations.values(), reverse=True, key=attrgetter('coverage'))

    def write_combinations(self, length, threshold, outp):
        for combination in self.get_combinations():
            if combination.coverage < threshold or len(combination) < length:
                continue
            outp.write(str(combination))
            outp.write('\n')

    def expanded_polymorphisms(self, threshold, rate, expanded_outp):
        expanded_set = set()

        for combination in self.combinations.values():
            expanded_set |= set(combination.mismatches)
            if combination.coverage < threshold:
                continue
            
            for pos_nt, coverage in sorted(combination.joint_polymorphisms(rate)):
                expanded_outp.write('%s\t%d\t%s\t%d\t%d\n'
                                    % (self.segment, *pos_nt, coverage,
                                       combination.coverage))
                expanded_set.add(pos_nt)
        return sorted(expanded_set)


def load_mismatches(f):
    while next(f).startswith('#'):
        pass

    segments = {}
    for line in f:
        split_line = line.strip().split('\t')
        segment = split_line[0]
        position = int(split_line[1])
        nt = split_line[2]

        if len(split_line) == 3 or split_line[3] == '+':
            if segment not in segments:
                segments[segment] = Segment(segment)
            segments[segment].add_mismatch(position, nt)
    return segments


def load_v_alignment(f, segments, j_hits, margins):
    left, right = margins

    while next(f).startswith('#'):
        pass

    prev_read = None

    for line in f:
        split_line = line.strip().split('\t')
        read = split_line[1]

        if read != prev_read:
            if prev_read and segment in segments and prev_read in j_hits:
                segments[segment].add_read(current_mismatches, j_hits[prev_read])
            prev_read = read
            current_mismatches = []

        segment = split_line[0]
        mism_type = split_line[2]
        pos = int(split_line[3])
        nt = split_line[4]

        if mism_type == 'mismatch' and left <= pos <= right:
            current_mismatches.append((pos, nt))
    
    if prev_read and segment in segments and prev_read in j_hits:
        segments[segment].add_read(current_mismatches, j_hits[prev_read])


def load_j_hits(f):
    j_hit = {}
    for read in nt_string.read_fasta(f):
        j_hit[read.name[1:]] = read.seq
    return j_hit


def write_combinations(segments, length, threshold, outp):
    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('segment\tcoverage\tlength\tcombination\tj_hits\n')
    for segment in segments.values():
        segment.write_combinations(length, threshold, outp)


def write_potential(segments, threshold, rate, potential_outp, expanded_outp):
    potential_outp.write('# %s\n' % ' '.join(sys.argv))
    expanded_outp.write('# %s\n' % ' '.join(sys.argv))
    potential_outp.write('segment\tposition\tnt\n')
    expanded_outp.write('segment\tposition\tnt\tpolym_cov\tcomb_cov\n')
    
    for name, segment in segments.items():
        for pos, nt in segment.expanded_polymorphisms(threshold, rate, expanded_outp):
            potential_outp.write('%s\t%d\t%s\n' % (name, pos, nt))


def main():
    parser = argparse.ArgumentParser(description='Combine potential mismatches into combinations',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -v File -j File -m File [-c File] [-p File [-e File]] [args]')
    input_files = parser.add_argument_group('Input files')
    input_files.add_argument('-v', '--v-alignment', help='Csv file containig V alignment',
                             type=argparse.FileType(), required=True, metavar='File')
    input_files.add_argument('-j', '--j-hit', help='Fasta file with cdrs',
                             type=argparse.FileType(), required=True, metavar='File')
    input_files.add_argument('-m', '--mismatches', help='Csv file contating potential mismatches',
                             type=argparse.FileType(), required=True, metavar='File')

    output_files = parser.add_argument_group('Output files')
    output_files.add_argument('-c', '--combinations', help='Csv combinations output',
                              type=argparse.FileType('w'), metavar='File')
    output_files.add_argument('-p', '--potential', type=argparse.FileType('w'),
                              required=False, metavar='File',
                              help='Csv file with expanded potential polymorphisms (optional)')
    output_files.add_argument('-e', '--expanded', type=argparse.FileType('w'),
                              default='/dev/null', metavar='File',
                              help='Only novel polymorphisms (only if -p File, optional)')

    det_args = parser.add_argument_group('Combinations detection (only if -c File)')
    det_args.add_argument('--length', default=4, type=int, metavar='Int',
                          help='Min combination length (default: 4)')
    det_args.add_argument('--detection-cov', default=15, type=int, metavar='Int',
                          help='Min coverage for combination detection (default: 15)')

    exp_args = parser.add_argument_group('Combinations expansion (only if -p File)')
    exp_args.add_argument('--range', nargs=2, type=int, metavar='Int', default=[60, 290],
                          help='Positions range (default: 60 290)')
    exp_args.add_argument('--expansion-cov', default=50, type=int, metavar='Int',
                          help='Min coverage for combination expansion (default: 50)')
    exp_args.add_argument('--rate', default=0.9, type=float, metavar='Float',
                          help='Polymorphism should appear together with combination\n'
                               'more often than <--rate> * <combination coverage>\n'
                               '(default: 0.9)')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other.add_argument('--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
    
    segments = load_mismatches(args.mismatches)
    j_hits = load_j_hits(args.j_hit)

    if not args.potential:
        args.range = [-float('inf'), float('inf')]

    load_v_alignment(args.v_alignment, segments, j_hits, args.range)
    if args.combinations:
        write_combinations(segments, args.length, args.detection_cov, args.combinations)

    if args.potential:
        write_potential(segments, args.expansion_cov, args.rate,
                        args.potential, args.expanded)
    elif args.expanded.name != '/dev/null':
        print('Warning: -e is set, while -p is unset')


if __name__ == '__main__':
    main()
