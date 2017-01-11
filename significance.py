#!/usr/bin/python3

import argparse
from operator import itemgetter
from os import path
import sys

from _version import __version__

class Combination:
    def __init__(self, line, name):
        split_line = line.strip('\n').split('\t')
        self.name = name
        self.segment = split_line[0]
        self.length = split_line[2]
        self.combination = split_line[3]

        if split_line[-1]:
            self.datasets = [(dataset, int(coverage)) for dataset, coverage in
                              [d.split(':') for d in split_line[-1].split(',')]]
        else:
            self.datasets = []
        self.datasets.append((name, int(split_line[1])))
        self.datasets.sort(key=itemgetter(1), reverse=True)

        self.j_segments = [(j_segment, int(coverage)) for j_segment, coverage in
                        [j.split(':') for j in split_line[-2].split(',')]]
        self.j_segments.sort(key=itemgetter(1), reverse=True)

        self.significance = max(sum(cov for _, cov in self.datasets[1:]),
                                sum(cov for _, cov in self.j_segments[1:]))

    def __str__(self):
        return '%s\t%s\t%d\t%s\t%s\t%s\t%s\n' % (self.name, self.segment, self.significance,
                                                 self.length, self.combination,
                                                 ','.join('%s:%d' % j for j in self.j_segments),
                                                 ','.join('%s:%d' % d for d in self.datasets))


def analyze_file(name, combinations_file, combinations_list):
    while next(combinations_file).startswith('#'):
        pass
    for line in combinations_file:
        combinations_list.append(Combination(line, name))


def count_significance(filelist):
    combinations_list = []
    for line in filelist:
        name, filename = line.strip().split()
        filename = path.join(path.dirname(filelist.name), filename)
        with open(filename) as combinations_file:
            analyze_file(name, combinations_file, combinations_list)
    return combinations_list


def write_combinations(combinations_list, min_significance, keep_duplicates, outp):
    combinations_list.sort(key=lambda x: -x.significance)
    outp.write('dataset\tsegment\tsignificance\tlength\tcombination\tj_hit\tdatasets\n')
    seen = set()

    for combination in combinations_list:
        mark = combination.segment + ':' + combination.combination
        if mark in seen:
            continue
        if not keep_duplicates:
            seen.add(mark)

        if combination.significance < min_significance:
            break
        outp.write(str(combination))


def main():
    parser = argparse.ArgumentParser(description='sort combinations by significance'
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    io_args = parser.add_argument_group('Input/output files')
    io_args.add_argument('-i', '--input', help='Input file with lines <name> <combinations.csv path>',
                         type=argparse.FileType('r'), required=True, metavar='File')
    io_args.add_argument('-o', '--output', help='Output in csv format',
                         type=argparse.FileType('w'), required=True, metavar='File')

    opt_args = parser.add_argument_group('Optional arguments')
    opt_args.add_argument('-m', '--min-significance', help='Minimum significance (default: 0)',
                          type=int, default=0, dest='min_significance', metavar='Int')
    opt_args.add_argument('-d', '--keep-duplicates', help='Keep duplicate combinations',
                          action='store_true', dest='keep_duplicates')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other.add_argument('--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
    combinations_list = count_significance(args.input)
    args.output.write('# %s\n' % ' '.join(sys.argv))
    write_combinations(combinations_list, args.min_significance, args.keep_duplicates, args.output)


if __name__ == '__main__':
    main()
