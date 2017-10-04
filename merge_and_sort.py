#!/usr/bin/env python3

import collections


def significance(labels, min_coverage):
    res = 0
    for label, count in labels.items():
        if count >= min_coverage:
            res += 1
    return res


def obj_significance(min_coverage):
    def inner_function(obj):
        return significance(obj.labels, min_coverage)
    return inner_function


class NovelSequence:
    def __init__(self, seq):
        self.seq = seq
        self.segm_comb = set()
        self.labels = collections.Counter()

    def add_entry(self, segment, combination, labels):
        self.segm_comb.add('%s(%s)' % (segment, combination))
        self.labels.update(labels)

    def to_string(self, min_coverage):
        return '{significance}\t{combination}\t{seq}\t{labels}' \
                    .format(significance=significance(self.labels, min_coverage),
                            combination=','.join(self.segm_comb),
                            seq=self.seq,
                            labels=','.join('%s:%d' % entry for entry in self.labels.items()))


class MergedSequence:
    def __init__(self):
        self.sequences = {}
        self.labels = collections.Counter()
        self.__sorted_with = None
        self.__sorted_sequences = None

    def add_line(self, line):
        self.__sorted_with = None

        line = line.strip()
        segment, combination, labels, seq = line.split('\t')
        if seq not in self.sequences:
            self.sequences[seq] = NovelSequence(seq)
        labels = [entry.split(':') for entry in labels.split(',')]
        labels = {label: int(count) for label, count in labels}

        self.sequences[seq].add_entry(segment, combination, labels)
        self.labels.update(labels)

    def sorted(self, min_coverage):
        if self.__sorted_with != min_coverage:
            self.__sorted_sequences = sorted(self.sequences.values(),
                                             key=obj_significance(min_coverage),
                                             reverse=True)
            self.__sorted_with = min_coverage
        return self.__sorted_sequences

    def to_string(self, min_coverage):
        sorted_sequences = self.sorted(min_coverage)
        return '{significance}\t{combination}\t{seq}\t{labels}' \
                    .format(significance=significance(self.labels, min_coverage),
                            combination=','.join(sorted_sequences[0].segm_comb),
                            seq=sorted_sequences[0].seq,
                            labels=','.join('%s:%d' % entry for entry in self.labels.items()))


def load_components(f):
    while next(f).startswith('#'):
        pass

    return dict(line.strip().split('\t') for line in f)


def merge_combinations(f, components, merged):
    while next(f).startswith('#'):
        pass

    for line in f:
        seq = line.strip().split('\t')[3]
        if seq not in components:
            continue
        merged[components[seq]].add_line(line)


def print_merged(merged, outp, short, min_coverage):
    import sys
    outp.write('# %s\n' % ' '.join(sys.argv))
    
    if not short:
        outp.write('super\t')
    outp.write('significance\tcombinations\tseq\tlabels\n')

    for merged_seq in sorted(merged.values(), key=obj_significance(min_coverage), reverse=True):
        if short:
            outp.write('%s\n' % merged_seq.to_string(min_coverage))
        else:
            outp.write('*\t%s\n' % merged_seq.to_string(min_coverage))
            for novel_seq in merged_seq.sorted(min_coverage):
                outp.write('\t%s\n' % novel_seq.to_string(min_coverage))


def main():
    import argparse
    from extra._version import __version__

    parser = argparse.ArgumentParser(description='Expand candidate polymorphisms',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -i File [File ...] -c File -o File File [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', help='Multiple files containing combined errors',
                         type=argparse.FileType(), required=True, metavar='File', nargs='+')
    io_args.add_argument('-c', '--components', help='File with lines <seq> <its component>',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-o', '--output', help='Two output csv files, short and long',
                         type=argparse.FileType('w'), required=True, metavar='File', nargs=2)

    sign_args = parser.add_argument_group('Significance evaluation')
    sign_args.add_argument('-m', '--min-coverage', help='Min label coverage (default: %(default)s)',
                           type=int, metavar='Int', default=5)
                
    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
    components = load_components(args.components)
    merged = {component: MergedSequence() for component in components.values()}

    for f in args.input:
        merge_combinations(f, components, merged)
    print_merged(merged, args.output[0], short=True, min_coverage=args.min_coverage)
    print_merged(merged, args.output[1], short=False, min_coverage=args.min_coverage)


if __name__ == '__main__':
    main()

