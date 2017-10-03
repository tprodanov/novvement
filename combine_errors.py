#!/usr/bin/env python3


def load_candidate(f):
    import collections
    candidate = collections.defaultdict(set)

    while next(f).startswith('#'):
        pass
    for line in f:
        segment, position, alt = line.strip().split()
        position = int(position)
        candidate[segment].add((position, alt))

    return candidate


class Combination:
    def __init__(self, segment, combination, seq):
        import collections
        self.segment = segment
        self.combination = combination
        self.seq = seq
        self.labels = collections.Counter()

    def add_label(self, label):
        self.labels[label] += 1

    def __str__(self):
        return '{segment}\t{combination}\t{labels}\t{seq}' \
                    .format(segment=self.segment,
                            combination=','.join('%d:%s' % item for item in self.combination)
                                        if self.combination else '*',
                            labels=','.join('%s:%d' % item for item in self.labels.most_common()),
                            seq=self.seq)
                           

class Segment:
    def __init__(self, name, seq, candidate):
        import collections
        self.name = name
        self.seq = seq
        self.combinations = {}
        self.candidate = candidate

    def __create_seq(self, combination):
        new_seq = list(self.seq)
        shift = -1
        for pos, alt in combination:
            if alt.startswith('i'):
                insert = list(alt[1:])
                new_seq = new_seq[ : pos + shift] + insert + new_seq[pos + shift : ]
                shift += len(insert)

            elif alt.startswith('d'):
                deletion_len = int(alt[1:])
                new_seq = new_seq[ : pos + shift] + new_seq[pos + shift + deletion_len : ]
                shift -= deletion_len

            else:
                new_seq[pos + shift] = alt

        return ''.join(new_seq)

    def get_combination(self, read_errors):
        combination = tuple(sorted(alt for alt in read_errors if alt in self.candidate))
        
        if combination in self.combinations:
            return self.combinations[combination]

        self.combinations[combination] = Combination(self.name,
                                                     combination,
                                                     self.__create_seq(combination))
        return self.combinations[combination]

    def get_combinations(self):
        return self.combinations.values()

def load_segments(segments_fa, candidate):
    segment_seqs = {}
    for entry in segments_fa:
        segment_seqs[entry.name[1:]] = entry.seq

    segments = {}
    for segment, segment_candidate in candidate.items():
        segments[segment] = Segment(segment, segment_seqs[segment], segment_candidate)
    return segments


def analyze_all_reads(f, segments, labels):
    import collections
    
    last_read = None
    last_segment = None
    read_errors = None
    
    while next(f).startswith('#'):
        pass
    for line in f:
        read, segment, position, alt = line.strip().split('\t')
        position = int(position)
        if read != last_read:
            if last_segment in segments and last_read in labels:
                segments[last_segment].get_combination(read_errors).add_label(labels[last_read])
            last_read = read
            last_segment = segment
            read_errors = []

        read_errors.append((position, alt))

    if last_segment in segments and last_read in labels:
        segments[last_segment].get_combination(read_errors).add_label(labels[last_read])


def load_labels(f):
    while next(f).startswith('#'):
        pass
    return dict(line.strip().split() for line in f)


def write_combinations(segments, outp):
    import sys
    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('segment\tcombination\tlabels\tseq\n')
    for segment in segments.values():
        for combination in segment.get_combinations():
            outp.write(str(combination))
            outp.write('\n')


def main():
    import argparse
    import os
    from extra._version import __version__
    from extra import nt_string

    parser = argparse.ArgumentParser(description='Convert combinations to sequences and count labels',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -f File -c File -l File -s File -o File')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-f', '--filtered', help='Input csv file with filtered errors',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-c', '--candidate', help='Input csv file with candidate polymorphisms',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-l', '--labels', help='Input csv file with read labels',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-s', '--segments', help='Fasta file with segments',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-o', '--output', help='Output csv file',
                         type=argparse.FileType('w'), required=True, metavar='File')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    candidate = load_candidate(args.candidate)
    segments = load_segments(nt_string.read_fasta(args.segments), candidate)
    labels = load_labels(args.labels)
    analyze_all_reads(args.filtered, segments, labels)
    
    write_combinations(segments, args.output)


if __name__ == '__main__':
    main()

