#!/usr/bin/env python3


import collections


def load_candidate(f):
    candidate = collections.defaultdict(set)

    while next(f).startswith('#'):
        pass
    for line in f:
        segment, position, alt = line.strip().split()
        position = int(position)
        candidate[segment].add((position, alt))

    return candidate


class ReadsSubset:
    def __init__(self, segment, candidate, reads, coverage=None):
        self.segment = segment
        self.candidate = candidate
        self.reads = reads
        if coverage is not None:
            self.coverage = coverage
        else:
            self.coverage = len(self.reads)

        self.categorized = collections.defaultdict(list)
        self.__categorize_errors()

    def __read_candidate_errors(self, errors):
        return tuple(sorted(error for error in errors if error in self.candidate))

    def __categorize_errors(self):
        self.categorized.clear()
        for name, errors in self.reads:
            read_category = self.__read_candidate_errors(errors)
            if read_category:
                self.categorized[read_category].append((name, errors))

    def __best_category(self):
        m = 0
        best = None
        for category, reads in self.categorized.items():
            if len(reads) > m and category:
                m = len(reads)
                best = category
        return best

    def __reduce_reads(self, category_reads):
        from operator import itemgetter
        names = set(map(itemgetter(0), category_reads))
        return [read for read in self.reads if read[0] not in names]

    def split_if_needed(self, coverage):
        category = self.__best_category()
        category_reads = self.categorized[category]

        if len(category_reads) < coverage:
            return None

        self.reads = self.__reduce_reads(category_reads)
        self.candidate -= set(category)
        self.coverage -= len(category_reads)
        self.__categorize_errors()

        return ReadsSubset(self.segment, self.candidate, category_reads)

    def consensus(self, min_ratio=0.51):
        frequences = collections.Counter()
        for _, read_errors in self.reads:
            for error in read_errors:
                frequences[error] += 1

        result = []
        for error, count in frequences.items():
            if count >= min_ratio * self.coverage:
                result.append(error)
        return sorted(result)


def create_seq(seq, combination):
    new_seq = list(seq)
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


def load_segments(segments_fa):
    return {entry.name[1:]: entry.seq for entry in segments_fa}


def load_reads(f, labels):
    reads = collections.defaultdict(list)
    
    last_read = None
    last_segment = None
    read_errors = None
    
    while next(f).startswith('#'):
        pass
    for line in f:
        read, segment, position, alt = line.strip().split('\t')
        position = int(position)
        if read != last_read:
            if last_read in labels:
                reads[last_segment].append((last_read, read_errors))
            last_read = read
            last_segment = segment
            read_errors = []

        read_errors.append((position, alt))

    if last_read in labels:
        reads[last_segment].append((last_read, read_errors))
    return reads

def load_labels(f):
    while next(f).startswith('#'):
        pass
    return dict(line.strip().split() for line in f)


def initialize_reads_subsets(reads, segments, candidate):
    reads_subsets = {}
    for segment in segments:
        if segment not in reads:
            continue
        reads_subsets[segment] = ReadsSubset(segment,
                                             candidate.get(segment, set()),
                                             reads[segment])
    return reads_subsets


def traverse_reads_subsets(reads_subsets, coverage):
    stack = list(reads_subsets.values())
    result = []

    while stack:
        current = stack.pop()
        if current.coverage < coverage:
            result.append(current)
            continue
        
        novel = current.split_if_needed(coverage)
        if novel:
            stack.append(current)
            stack.append(novel)
        else:
            result.append(current)

    return result


def reads_labels(reads, labels):
    reads_labels = collections.Counter()
    for name, _ in reads:
        reads_labels[labels[name]] += 1
    return ','.join('%s:%d' % item for item in reads_labels.most_common())


def write_subsets(reads_subsets, labels, segments, outp, coverage, consensus_ratio):
    import sys
    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('segment\tcombination\tlabels\tseq\n')
    for subset in reads_subsets:
        if subset.coverage < coverage:
            continue

        consensus = subset.consensus(consensus_ratio)
        seq = create_seq(segments[subset.segment], consensus)

        outp.write('{segment}\t{combination}\t{labels}\t{seq}\n'
                .format(segment=subset.segment,
                        combination=','.join('%d:%s' % item for item in consensus) if consensus else '*',
                        labels=reads_labels(subset.reads, labels),
                        seq=seq))


def main():
    import argparse
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

    sc_args = parser.add_argument_group('Split and consensus arguments')
    sc_args.add_argument('--coverage', type=int, metavar='Int', default=100,
                         help='Minimal coverage of a reads subset (default: %(default)s')
    sc_args.add_argument('--cons-ratio', type=float, metavar='Float', default=0.51,
                         help='Ratio for a new polymorphisms to be chosen.\n'
                         'Should be higher than 0.5 (default: %(default)s)')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    candidate = load_candidate(args.candidate)
    segments = load_segments(nt_string.read_fasta(args.segments))
    labels = load_labels(args.labels)
    reads = load_reads(args.filtered, labels)

    reads_subsets = initialize_reads_subsets(reads, segments, candidate)
    reads_subsets = traverse_reads_subsets(reads_subsets, args.coverage)
    write_subsets(reads_subsets, labels, segments, args.output,
                  args.coverage, args.cons_ratio)


if __name__ == '__main__':
    main()

