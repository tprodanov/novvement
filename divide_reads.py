#!/usr/bin/env python3


import math
import collections
import itertools


def tuple_array_search(arr, x, ind=0):
    for v in arr:
        if v[ind] == x:
            return v
    return None


class Read:
    def __init__(self, name, errors, label):
        self.name = name
        self.label = label

        if errors == '*':
            self.original_errors = []
        else:
            errors = [error.split(':') for error in errors.split(',')]
            self.original_errors = [(int(pos), alt) for pos, alt in errors]
        self.update([])

    def update(self, consensus):
        errors = []
        for pos, nt in consensus:
            if tuple_array_search(self.original_errors, pos) is None:
                errors.append((pos, '*'))
        
        for error in self.original_errors:
            if error not in consensus:
                errors.append(error)
        self.current_errors = sorted(errors)

    def str(self, segment, subset):
        errors = ','.join('%d:%s' % error for error in self.original_errors) if self.original_errors else '*'
        return '%s\t%s\t%s\t%s\t%s' % (self.name, segment, errors, self.label, subset)



def log_factorial(n):
    if n >= len(log_factorial.values):
        for i in range(len(log_factorial.values), n + 1):
            log_factorial.values.append(log_factorial.values[-1] + math.log(i))
    return log_factorial.values[n]

log_factorial.values = [math.log(math.factorial(i)) for i in range(10)]


def consensus(reads):
    counter = collections.Counter()
    for read in reads:
        for error in read.original_errors:
            counter[error] += 1
    s = len(reads) // 2

    res = []
    for error, count in counter.items():
        if count > s:
            res.append(error)
    return sorted(res)


class ReadsSubset:
    def __init__(self, segment, reads):
        self.segment = segment
        self.reads = reads
        self.consensus = set(consensus(self.reads))
        for read in self.reads:
            read.update(self.consensus)

    def top_pairs(self, min_subset_size):
        counter = collections.Counter()
        for read in self.reads:
            for err1, err2 in itertools.combinations(read.current_errors, 2):
                counter[(err1, err2)] += 1

        for pair, count in counter.most_common():
            if count < min_subset_size or len(self.reads) - count < min_subset_size:
                continue
            yield pair, count

    def significance(self, n1, n2, n12):
        n = len(self.reads)
        base = sum(map(log_factorial, [n1, n - n1, n2, n - n2])) - log_factorial(n)

        s = 0
        for i in range(n12, min(n1, n2) + 1):
            current = sum(map(log_factorial, [i, n1 - i, n2 - i, n - n1 - n2 + i]))
            s += math.exp(base - current)
        return s

    def count_reads_with(self, error):
        return sum(error in read.current_errors for read in self.reads)

    def divide_by(self, err1, err2):
        reads1 = []
        reads2 = []
        for read in self.reads:
            if err1 in read.current_errors and err2 in read.current_errors:
                reads1.append(read)
            else:
                reads2.append(read)
        return reads1, reads2

    def divide_if_possible(self, significance, min_subset_size, max_pairs):
        pairs = list(itertools.islice(self.top_pairs(min_subset_size), max_pairs))
        for (err1, err2), count in pairs:
            n1 = self.count_reads_with(err1)
            n2 = self.count_reads_with(err2)
            s = self.significance(n1, n2, count)
            if s > significance / len(pairs):
                continue

            reads1, reads2 = self.divide_by(err1, err2)
            assert len(reads1) == count

            return ReadsSubset(self.segment, reads1), ReadsSubset(self.segment, reads2)
        return None

    def labels(self):
        lengths = set()
        for read in self.reads:
            lengths.add(len(read.label))
        return len(lengths)

    def write(self, index, read_outp, summary_outp):
        name = '%s:%d' % (self.segment, index)
        labels = self.labels()
        long_name = '"%s:%2d. Coverage: %d. Labels: %d"' % (self.segment, index, len(self.reads), labels)

        for read in self.reads:
            read_outp.write(read.str(self.segment, long_name))
            read_outp.write('\n')

        cons = sorted(self.consensus)
        cons = ','.join('%d:%s' % error for error in cons) if cons else '*'
        summary_outp.write('%s\t%d\t%d\t%s\n' % (name, len(self.reads), labels, cons))


def divide_while_possible(segment, reads, significance, min_subset_size, max_pairs):
    subsets = [ReadsSubset(segment, reads)]
    result = []
    while subsets:
        subset = subsets.pop()
        division = subset.divide_if_possible(significance, min_subset_size, max_pairs)
        if division is None:
            result.append(subset)
        else:
            subsets.extend(division)
    return result


def write_subsets(segment, subsets, read_outp, summary_outp):
    for i, subset in enumerate(sorted(subsets, key=lambda x: -len(x.reads))):
        subset.write(i + 1, read_outp, summary_outp)


def load_errors(f):
    while next(f).startswith('#'):
        pass

    reads = collections.defaultdict(list)
    for line in f:
        read, segment, errors, label = line.strip().split('\t')
        reads[segment].append(Read(read, errors, label))
    return reads


def main():
    import argparse
    from extra._version import __version__
    import sys

    parser = argparse.ArgumentParser(description='Divide read sets into subsets',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='-i File -r File -s File [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', type=argparse.FileType(), metavar='File',
                         help='Csv file with reads merged by label')
    io_args.add_argument('-r', '--reads', type=argparse.FileType('w'), metavar='File',
                         help='Output csv file with information about reads')
    io_args.add_argument('-s', '--summary', type=argparse.FileType('w'), metavar='File',
                         help='Output csv file with subsets summary')

    div_args = parser.add_argument_group('Division arguments')
    div_args.add_argument('-S', '--significance', type=float, metavar='Float', default=.001,
                          help='Division significance (default: %(default)s)')
    div_args.add_argument('-p', '--pairs', type=int, metavar='Int', default=5,
                          help='Maximum number of evaluated pairs of errors (default: %(default)s)')
    div_args.add_argument('-c', '--coverage', type=int, metavar='Int', default=50,
                          help='Minimum subset size (default: %(default)s)')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    args.reads.write('# %s\n' % ' '.join(sys.argv))
    args.reads.write('read\tsegment\terrors\tlabel\tsubset\n')
    args.summary.write('# %s\n' % ' '.join(sys.argv))
    args.summary.write('subset\tcoverage\tlabels\tconsensus\n')

    reads = load_errors(args.input)
    for segment, segm_reads in reads.items():
        subsets = divide_while_possible(segment, segm_reads, args.significance, args.coverage, args.pairs)
        write_subsets(segment, subsets, args.reads, args.summary)


if __name__ == '__main__':
    main()
