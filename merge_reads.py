#!/usr/bin/env python3


from collections import Counter


class MergedRead:
    def __init__(self, segment, label):
        self.segment = segment
        self.label = label
        self.coverage = 0
        self.errors = Counter()

    def add_read(self, errors):
        self.coverage += 1
        if errors == '*':
            return

        for error in errors.split(','):
            pos, nt = error.split(':')
            pos = int(pos)
            self.errors[(pos, nt)] += 1

    def consensus(self):
        threshold = self.coverage // 2
        res = []
        for error, count in self.errors.items():
            if count > threshold:
                res.append(error)
        return sorted(res)

    def str(self, index):
        name = 'read_%d___size_%d' % (index, self.coverage)
        cons = self.consensus()
        errors = ','.join('%d:%s' % error for error in cons) if cons else '*'
        return '%s\t%s\t%s\t%s' % (name, self.segment, errors, self.label)


def load_labels(f):
    while next(f).startswith('#'):
        pass

    return dict(line.strip().split('\t') for line in f)


def load_errors(f, labels):
    while next(f).startswith('#'):
        pass

    merged_reads = {}
    for line in f:
        read, segment, errors = line.strip().split('\t')

        if read not in labels:
            continue
        label = labels[read]

        if label not in merged_reads:
            merged_reads[label] = MergedRead(segment, label)
        merged_reads[label].add_read(errors)
    return merged_reads


def write_merged_reads(merged_reads, outp):
    import sys
    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('read\tsegment\terrors\tlabel\n')

    for i, read in enumerate(merged_reads.values()):
        outp.write(read.str(i + 1))
        outp.write('\n')


def main():
    import argparse
    from extra._version import __version__

    parser = argparse.ArgumentParser(description='Merge reads with the same CDR3 labels',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='-e File -l File -o File')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-e', '--errors', metavar='File', type=argparse.FileType(),
                         help='Csv file with errors')
    io_args.add_argument('-l', '--labels', metavar='File', type=argparse.FileType(),
                         help='Csv file with read labels')
    io_args.add_argument('-o', '--output', metavar='File', type=argparse.FileType('w'),
                         help='Output csv file')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    labels = load_labels(args.labels)
    merged_reads = load_errors(args.errors, labels)
    write_merged_reads(merged_reads, args.output)
 

if __name__ == '__main__':
    main()
