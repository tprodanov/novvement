#!/usr/bin/env python3


import collections


def reformat_consensus(consensus):
    if consensus == '*':
        return consensus
    consensus = [s.split(':') for s in consensus.split(',')[::-1]]
    result = []
    inserts = collections.defaultdict(list)

    consensus = iter(consensus)
    for pos, alt in consensus:
        pos = int(pos)

        deletion = 0
        try:
            while alt == '-':
                deletion += 1
                new_pos, alt = next(consensus)
        except StopIteration:
            new_pos = -1

        if deletion > 0:
            result.append((pos, 'd%d' % deletion))
            if new_pos == -1:
                break
            pos = int(new_pos)

        if alt.lower() == alt:
            inserts[pos].append(alt)
        else:
            result.append((pos, alt))

    for pos, alt in inserts.items():
        resutl.append((pos, ''.join(alt).upper()))
    result.sort()
    return ','.join('%d:%s' % item for item in result)


class FitResult:
    def __init__(self, line):
        self.subset, self.segment, score, mismatches, indels, consensus = line.strip().split('\t')
        self.score = float(score)
        self.mismatches = int(mismatches)
        self.indels = int(indels)
        self.consensus = reformat_consensus(consensus)


class SubsetSummary:
    def __init__(self, line):
        self.subset, self.coverage, self.labels, self.consensus = line.strip().split('\t')
        self.coverage = int(self.coverage)
        self.labels = int(self.labels)


def process_fit(f, summaries, full_seqs, cropped_seqs, outp, *, differences, indels_enough):
    import sys
    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('subset\tsegment\tcoverage\tlabels\tscore\tmismatches\tindels\tconsensus\tfull_seq\tcropped_seq\n')

    while next(f).startswith('#'):
        pass
    for line in f:
        fit = FitResult(line)
        if fit.indels < indels_enough and fit.indels + fit.mismatches < differences:
            continue

        summary = summaries[fit.subset]
        outp.write('{subset}\t{segment}\t{coverage}\t{labels}\t{score:.1f}\t{mismatches}\t{indels}\t{consensus}\t'
                   '{full_seq}\t{cropped_seq}\n'
                   .format(subset=fit.subset, segment=fit.segment, coverage=summary.coverage, labels=summary.labels,
                           score=fit.score, mismatches=fit.mismatches, indels=fit.indels,
                           consensus=fit.consensus, full_seq=full_seqs[fit.subset], cropped_seq=cropped_seqs[fit.subset]))


def load_summaries(f):
    while next(f).startswith('#'):
        pass
    summaries = [SubsetSummary(line) for line in f]
    return {summary.subset: summary for summary in summaries}


def load_sequences(f):
    from extra import nt_string
    return {s.name: s.seq for s in nt_string.read_fasta(f)}


def main():
    import argparse
    from extra._version import __version__

    parser = argparse.ArgumentParser(description='Create summary over novel segments',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -s File -f File -S File File -o File [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-s', '--subsets', type=argparse.FileType(), metavar='File',
                         help='Input csv file with summary over subsets')
    io_args.add_argument('-f', '--fit', type=argparse.FileType(), metavar='File',
                         help='Input csv file with the results of the fitting alignment')
    io_args.add_argument('-S', '--sequences', type=argparse.FileType(), metavar=('File', 'File'), nargs=2,
                         help='Two input fasta files, with full and cropped sequences')
    io_args.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='File',
                         help='Output csv file')

    filt_args = parser.add_argument_group('Filtering arguments')
    filt_args.add_argument('-D', '--differences', type=int, metavar='Int', default=2,
                           help='Minimal required number of mismatches + indels (default: %(default)s)')
    filt_args.add_argument('-I', '--indels-enough', type=int, metavar='Int', default=1,
                           help='Number of indels enough for selecting a subset (default: %(default)s)')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    full_seqs = load_sequences(args.sequences[0])
    cropped_seqs = load_sequences(args.sequences[1])
    summaries = load_summaries(args.subsets)
    process_fit(args.fit, summaries, full_seqs, cropped_seqs, args.output,
                differences=args.differences, indels_enough=args.indels_enough)
 

if __name__ == '__main__':
    main()
