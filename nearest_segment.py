#!/usr/bin/python3

import sys
import argparse
from extra import nt_string


def compare(seq1, seq2, pos_range):
    l1 = len(seq1)
    l2 = len(seq2)
    min_l = min(l1, l2)
    max_l = max(l1, l2)
    
    l, r = pos_range
    if min_l <= l:
        return r - l + 1

    s = sum(seq1[i] != seq2[i] for i in range(l, min(min_l, r + 1)))
    if l1 <= l2:
        return s
    else:
        return s + min(r + 1, l1) - min(r + 1, l2)


def find_nearest(seq, source, pos_range):
    best = None
    best_diff = pos_range[1] - pos_range[0] + 2
    for read in source:
        diff = compare(seq, read.seq, pos_range)
        if diff < best_diff:
            best = read.name
            best_diff = diff
    return best_diff, best


def analyze(target, source, pos_range, outp):
    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('target\tlen\tnearest\tdiff\n')
    for read in target:
        best_diff, best = find_nearest(read.seq, source, pos_range)
        outp.write('%s\t%d\t%s\t%d\n' % (read.name[1:], len(read.seq), best[1:], best_diff))


def main():
    parser = argparse.ArgumentParser(description='for every target segment find nearest source segment')
    parser.add_argument('-t', '--target', help='fasta file with segments to analyze',
                        type=argparse.FileType(), required=True, metavar='FILE')
    parser.add_argument('-s', '--source', help='fasta file with source segments',
                        type=argparse.FileType(), required=True, metavar='FILE')
    parser.add_argument('-o', '--output', help='csv output file', type=argparse.FileType('w'),
                        required=True, metavar='FILE')
    parser.add_argument('--range', help='positions range (default: [60, 280])',
                        metavar=('INT', 'INT'), nargs=2, default=[60, 280], type=int)
    args = parser.parse_args()

    target = list(nt_string.read_fasta(args.target))
    source = list(nt_string.read_fasta(args.source))
    analyze(target, source, args.range, args.output)



if __name__ == '__main__':
    main()

