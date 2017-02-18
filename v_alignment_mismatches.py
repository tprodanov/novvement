#!/usr/bin/python3

import argparse
import re
import sys

from extra import nt_string
from extra._version import __version__


def polymorphism_type(nt1, nt2, j):
    if nt1 == '-':
        return 'deletion', j
    if nt2 == '-':
        return 'insert', j - 1
    if nt1 == 'N' or nt2 == 'N':
        return 'degenerate', j
    return 'mismatch', j



def analyze_v_alignment(inp, outp, mismatches_only):
    outp.write('segment\tread\ttype\tposition\tnt\n')
    if mismatches_only:
        not_mism_write = lambda s: None
    else:
        not_mism_write = lambda s: outp.write(s)

    fasta_iter = nt_string.read_fasta(inp)
    for read in fasta_iter:
        segment = next(fasta_iter)

        seq1 = read.seq
        seq2 = segment.seq
        read_name = re.search(r'READ:([^|\n]*)', read.name).group(1)
        segment_name = re.search(r'GENE:([^|\n]*)', segment.name).group(1)

        assert (len(seq1) == len(seq2))
        j = 0
        for i, (nt1, nt2) in enumerate(zip(seq1, seq2)):
            j += 1
            if nt1 == nt2:
                continue
            polym, j = polymorphism_type(nt1, nt2, j)
            if mismatches_only and polym[0] != 'm':
                continue

            outp.write('%s\t%s\t%s\t%d\t%s\n' % (segment_name, read_name, polym, j, nt1))


def main():
    parser = argparse.ArgumentParser(description='V alignment mismatches and indels statistics',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', help='Input: v_alignment.fasta',
                         type=argparse.FileType('r'), required=True, metavar='File')
    io_args.add_argument('-o', '--output', help='Output: v_alignment.csv',
                         type=argparse.FileType('w'), required=True, metavar='File')
    io_args.add_argument('-m', '--mismatches-only', help='Print mismatches only',
                         action='store_true', dest='mismatches_only')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other.add_argument('--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
    args.output.write('# %s\n' % ' '.join(sys.argv))
    analyze_v_alignment(args.input, args.output, args.mismatches_only)


if __name__ == '__main__':
    main()
