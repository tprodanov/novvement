#!/usr/bin/python3

import argparse
import re
import sys

from extra import nt_string


def analyze_v_alignment(inp, outp):
    fasta_iter = nt_string.read_fasta(inp)
    outp.write('segment\tread\ttype\tposition\tnt\n')

    try:
        while True:
            read = next(fasta_iter)
            segment = next(fasta_iter)

            seq1 = read.seq
            seq2 = segment.seq
            read_name = re.search(r'READ:([^|\n]*)', read.name).group(1)
            segment_name = re.search(r'GENE:([^|\n]*)', segment.name).group(1)

            assert (len(seq1) == len(seq2))
            j = 1
            for i in range(len(seq1)):
                if seq1[i] != seq2[i]:
                    if seq1[i] == '-':
                        outp.write('%s\t%s\tdeletion\t%d\t-\n' % (segment_name, read_name, j))
                    elif seq2[i] == '-':
                        outp.write('%s\t%s\tinsert\t%d\t%s\n' % (segment_name, read_name, j, seq1[i]))
                        j -= 1
                    elif seq1[i] == 'N' or seq2[i] == 'N':
                        outp.write('%s\t%s\degenerate\t%d\t%s\n' % (segment_name, read_name, j, seq1[i]))
                    else:
                        outp.write('%s\t%s\tmismatch\t%d\t%s\n' % (segment_name, read_name, j, seq1[i]))
                j += 1
    except StopIteration:
        pass


def main():
    parser = argparse.ArgumentParser(description='v alignment mismatches and indels statistics')
    parser.add_argument('-i', '--input', help='input: v_alignment.fasta',
                        type=argparse.FileType('r'), required=True, metavar='FILE')
    parser.add_argument('-o', '--output', help='output: v_alignment.csv',
                        type=argparse.FileType('w'), required=True, metavar='FILE')

    args = parser.parse_args()
    args.output.write('# %s\n' % ' '.join(sys.argv))
    analyze_v_alignment(args.input, args.output)


if __name__ == '__main__':
    main()
