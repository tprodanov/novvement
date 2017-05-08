#!/usr/bin/python3

import argparse
import re
from operator import itemgetter
from collections import defaultdict
from extra import nt_string
from extra._version import __version__


def read_segments(inp, outp=None):
    segments = dict()
    for segment in nt_string.read_fasta(inp):
        segment.name = segment.name.lstrip('>')
        segments[segment.name] = segment

        if outp:
            outp.write('>%s\n' % segment.name)
            for i in range(0, len(segment.seq), 80):
                outp.write('%s\n' % segment.seq[i:i + 80])
    return segments


def generate_novel_segment(segment_seq, combination):
    mismatches = [(int(position), nt) for position, nt in
                 [mismatch.split(':') for mismatch in combination.split(',')]]
    mismatches.sort(key=itemgetter(0))

    new_seq = list(segment_seq)
    for position, nt in mismatches:
        if nt == '-':
            new_seq = new_seq[:(position -1)]
            break
        new_seq[position - 1] = nt
    return ''.join(new_seq)


def generate_and_write(segment, combination, i, outp):
    outp.write('>%s-M%d\n' % (segment.name, i))
    new_seq = generate_novel_segment(segment.seq, combination)
    for i in range(0, len(new_seq), 80):
        outp.write('%s\n' % new_seq[i:i + 80])


def generate(segments, combinations, outp):
    while next(combinations).startswith('#'):
        pass

    i = 1
    for line in combinations:
        split_line = line.strip().split('\t')
        segment = split_line[0]
        combination = split_line[1]
        generate_and_write(segments[segment], combination, i, outp)
        i += 1


def main():
    parser = argparse.ArgumentParser(description='Convert combinations into sequences',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', help='Input file with lines <segment> <mismatches>',
                         type=argparse.FileType('r'), required=True, metavar='File')
    io_args.add_argument('-v', '--v-segments', help='V segments fasta file', type=argparse.FileType('r'),
                         dest='v_segments', required=True, metavar='File')
    io_args.add_argument('-o', '--output', help='Output fasta', type=argparse.FileType('w'),
                         required=True, metavar='File')
    io_args.add_argument('--discard-original', help='Do not include original V segments',
                         action='store_true', dest='discard_original')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other.add_argument('--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    segments = read_segments(args.v_segments, None if args.discard_original else args.output)
    generate(segments, args.input, args.output)


if __name__ == '__main__':
    main()
