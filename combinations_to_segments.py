#!/usr/bin/python3

import argparse
import re
from operator import itemgetter
from collections import defaultdict
from extra import nt_string


def read_segments(segments, outp=None):
    segments = dict()
    for segment in nt_string.read_fasta(segments):
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
        new_seq[position - 1] = nt
    return ''.join(new_seq)


def generate_and_write(segment, combination, i, outp):
    outp.write('>%s-M%d\n' % (segment.name, i))
    new_seq = generate_novel_segment(segment.seq, combinations)
    for i in range(0, len(new_seq), 80):
        outp.write('%s\n' % new_seq[i:i + 80])


def generate(segments, combinations, outp):
    while next(combinations).startswith('#'):
        pass

    i = 1
    for line in combinations:
        split_line = line.strip().split('\t')
        segment = split_line[0]
        combination = split_line[3]
        generate_and_write(segments[segment], combination, i, outp)
        i += 1


def main():
    parser = argparse.ArgumentParser(description='convert combinations into sequences')
    parser.add_argument('-c', '--combinations', help='combinations csv file',
                        type=argparse.FileType('r'), required=True, metavar='FILE')
    parser.add_argument('-v', '--v-segments', help='V segments fasta file', type=argparse.FileType('r'),
                        dest='v_segments', required=True, metavar='FILE')
    parser.add_argument('-o', '--output', help='output fasta', type=argparse.FileType('w'),
                        required=True, metavar='FILE')
    parser.add_argument('--discard-original', help='do not include original IGHV fasta entries',
                        action='store_true', dest='discard_original')
    args = parser.parse_args()

    segments = read_segments(args.v_segments, None if args.discard_original else args.output)
    generate(segments, args.combinations, args.output, segment_combinations)


if __name__ == '__main__':
    main()
