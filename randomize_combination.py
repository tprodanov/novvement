#!/usr/bin/python3

import sys
import argparse
import random

import similarity_filter
import combinations_to_segments
from extra._version import __version__

def random_combination(l, r, segment_seq, length, min_dist, source):
    while True:
        positions = sorted(random.sample(range(l, r), length))
        nts = [random.choice('ACGT') for _ in range(length)]
        comb = ','.join('%d:%s' % x for x in zip(positions, nts))
        seq = combinations_to_segments.generate_novel_segment(segment_seq, comb)[l:r]

        if similarity_filter.min_distance(seq, source)[0] >= min_dist:
            return comb


def randomize_combinations(args):
    segments = combinations_to_segments.read_segments(args.v_segments)
    inp = args.input
    outp = args.output

    l, r = args.range
    r += 1

    source = [(segment.seq[l:r], segment.name) for segment in segments.values()]

    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('segment\tsignificance\tlength\tcombination\n')

    while next(inp).startswith('#'):
        pass

    for line in inp:
        line = line.strip().split('\t')
        segment = line[0]
        length = int(line[2])
        segment_seq = segments[segment].seq
        comb = random_combination(l, min(r, len(segment_seq) + 1), segment_seq, length, args.min_dist, source)
        outp.write('%s\t0\t%d\t%s\n' % (segment, length, comb))


def main():
    parser = argparse.ArgumentParser(description='Randomize combinations',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    io_args = parser.add_argument_group('Input/output files')
    io_args.add_argument('-i', '--input', help='Input combinations csv file', type=argparse.FileType(),
                         metavar='File', required=True)
    io_args.add_argument('-v', '--v-segments', help='V segments fasta', type=argparse.FileType(),
                         metavar='File', required=True, dest='v_segments')
    io_args.add_argument('-o', '--output', help='Output combinations csv file', type=argparse.FileType('w'),
                         metavar='File', required=True)

    opt_args = parser.add_argument_group('Optional arguments')
    opt_args.add_argument('--range', help='Positions range (default: [60, 290])', metavar=('Int', 'Int'),
                          nargs=2, default=[60, 290], type=int)
    opt_args.add_argument('--min-dist', help='Min distance from any V segment (default: 2)',
                          type=int, metavar='Int', default=2, dest='min_dist')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other.add_argument('--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
    randomize_combinations(args)


if __name__ == '__main__':
    main()
