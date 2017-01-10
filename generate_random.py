#!/usr/bin/python3

import sys
import argparse
import random

import similarity_filter
import combinations_to_segments

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
    l -= 1

    source = [(segment.seq[l:r], segment.name) for segment in segments.values()]

    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('segment\tsignificance\tlength\tcombinations\tj_hit\tdatasets\n')

    while next(inp).startswith('#'):
        pass

    for line in inp:
        line = line.strip().split('\t')
        segment = line[0]
        length = int(line[2])
        comb = random_combination(l, r, segments[segment].seq, length, args.min_dist, source)
        outp.write('%s\t0\t%d\t%s\t*\t*\n' % (segment, length, comb))


def main():
    parser = argparse.ArgumentParser(description='generate random combinations')
    parser.add_argument('-i', '--input', help='input combinations csv file', type=argparse.FileType(),
                        metavar='FILE', required=True)
    parser.add_argument('-v', '--v-segments', help='V segments fasta', type=argparse.FileType(),
                        metavar='FILE', required=True, dest='v_segments')
    parser.add_argument('-o', '--output', help='output combinations csv file', type=argparse.FileType('w'),
                        metavar='FILE', required=True)
    parser.add_argument('--range', help='positions range (default: [60, 290])', metavar=('INT', 'INT'),
                        nargs=2, default=[60, 290], type=int)
    parser.add_argument('--min-dist', help='min distance from any V segment (default: 2)',
                        type=int, metavar='INT', default=2, dest='min_dist')
    args = parser.parse_args()

    randomize_combinations(args)


if __name__ == '__main__':
    main()
