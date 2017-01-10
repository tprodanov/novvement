#!/usr/bin/python3

import sys
import argparse
import random


def random_combination(args, segment, length, generated):
    while True:
        positions = sorted(random.sample(range(args.range[0], args.range[1] + 1), length))
        nts = [random.choice('ACGT') for _ in range(length)]
        comb = ','.join('%d:%s' % x for x in zip(positions, nts))

        if (segment, comb) not in generated:
            generated.add((segment, comb))
            return comb


def shuffle_combinations(args):
    inp = args.input
    outp = args.output

    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('segment\tsignificance\tlength\tcombinations\tj_hit\tdatasets\n')

    while next(inp).startswith('#'):
        pass

    generated = set()
    for line in inp:
        line = line.strip().split('\t')
        segment = line[0]
        length = int(line[2])
        outp.write('%s\t0\t%d\t%s\t*\t*\n' % (segment, length, random_combination(args, segment, length, generated)))


def main():
    parser = argparse.ArgumentParser(description='generate random combinations')
    parser.add_argument('-i', '--input', help='input combinations csv file', type=argparse.FileType(),
                        metavar='FILE', required=True)
    parser.add_argument('-o', '--output', help='output combinations csv file', type=argparse.FileType('w'),
                        metavar='FILE', required=True)
    parser.add_argument('--range', help='positions range (default: [60, 290])', metavar=('INT', 'INT'),
                        nargs=2, default=[60, 290], type=int)
    args = parser.parse_args()

    shuffle_combinations(args)


if __name__ == '__main__':
    main()
