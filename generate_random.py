#!/usr/bin/python3

import sys
import argparse
import random


def random_combination(args, gene, length, generated):
    while True:
        positions = sorted(random.sample(range(args.range[0], args.range[1] + 1), length))
        nts = [random.choice('ACGT') for _ in range(length)]
        comb = ','.join('%d:%s' % x for x in zip(positions, nts))
    
        if (gene, comb) not in generated:
            generated.add((gene, comb))
            return comb


def shuffle_combinations(args):
    inp = args.input
    outp = args.output

    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('gene\tsignificance\tlength\tcombinations\tj_hit\tdatasets\n')

    while next(inp).startswith('#'):
        pass

    generated = set()
    for line in inp:
        line = line.strip().split('\t')
        gene = line[0]
        length = int(line[2])
        outp.write('%s\t*\t%d\t%s\t*\t*\n' % (gene, length, random_combination(args, gene, length, generated)))


def main():
    parser = argparse.ArgumentParser(description='generate random combinations')
    parser.add_argument('-i', '--input', help='input combinations csv file', type=argparse.FileType(), metavar='FILE', required=True)
    parser.add_argument('-o', '--output', help='output combinations csv file', type=argparse.FileType('w'), metavar='FILE', required=True)
    parser.add_argument('--range', help='positions range (default: [60, 280])', metavar=('INT', 'INT'),
                        nargs=2, default=[60, 280], type=int)
    args = parser.parse_args()

    shuffle_combinations(args)


if __name__ == '__main__':
    main()

