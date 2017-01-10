#!/usr/bin/python3

import argparse
import re
from operator import itemgetter
from collections import defaultdict
from extra import nt_string


def read_igv(igv, outp=None):
    genes = dict()
    for gene in nt_string.read_fasta(igv):
        gene.name = re.match(r'.*(IG[^| ]*).*', gene.name).group(1)
        gene.seq = gene.seq.upper()
        genes[gene.name] = gene

        if outp:
            outp.write('>%s\n' % gene.name)
            for i in range(0, len(gene.seq), 80):
                outp.write('%s\n' % gene.seq[i:i + 80])
    return genes


def generate_novel_segment(gene, combination):
    mutations = [(int(position), nt) for position, nt in
                 [mutation.split(':') for mutation in combination.split(',')]]
    mutations.sort(key=itemgetter(0))

    new_seq = list(gene.seq)
    for position, nt in mutations:
        new_seq[position - 1] = nt
    return ''.join(new_seq)


def generate_and_write(gene, combination, i, outp):
    outp.write('>%s-M%d\n' % (gene.name, i))
    new_seq = generate_novel_segment(gene, combinations)
    for i in range(0, len(new_seq), 80):
        outp.write('%s\n' % new_seq[i:i + 80])


def generate(genes, combinations, outp):
    while next(combinations).startswith('#'):
        pass

    i = 1
    for line in combinations:
        split_line = line.strip().split('\t')
        gene = split_line[0]
        combination = split_line[3]
        generate_and_write(genes[gene], combination, i, outp)
        i += 1


def main():
    parser = argparse.ArgumentParser(description='generate new IGHV using common mutations')
    parser.add_argument('-c', '--combinations', help='combinations csv file',
                        type=argparse.FileType('r'), required=True, metavar='C')
    parser.add_argument('-i', '--igv', help='V genes fasta file', type=argparse.FileType('r'),
                        required=True, metavar='I')
    parser.add_argument('-o', '--output', help='output fasta', type=argparse.FileType('w'),
                        required=True, metavar='O')
    parser.add_argument('--discard-original', help='do not include original IGHV fasta entries',
                        action='store_true', dest='discard_original')
    args = parser.parse_args()

    genes = read_igv(args.igv, None if args.discard_original else args.output)
    generate(genes, args.combinations, args.output, gene_combinations)


if __name__ == '__main__':
    main()
