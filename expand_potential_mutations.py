#!/usr/bin/python3

import argparse
import sys
from collections import Counter
from collections import defaultdict
from operator import itemgetter
import potential_mutations
import combinations


class Parameters:
    rate = 0.9
    left_margin = 10
    right_margin = 290
    combination_coverage = 50


class Combination:
    def __init__(self, line, mutations_ind):
        split_line = line.strip().split('\t')
        self.gene = split_line[0]
        self.coverage = int(split_line[1])

        combination = [(int(position), nt) for position, nt in
                       [mutation.split(':') for mutation in split_line[3].split(',')]]
        self.hash = combinations.combination_hash(combination, mutations_ind[self.gene])
        self.initial = set(combination)
        self.new_mutations = Counter()


def load_combinations(combinations_file, mutations_ind):
    combinations_dict = defaultdict(dict)
    while next(combinations_file).startswith('#'):
        pass

    for line in combinations_file:
        combination = Combination(line, mutations_ind)
        if combination.coverage < Parameters.combination_coverage:
            continue
        combinations_dict[combination.gene][combination.hash] = combination
    return combinations_dict


def process_read(read_mutations, gene_combinations, gene_mutations, gene_mutations_ind):
    potential = []
    for pos, nt in read_mutations:
        if (pos, nt) in gene_mutations:
            potential.append((pos, nt))
    potential.sort(key=itemgetter(0))
    h = combinations.combination_hash(potential, gene_mutations_ind)

    if h not in gene_combinations:
        return
    read_combination = gene_combinations[h]
    for pos, nt in read_mutations:
        if (pos, nt) not in read_combination.initial:
            read_combination.new_mutations[(pos, nt)] += 1


def process_v_alignment(v_alignment, combinations_dict, mutations, mutations_ind):
    mutations_coverage = defaultdict(Counter)
    genes_coverage = Counter()

    while next(v_alignment).startswith('#'):
        pass

    prev_read = None
    prev_gene = None
    read_mutations = None

    for line in v_alignment:
        gene, read, mut_type, position, nt = line.strip().split()
        position = int(position)

        if prev_read != read:
            if read_mutations and prev_gene in combinations_dict:
                process_read(read_mutations, combinations_dict[prev_gene],
                             mutations[prev_gene], mutations_ind[prev_gene])
            prev_read = read
            prev_gene = gene
            genes_coverage[prev_gene] += 1
            read_mutations = []

        if mut_type == 'mismatch' and Parameters.left_margin <= position <= Parameters.right_margin:
            mutations_coverage[gene][(position, nt)] += 1
            read_mutations.append((position, nt))

    if read_mutations and prev_gene in combinations_dict:
        process_read(read_mutations, combinations_dict[prev_gene],
                     mutations[prev_gene], mutations_ind[prev_gene])
    genes_coverage[prev_gene] += 1
    return mutations_coverage, genes_coverage


def print_mutations(outp, expanded_output, mutations, combinations_dict,
                    mutations_coverage, genes_coverage):
    outp.write('gene\tposition\tnt\tmut_coverage\tgene_coverage\n')
    if expanded_output:
        expanded_output.write('gene\tposition\tnt\t'
                              'cov_in_comb\tcomb_coverage\tmut_coverage\tgene_coverage\n')

    for gene, gene_coverage in genes_coverage.items():
        if gene not in mutations:
            continue
        gene_mutations = set(mutations[gene])

        for combination in combinations_dict[gene].values():
            for new_mut, cov in combination.new_mutations.items():
                if cov >= Parameters.rate * combination.coverage:
                    gene_mutations.add(new_mut)
                    if expanded_output:
                        expanded_output.write('%s\t%d\t%s\t%d\t%d\t%d\t%d\n' %
                                              (gene, *new_mut, cov, combination.coverage,
                                               mutations_coverage[gene][new_mut], gene_coverage))

        for mut in gene_mutations:
            outp.write('%s\t%d\t%s\t%d\t%d\n' % (gene, *mut, mutations_coverage[gene][mut],
                                                 gene_coverage))


def main():
    parser = argparse.ArgumentParser(description='expand potential mutations', add_help=False)
    input_files = parser.add_argument_group('input files')
    input_files.add_argument('-v', '--v-alignment', help='v_alignment.csv', dest='v_alignment',
                             type=argparse.FileType('r'), required=True, metavar='V')
    input_files.add_argument('-m', '--mutations', help='potential_mutations.csv',
                             type=argparse.FileType('r'), required=True, metavar='M')
    input_files.add_argument('-c', '--combinations', help='combinations.csv',
                             type=argparse.FileType('r'), required=True, metavar='C')

    output_files = parser.add_argument_group('output files')
    output_files.add_argument('-o', '--output', help='expanded potential mutations in csv format',
                              type=argparse.FileType('w'), required=True, metavar='O')
    output_files.add_argument('-e', '--expanded', help='only expanded mutations (optional)',
                              type=argparse.FileType('w'), metavar='E')

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('--range', help='positions range (default: [60, 280])',
                          metavar=('L', 'R'), nargs=2, default=[60, 280])
    optional.add_argument('--coverage', help='combination coverage threshold (default: 50)',
                          type=int, default=50)
    optional.add_argument('--rate', help='mutation coverage threshold (ratio to combination coverage)'
                                          ' (default: 0.9)',
                          type=float, default=0.9)

    other = parser.add_argument_group('other arguments')
    other.add_argument('-h', '--help', action='help', help='show this help message and exit')

    args = parser.parse_args()

    Parameters.left_margin = int(args.range[0])
    Parameters.right_margin = int(args.range[1])
    Parameters.combination_coverage = args.coverage
    Parameters.rate = args.rate

    mutations = potential_mutations.load_potential_mutations(args.mutations)
    mutations_ind = combinations.index_mutations(mutations)[1]

    combinations_dict = load_combinations(args.combinations, mutations_ind)
    mutations_coverage, genes_coverage = process_v_alignment(args.v_alignment, combinations_dict,
                                                             mutations, mutations_ind)
    args.output.write('# %s\n' % ' '.join(sys.argv))
    if args.expanded:
        args.expanded.write('# %s\n' % ' '.join(sys.argv))
    print_mutations(args.output, args.expanded, mutations, combinations_dict,
                    mutations_coverage, genes_coverage)


if __name__ == '__main__':
    main()
