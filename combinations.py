#!/usr/bin/python3

import argparse
from collections import Counter
from collections import defaultdict
from operator import itemgetter
from extra import utilities
import potential_mutations
import os
import re
import sys


def index_mutations(mutations):
    mutations_list = dict()
    mutations_ind = dict()
    for gene, gene_mutations in mutations.items():
        mutations_list[gene] = list(sorted(gene_mutations, key=itemgetter(0)))
        mutations_ind[gene] = dict()
        for i in range(len(mutations_list[gene])):
            mutations_ind[gene][mutations_list[gene][i]] = i
    return mutations_list, mutations_ind


def combination_hash(combination, mutations_ind):
    h = 0
    for mut in combination:
        h += 2 ** (mutations_ind[mut])
    return h


class Combination:
    def __init__(self, combination):
        self.combination = combination
        self.coverage = 0
        self.j_hit = Counter()
        self.datasets = Counter()


class Combinations:
    def __init__(self, v_alignment, mutations, j_hit):
        # self.mutations, self.mutations_ind = index_mutations(mutations)
        self.mutations = mutations
        self.mutations_ind = index_mutations(mutations)[1]

        self.mut_combinations = defaultdict(dict)
        for read, gene, read_mutations in utilities.v_alignment_to_mutations(v_alignment,
                                                                             self.mutations):
            h = combination_hash(read_mutations, self.mutations_ind[gene])
            if h not in self.mut_combinations[gene]:
                self.mut_combinations[gene][h] = Combination(read_mutations)

            self.mut_combinations[gene][h].coverage += 1
            self.mut_combinations[gene][h].j_hit[j_hit[read]] += 1

    def write(self, outp, length, cov1, covM, validation):
        outp.write('gene\tcoverage\tlength\tcombination\tj_hit')
        if validation:
            outp.write('\tdatasets')
        outp.write('\n')

        for gene, gene_combinations in self.mut_combinations.items():
            for combination in gene_combinations.values():
                if len(combination.combination) < length:
                    continue
                if len(combination.j_hit) == 1 and combination.coverage < cov1:
                    continue
                if len(combination.j_hit) > 1 and combination.coverage < covM:
                    continue
                outp.write('%s\t%d\t%d\t%s\t%s' %
                           (gene, combination.coverage, len(combination.combination),
                            ','.join('%d:%s' % mut for mut in combination.combination),
                            ','.join('%s:%d' % j for j in combination.j_hit.items())))
                if validation:
                    outp.write('\t')
                    outp.write(','.join('%s:%d' % d for d in combination.datasets.items()))
                outp.write('\n')

    def write_human(self, outp, length, cov1, covM):
        for gene, gene_combinations in self.mut_combinations.items():
            outp.write('-------\n%s\n' % gene)
            for combination in gene_combinations.values():
                if len(combination.combination) < length:
                    continue
                if len(combination.j_hit) == 1 and combination.coverage < cov1:
                    continue
                if len(combination.j_hit) > 1 and combination.coverage < covM:
                    continue

                outp.write('\t%s:\n' % ', '.join('%3d %s' % mut for mut in combination.combination))
                outp.write('\t\tCoverage: %d\n' % combination.coverage)
                outp.write('\t\tLength: %d\n' % len(combination.combination))
                if combination.j_hit:
                    outp.write('\t\tJ hits:\n\t\t\t')
                    outp.write('\n\t\t\t'.join('%s: %d' % j for j in combination.j_hit.items()))
                    outp.write('\n')
                if combination.datasets:
                    outp.write('\t\tDatasets:\n\t\t\t')
                    outp.write('\n\t\t\t'.join('%s: %d' % d for d in combination.datasets.items()))
                    outp.write('\n')

    def validate_v_alignment(self, v_alignment, name):
        for _, gene, read_mutations in utilities.v_alignment_to_mutations(v_alignment, self.mutations):
            h = combination_hash(read_mutations, self.mutations_ind[gene])
            if h not in self.mut_combinations[gene]:
                continue

            self.mut_combinations[gene][h].datasets[name] += 1


def load_j_gene(cdr_details):
    j_hit = dict()
    for line in cdr_details:
        read_name, j_gene = line.strip().split('\t')
        j_hit[read_name] = j_gene
    return j_hit


def validate(mut_combinations, validation_f, my_v_alignment):
    for line in validation_f:
        line = line.strip()
        if line[0] == '#':
            continue
        name, v_alignment = re.findall(r"[^ \t:,]+", line)
        full_v_alignment = os.path.dirname(validation_f.name)
        if full_v_alignment:
            full_v_alignment += '/'
        full_v_alignment += v_alignment

        if os.path.abspath(full_v_alignment) == os.path.abspath(my_v_alignment):
            continue
        # print(full_v_alignment)

        with open(full_v_alignment) as v_alignment_f:
            mut_combinations.validate_v_alignment(v_alignment_f, name)


def main():
    parser = argparse.ArgumentParser(description='suspicious mutation combinations', add_help=False)
    input_files = parser.add_argument_group('input files')
    input_files.add_argument('-v', '--v-alignment', help='v_alignment.csv', dest='v_alignment',
                             type=argparse.FileType('r'), required=True, metavar='V')
    input_files.add_argument('-j', '--j-hit', help='file with lines <read name> <j gene>',
                             dest='j_hit', type=argparse.FileType('r'), required=True, metavar='J')
    input_files.add_argument('-m', '--mutations', help='potential_mutations.csv',
                             type=argparse.FileType('r'), required=True, metavar='M')
    input_files.add_argument('-d', '--datasets', type=argparse.FileType('r'),
                             help='file with lines <name> <v_alignment.csv> '
                                  'for the validation (optional)', metavar='D')

    output_files = parser.add_argument_group('output files')
    output_files.add_argument('-o', '--output', help='csv combinations output',
                              type=argparse.FileType('w'), required=True, metavar='O')
    output_files.add_argument('-u', '--human', help='human readable combinations output (optional)',
                              type=argparse.FileType('w'), metavar='H')

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('--length', help='min combination length (default: 2)',
                          default=2, type=int, metavar='INT')
    optional.add_argument('--cov-single-j', help='min coverage with the single J hit (default: 15)',
                          default=15, type=int, metavar='INT', dest='cov1')
    optional.add_argument('--cov-mult-j', help='min coverage with multiple J hits (default: 5)',
                          default=5, type=int, metavar='INT', dest='cov')

    other = parser.add_argument_group('other arguments')
    other.add_argument('-h', '--help', action='help', help='show this help message and exit')

    args = parser.parse_args()

    mutations = potential_mutations.load_potential_mutations(args.mutations)
    j_hit = load_j_gene(args.j_hit)
    mut_combinations = Combinations(args.v_alignment, mutations, j_hit)

    if args.datasets:
        validate(mut_combinations, args.datasets, args.v_alignment.name)

    args.output.write('# %s\n' % ' '.join(sys.argv))
    mut_combinations.write(args.output, args.length, args.cov, args.cov, bool(args.datasets))
    if args.human:
        mut_combinations.write_human(args.human, args.length, args.cov1, args.cov)


if __name__ == '__main__':
    main()
