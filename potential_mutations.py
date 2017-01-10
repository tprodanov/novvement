#!/usr/bin/python3

import argparse
import sys
from collections import Counter
from collections import defaultdict


class Parameters:
    left_margin = 60
    right_margin = 280
    coverage_threshold = 200
    mutation_rate = 0.1


def read_v_alignment(v_alignment, gene_coverages):
    mutations = defaultdict(Counter)
    
    while next(v_alignment).startswith('#'):
            pass

    for line in v_alignment:
        gene, read, error_type, position, nt = line.strip().split('\t')
        position = int(position)
        if error_type == 'mismatch' and gene_coverages[gene] >= Parameters.coverage_threshold and \
                Parameters.left_margin <= position <= Parameters.right_margin:
            mutations[gene][(position, nt)] += 1
    return mutations


def filter_mutations(gene, gene_mutations, gene_coverages, outp):
    i = 1
    gene_important_mutations = dict()
    gene_mutation_indexes = dict()
    for pos_nt, coverage in gene_mutations.items():
        if coverage >= gene_coverages * Parameters.mutation_rate:
            outp.write('%s\t%d\t%s\t%d\t%d\n' % (gene, pos_nt[0], pos_nt[1], coverage, gene_coverages))
            gene_important_mutations[pos_nt] = i
            gene_mutation_indexes[i] = pos_nt
            i *= 2
    return gene_important_mutations, gene_mutation_indexes


def extract_mutations(mutations, gene_coverages, outp):
    outp.write('gene\tposition\tnt\tmut_coverage\tgene_coverage\n')
    for gene, gene_mutations in mutations.items():
        filter_mutations(gene, gene_mutations, gene_coverages[gene], outp)


def calculate_potential_mutations(gene_coverage, v_alignment, outp):
    mutations = read_v_alignment(v_alignment, gene_coverage)
    extract_mutations(mutations, gene_coverage, outp)


def load_potential_mutations(f):
    potential_mutations = defaultdict(Counter)
    
    while next(f).startswith('#'):
        pass

    for line in f:
        split_line = line.strip().split()
        gene = split_line[0]
        position = int(split_line[1])
        nt = split_line[2]
        mut_coverage = int(split_line[3])
        # gene_coverage = int(split_line[4])

        potential_mutations[gene][(position, nt)] = mut_coverage
    return potential_mutations


def load_gene_coverage(gene_coverage_f):
    gene_coverage = dict()
    for line in gene_coverage_f:
        gene, coverage = line.strip().split()
        gene_coverage[gene] = int(coverage)
    return gene_coverage


def main():
    parser = argparse.ArgumentParser(description='potential mutations extraction', add_help=False)
    input_files = parser.add_argument_group('input files')
    input_files.add_argument('-g', '--gene-coverage', help='file with lines <gene> <coverage>',
                             type=argparse.FileType('r'), required=True, metavar='G',
                             dest='gene_coverage')
    input_files.add_argument('-v', '--v-alignment', help='v_alignment.csv', dest='v_alignment',
                             type=argparse.FileType('r'), required=True, metavar='V')

    output_files = parser.add_argument_group('output files')
    output_files.add_argument('-o', '--output', help='potential mutations in csv format',
                              type=argparse.FileType('w'), required=True, metavar='O')

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('--range', help='positions range (default: [60, 280])',
                          metavar=('L', 'R'), nargs=2, default=[60, 280], type=int)
    optional.add_argument('--coverage', help='gene coverage threshold (default: 200)',
                          type=int, default=200)
    optional.add_argument('--rate', help='mutation rate threshold (default: 0.1)',
                          type=float, default=0.1)

    other = parser.add_argument_group('other arguments')
    other.add_argument('-h', '--help', action='help', help='show this help message and exit')
    args = parser.parse_args()

    Parameters.left_margin = args.range[0]
    Parameters.right_margin = args.range[1]
    Parameters.coverage_threshold = args.coverage
    Parameters.mutation_rate = args.rate

    gene_coverage = load_gene_coverage(args.gene_coverage)
    args.output.write('# %s\n' % ' '.join(sys.argv))
    calculate_potential_mutations(gene_coverage, args.v_alignment, args.output)


if __name__ == '__main__':
    main()
