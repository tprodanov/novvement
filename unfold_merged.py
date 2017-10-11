#!/usr/bin/env python3


def splitter(delimeter):
    def inner(seq):
        return seq.split(delimeter)
    return inner


def novel_entries(combinations, labels):
    from operator import itemgetter
    size = min(map(int, map(itemgetter(1), map(splitter('|'), combinations.split(';')))))
    labels = list(map(int, map(itemgetter(1), map(splitter(':'), labels.split(',')))))
    for label in labels:
        yield '%d\t%d' % (size, label)

def hamming_distance(seq1, seq2):
    from itertools import zip_longest
    return sum(nt1 != nt2 for nt1, nt2 in zip_longest(seq1, seq2))

def write_component(component_num, lines, outp, outp_edges,
                    min_significance, left_margin, right_margin, max_distance):
    seqs = []
    novel = []
    for line in lines:
        significance, combinations, seq, labels = line.strip().split('\t')
        if int(significance) < min_significance:
            continue
        seqs.append(seq[left_margin : right_margin])
        novel.append((combinations, labels))

    for i, seq1 in enumerate(seqs):
        for j in range(i + 1, len(seqs)):
            seq2 = seqs[j]
            if hamming_distance(seq1, seq2) <= max_distance:
                outp_edges.write('%d\t%d\t%d\t1\n' % (component_num, i, j))

    for i, (combinations, labels) in enumerate(novel):
        for novel_str in novel_entries(combinations, labels):
            outp.write('%d\t%d\t%s\n' % (component_num, i, novel_str))


def analyze_file(f, outp, outp_edges,
                 **kwargs):
    outp.write('component\tnovel\tsize\tlabel\n')
    outp_edges.write('component\tx1\tx2\ty\n')
    while next(f).startswith('#'):
        pass
    
    current_lines = []
    current_number = -1
    
    for line in f:
        if line.startswith('*'):
            if current_lines:
                write_component(current_number, current_lines, outp, outp_edges,
                                **kwargs)
            current_lines = []
            current_number += 1
        else:
            current_lines.append(line)
    if current_lines:
        write_component(current_number, current_lines, outp, outp_edges, **kwargs)



def main():
    import argparse
    from extra._version import __version__

    parser = argparse.ArgumentParser('Unfold merged_full csv',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -i File -o File -e File [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', help='Input merged_full.csv file',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-o', '--output', help='Output csv file',
                         type=argparse.FileType('w'), required=True, metavar='File')
    io_args.add_argument('-e', '--edges', help='Output csv file with edges between sequences',
                         type=argparse.FileType('w'), required=True, metavar='File')

    alg_args = parser.add_argument_group('Novvement arguments')
    alg_args.add_argument('-s', '--significance', type=int, metavar='Int', default=1,
                          help='Minimal significance (default: %(default)s)')
    alg_args.add_argument('-r', '--range', type=int, metavar=('Int', 'Int'), nargs=2, default=[40, 290],
                          help='Positions range (default: %(default)s)')
    alg_args.add_argument('-d', '--distance', type=int, metavar='Int', default=1,
                          help='Hamming distance between sequences (default: %(default)s)')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
    analyze_file(args.input, args.output, args.edges,
                 min_significance=args.significance,
                 left_margin=args.range[0] - 1,
                 right_margin=args.range[1],
                 max_distance=args.distance)


if __name__ == '__main__':
    main()

