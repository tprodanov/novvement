#!/usr/bin/python3

import argparse
import sys
import itertools

import generate_possible_igv


def distance(seq1, seq2):
    return sum(nt1 != nt2 for nt1, nt2 in itertools.zip_longest(seq1, seq2))


def min_distance(seq, source):
    nearest = None
    dist = None

    for entry in source:
        seq2 = entry[0]
        d = distance(seq, seq2)

        if dist is None or d < dist:
            dist = d
            nearest = entry[1:]
    return (dist, *nearest) if dist is not None else (None, )


def keep(name, seq, significance, source, target, args):
    if args.source_dist:
        dist, gene = min_distance(seq, source)
        if dist < args.source_dist:
            args.log.write('%s is too close to %s (distance %d)\n' % (name, gene, dist))
            return False

    if not args.target_dist:
        return True
    
    if not target:
        return True

    dist, nearest_significance, oth = min_distance(seq, target)
    if dist >= args.target_dist:
        return True
    if significance * args.target_mf < nearest_significance:
        args.log.write('%s is too close to %s (distance %d, mf %d / %d = %.2f)\n'
                       % (name, oth, dist, nearest_significance, significance, nearest_significance / significance))
        return False
    return True


def run(args):
    segments = generate_possible_igv.read_igv(args.v_segments)
    inp = args.combinations
    outp = args.output
    l, r = args.range
    l -= 1

    source = [(gene.seq[l:r], gene.name) for gene in segments.values()]
    target = []

    line = next(inp)
    while line.startswith('#'):
        line = next(inp)
    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write(line)

    for i, line in enumerate(inp):
        split_line = line.strip().split('\t')
        gene = split_line[0]
        significance = int(split_line[1])
        combination = split_line[3]
        seq = generate_possible_igv.generate_novel_segment(segments[gene].seq, combination)[l:r]

        name = '%s-M%d' % (gene, i + 1)

        if keep(name, seq, significance, source, target, args):
            outp.write(line)
            target.append((seq, significance, name))


def main():
    parser = argparse.ArgumentParser(description='filter out too similar combinations', add_help=False)
    io_args = parser.add_argument_group('input/output arguments')
    io_args.add_argument('-c', '--combinations', help='input: target combinations csv file',
                         type=argparse.FileType(), required=True, metavar='FILE')
    io_args.add_argument('-v', '--v-segments', help='input: source v segments',
                         type=argparse.FileType(), required=True, metavar='FILE', dest='v_segments')
    io_args.add_argument('-o', '--output', help='output: filtered combinations csv file',
                         type=argparse.FileType('w'), metavar='FILE', required=True)
    io_args.add_argument('-l', '--log', help='output: log file (optional)',
                         type=argparse.FileType('w'), metavar='FILE', default='/dev/null')

    filter_args = parser.add_argument_group('filter arguments')
    filter_args.add_argument('--range', help='positions range (default: [60, 280])',
                             metavar=('INT', 'INT'), type=int, nargs=2, default=[60, 280])
    filter_args.add_argument('--source-dist', help='minimum required distance to '
                                                   'any source segment (default: 2)',
                             metavar='INT', dest='source_dist', type=int, default=2)
    filter_args.add_argument('--target-dist', help='minimum reliable distance to '
                                                   'any target combination (default: 2)',
                             metavar='INT', dest='target_dist', type=int, default=2)
    filter_args.add_argument('--target-mf', help='significance multiplication factor '
                                                 'to filter out combination with unreliable '
                                                 'target distance (default: 3)',
                             metavar='FLOAT', dest='target_mf', type=float, default=3)

    other = parser.add_argument_group('other arguments')
    other.add_argument('-h', '--help', action='help', help='show this help message and exit')

    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
