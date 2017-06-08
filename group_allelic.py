#!/usr/bin/python3

import argparse
import sys
import os
import itertools
from operator import itemgetter, attrgetter
from collections import Counter, deque
from tqdm import tqdm

import combinations_to_segments
from extra._version import __version__


def distance(seq1, seq2):
    return sum(nt1 != nt2 for nt1, nt2 in zip(seq1, seq2)) + 1.2 * (abs(len(seq1) - len(seq2)))


def min_distance(seq, source):
    best = None
    dist = None

    for entry in source:
        seq2 = entry[0]
        d = distance(seq, seq2)

        if dist is None or d < dist:
            dist = d
            nearest = entry[1:]
    return (dist, *nearest) if dist is not None else None


def letter(number):
    from string import ascii_letters
    m = len(ascii_letters)
    if number < m:
        return ascii_letters[number]
    return ascii_letters[number // m - 1] + ascii_letters[number % m]


class UnknownSegment(Exception):
    pass


def count_covered(it, threshold):
    return sum(count >= threshold for count in it)


class SegmentInDataset:
    def __init__(self, name, line=None):
        self.name = name

        if line is None:
            self.coverage = 0
            self.j_hits = Counter()
            self.covered_j = 0
            self.seq = ''
            self.original = []
            return

        segment, coverage, length, combination, j_hits = line.strip().split('\t')
        
        if segment not in NovelSegment.segments:
            raise UnknownSegment('Segment %s unknown' % segment)

        self.name = name
        self.coverage = int(coverage)
        j_hits = [j_hit.split(':') for j_hit in j_hits.split(',')]

        self.j_hits = Counter()
        for j_segm, count in j_hits:
            self.j_hits[j_segm] = int(count)

        self.covered_j = count_covered(self.j_hits.values(), NovelSegment.j_coverage)

        self.seq = combinations_to_segments.generate_novel_segment(NovelSegment.segments[segment].seq,
                                                                   combination)
        self.seq = self.seq[NovelSegment.left_margin : NovelSegment.right_margin]

        self.original = [(segment, combination)]

    def distance(self, other):
        return distance(self.seq, other.seq)

    def merge(self, other):
        assert isinstance(other, SegmentInDataset)
        assert self.name == other.name
        
        new = SegmentInDataset(self.name)
        new.seq = self.seq
        new.coverage = self.coverage + other.coverage
        new.j_hits = self.j_hits + other.j_hits
        new.covered_j = count_covered(new.j_hits.values(), NovelSegment.j_coverage)
        new.original = self.original + other.original
        return new

    def decoration(self):
        return self.covered_j, self.coverage

    def __str__(self):
        return '%s:%d:%d' % (self.name, self.covered_j, self.coverage)

    def full_info(self, include_original):
        res = '%-8s J:%d, Cov:%d  |  ' % (self.name, self.covered_j, self.coverage)
        res += ' '.join('%s:%d' % t for t in self.j_hits.most_common())
        if include_original:
            res += '  |  '
            res += ' & '.join('%s.%s' % t for t in self.original)
        return res


class NovelSegment:
    left_margin = None
    right_margin = None
    j_coverage = None
    segments = None
    source = None
    nicknames = Counter()
    hamming_threshold = None
    shared_polymorphism_rate = None
    show_neighbors = True
    
    def __init__(self, other):
        assert isinstance(other, NovelSegment) \
            or isinstance(other, SegmentInDataset)

        if isinstance(other, NovelSegment):
            self.datasets = other.datasets.copy()
            self.seq = other.seq
            self.neighbors = None
            self.nearest_dist = other.nearest_dist
            self.nearest_segment = other.nearest_segment
            self.nickname = 'Component_%d' % NovelSegment.nicknames[None]
            NovelSegment.nicknames[None] += 1
            return

        self.datasets = {other.name: other}
        self.seq = other.seq
        self.neighbors = []
        
        nearest = min_distance(self.seq, NovelSegment.source)
        self.nearest_dist, self.nearest_segment = nearest
        self.nickname = '%s|%d%s' % (self.nearest_segment, self.nearest_dist,
                                     letter(NovelSegment.nicknames[nearest]))
        NovelSegment.nicknames[nearest] += 1

    def covered_j(self):
        return sum(map(attrgetter('covered_j'), self.datasets.values()))

    def coverage(self):
        return sum(map(attrgetter('coverage'), self.datasets.values()))

    def add(self, other):
        assert isinstance(other, SegmentInDataset)

        if other.name in self.datasets:
            self.datasets[other.name] = self.datasets[other.name].merge(other)
        else:
            self.datasets[other.name] = other

    def merge(self, other):
        assert isinstance(other, NovelSegment)
        for dataset in other.datasets.values():
            self.add(dataset)

    def distance(self, other):
        return distance(self.seq, other.seq)

    def add_edge_if_needed(self, other, dist):
        assert isinstance(other, NovelSegment)

        import math
        
        if dist > NovelSegment.hamming_threshold \
                and dist > math.ceil((1 - NovelSegment.shared_polymorphism_rate) \
                               * max(self.nearest_dist, other.nearest_dist)):
            return False
        self.neighbors.append((other, dist))
        other.neighbors.append((self, dist))
        return True

    def combination(self):
        l = NovelSegment.left_margin
        r = NovelSegment.right_margin
        segm_seq = NovelSegment.segments[self.nearest_segment].seq
        res = []

        for i in range(l, min(r, len(segm_seq))):
            if i - l >= len(self.seq):
                res.append('%d:-' % (i + 1))
            elif self.seq[i - l] != segm_seq[i]:
                res.append('%d:%s' % (i + 1, self.seq[i - l]))
        return ','.join(res)

    def decoration(self):
        return self.covered_j()

    def __str__(self):
        res = []
        res.append(self.nickname)
        res.append(str(self.decoration()))
        res.append(str(self.coverage()))
        res.append(','.join(map(str, sorted(self.datasets.values(),
                                            key=SegmentInDataset.decoration, reverse=True))))

        if NovelSegment.show_neighbors:
            if self.neighbors:
                res.append(','.join('%s:%d' % (neigh.nickname, dist)
                                    for neigh, dist in sorted(self.neighbors, key=itemgetter(1))))
            else:
                res.append('*')
        elif self.neighbors is not None:
            res.append('N:%d' % len(self.neighbors))

        combination = self.combination()
        res.append('%s(%d)%s' % (self.nearest_segment, combination.count(':'), combination))
        return '\t'.join(res)

    def full_info(self, include_original):
        res = [str(self)]
        res.extend('  %s' % d.full_info(include_original) for d in 
                    sorted(self.datasets.values(), key=SegmentInDataset.decoration, reverse=True))
        return '\n'.join(res) + '\n'
        

def load_novel_segments(f, dataset_name, novel, min_dist):
    print(dataset_name, end=' ')
    sys.stdout.flush()

    line = next(f)
    while line.startswith('#'):
        line = next(f)

    for line in f:
        try:
            novel_segment = SegmentInDataset(dataset_name, line)
        except UnknownSegment as err:
            # print('Error')
            # TODO Warning
            continue
        distances = []

        for other in novel:
            distances.append(other.distance(novel_segment))
            if distances[-1] == 0:
                other.add(novel_segment)
                break
        else:
            novel_segment = NovelSegment(novel_segment)
            if novel_segment.nearest_dist < min_dist:
                continue
            for i, other in enumerate(novel):
                other.add_edge_if_needed(novel_segment, distances[i])
            novel.append(novel_segment)


def load_datasets(args):
    novel = []
    for line in args.datasets:
        if line.startswith('#'):
            continue

        name, path = line.strip().split()
        filename = os.path.join(os.path.dirname(args.datasets.name), path)
        with open(filename) as f:
            load_novel_segments(f, name, novel, args.source_dist)
    return novel    


def traversal(start, visited):
    component = []
    d = deque()
    d.append(start)
    while d:
        current = d.popleft()
        if current in visited:
            continue
        visited.add(current)
        component.append(current)

        for neighbor, _ in current.neighbors:
            if neighbor not in visited:
                d.append(neighbor)
    return component


def components(novel):
    res = []
    visited = set()
    for segment in novel:
        if segment in visited:
            continue
        res.append(traversal(segment, visited))
    return res


def write_group(summed_group, group, args):
    args.output.write('\n')
    args.output.write(str(summed_group))
    args.output.write('\n')
    for segm in group:
        args.output.write(str(segm))
        args.output.write('\n')

    args.full.write('\n')
    args.full.write(str(summed_group))
    args.full.write('\n')
    for segm in group:
        args.full.write(segm.full_info(True))
    args.full.write('\n')
    
    args.components.write('\n')
    args.components.write(summed_group.full_info(False))


def write_groups(groups, args):
    args.output.write('# %s\n# Short info\n' % ' '.join(sys.argv))
    args.full.write('# %s\n# Full info\n' % ' '.join(sys.argv))
    args.components.write('# %s\n# Components info\n' % ' '.join(sys.argv))
    
    new_groups = []
    for group in groups:
        group.sort(key=NovelSegment.decoration, reverse=True)
        summed_group = NovelSegment(group[0])
        for segm in group[1:]:
            summed_group.merge(segm)
        new_groups.append((summed_group, group))

    for summed_group, group in sorted(new_groups, key=lambda t: t[0].decoration(),
                                      reverse=True):
        write_group(summed_group, group, args)


def run(args):
    import itertools
    segments = combinations_to_segments.read_segments(itertools.chain(*args.v_segments))
    l, r = args.range
    l -= 1

    source = [(segment.seq[l:r], segment.name) for segment in segments.values()]

    NovelSegment.left_margin = l
    NovelSegment.right_margin = r
    NovelSegment.j_coverage = args.j_coverage
    NovelSegment.segments = segments
    NovelSegment.source = source
    NovelSegment.hamming_threshold = args.hamming
    NovelSegment.shared_polymorphism_rate = args.shared
    NovelSegment.show_neighbors = args.show_neighbors

    novel = load_datasets(args)
    groups = components(novel)

    write_groups(groups, args)


def main():
    parser = argparse.ArgumentParser(description='Group combinations by alleles and decorate them',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -c File -v File -o File [-l File] [args]')
    in_args = parser.add_argument_group('Input arguments')
    in_args.add_argument('-d', '--datasets', help='Input: file with lines <name> <combinations.csv>',
                         type=argparse.FileType(), required=True, metavar='File')
    in_args.add_argument('-v', '--v-segments', help='Input: source V segments', nargs='+',
                         type=argparse.FileType(), required=True, metavar='File')

    out_args = parser.add_argument_group('Output arguments')
    out_args.add_argument('-o', '--output', help='Short output',
                          type=argparse.FileType('w'), required=True, metavar='File')
    out_args.add_argument('-f', '--full', help='Full output (optional)', default='/dev/null',
                          type=argparse.FileType('w'), metavar='File')
    out_args.add_argument('-c', '--components', help='Output with components description (optional)',
                          type=argparse.FileType('w'), metavar='File', default='/dev/null')
    out_args.add_argument('--show-neighbors', help='Show neighbors in output files',
                          action='store_true')

    filter_args = parser.add_argument_group('Grouping and filtering arguments')
    filter_args.add_argument('--range', type=int, nargs=2, default=(45, 290), metavar='Int',
                             help='Positions range (default: 45 290)')
    filter_args.add_argument('--source-dist', type=int, default=4, metavar='Int',
                             help='Minimal distance to source V segments (default: 4)')
    filter_args.add_argument('--j-coverage', type=int, default=50, metavar='Int',
                             help='With this coverage J segment will be\n'
                                  'counted as significant (default: 50)')
    filter_args.add_argument('--hamming', type=int, default=2, metavar='Int',
                             help='If distance between two novel segments is at\n'
                                  'most <--hamming> - these segments will be\n'
                                  'automatically placed in one group (default: 2)')
    filter_args.add_argument('--shared', type=float, default=0.9, metavar='Float',
                             help='Two novel segments will be placed in one group\n'
                                  'if they share <--shared> per cent of\n'
                                  'mismatches (default: 0.9)')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other.add_argument('--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
