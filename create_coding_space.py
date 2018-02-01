#!/usr/bin/env python3


class Vertex:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.edges = set()

    def add_edge(self, other):
        self.edges.add(other)
        other.edges.add(self)

    def remove_edge(self, other):
        self.edges.remove(other)
        other.edges.remove(self)

    def remove_neighbours(self):
        neighbours = list(self.edges)
        for other in neighbours:
            other.edges.remove(self)
        self.edges.clear()
        return neighbours

    def distance(self, other):
        return sum(nt1 != nt2 for nt1, nt2 in zip(self.seq, other.seq))

    def group(self):
        return self.name.split('*')[0]

    def allele(self):
        return int(self.name.split('*')[1])

    def len_edges(self):
        return len(self.edges)

    def __str__(self):
        return '{name} : {edges} edges'.format(name=self.name, edges=self.len_edges())


class Graph:
    removed_total = 0
    def __init__(self, vertices=None):
        self.vertices = set()
        if vertices:
            for vertex in vertices:
                self.vertices.add(vertex)

    def add_seq(self, name, seq, tau):
        current = Vertex(name, seq)
        for vertex in self.vertices:
            if vertex.distance(current) <= tau:
#                if current.allele() == vertex.allele() == 1:
#                    print(current.name, vertex.name, vertex.distance(current), sep='\t')
                vertex.add_edge(current)

        self.vertices.add(current)

    @staticmethod
    def traverse(start, visited):
        import collections
        to_visit = collections.deque()

        to_visit.append(start)
        visited.add(start)

        while to_visit:
            current = to_visit.popleft()
            for other in current.edges:
                if other not in visited:
                    to_visit.append(other)
                    visited.add(other)
            yield current

    def __len__(self):
        return len(self.vertices)

    def components(self):
        visited = set()
        for vertex in self.vertices:
            if vertex not in visited:
                yield list(self.traverse(vertex, visited))

    def remove_by(self, selector):
        for vertex in selector(self.vertices):
            # print('Removing', vertex)
            vertex.remove_neighbours()
            self.vertices.remove(vertex)
            Graph.removed_total += 1

    def remove_first_allele_neighbours(self):
        to_delete = set()

        for vertex in self.vertices:
            if vertex in to_delete:
                continue
            if vertex.allele() == 1:
                neighbours = vertex.remove_neighbours()
                to_delete |= set(neighbours)
        
        print('Deleting neighbours of the dominant alleles:', len(to_delete))
        for vertex in to_delete:
            vertex.remove_neighbours()
            self.vertices.remove(vertex)


def decompose_graph(graph, selector):
    import collections
    queue = collections.deque()

    for component in graph.components():
        if len(component) == 1:
            yield component[0]
        else:
            queue.append(Graph(component))

    while queue:
        current = queue.popleft()

        current.remove_by(selector)
        for component in current.components():
            if len(component) == 1:
                yield component[0]
            else:
                queue.append(Graph(component))


def argmax(vertices, fun):
    m = -1000000
    selected = None
    for vertex in vertices:
        if fun(vertex) > m:
            m = fun(vertex)
            selected = vertex
    return selected


def select_by_edges(vertices):
    return [argmax(vertices, Vertex.len_edges)]

def select_opposite_edges(vertices):
    return list(argmax(vertices, Vertex.len_edges).edges)


def select_by_coverage(coverages):
    fitness = lambda vertex: -coverages.get(vertex.name, 0)
    def inner(vertices):
        return [argmax(vertices, fitness)]
    return inner


def select_by_alleles(max_alleles):
    fitness = lambda vertex: (vertex.allele() - 1) / max_alleles[vertex.group()]
    def inner(vertices):
        return [argmax(vertices, fitness)]
    return inner


def prepare_selector(args, all_vertices):
    import collections
    if args.resolution == 'edges':
        return select_by_edges
    if args.resolution == 'nedges':
        return select_opposite_edges

    if args.resolution == 'coverage':
        coverages = {}
        for line in args.coverage:
            line = line.strip().split()
            coverages[line[0]] = int(line[1])
        return select_by_coverage(coverages)

    if args.resolution == 'alleles':
        alleles = collections.defaultdict(int)
        for vertex in all_vertices:
            alleles[vertex.group()] = max(alleles[vertex.group()], vertex.allele())
        return select_by_alleles(alleles)


def read_input(f, tau):
    graph = Graph()
    seqs = {}
    for line in f:
        line = line.strip()
        if not line:
            continue
        name, seq = line.split()
        positions = [i for i, nt in enumerate(seq) if nt != '-' if nt.upper() == nt]
        assert positions

        graph.add_seq(name, seq[min(positions) : max(positions) + 1], tau)
        seqs[name] = seq
    return graph, seqs


def write_output(vertices, seqs, outp_fa, outp_aln):
    from operator import attrgetter
    from extra import nt_string
    vertices.sort(key=attrgetter('name'))

    if outp_aln:
        for i, vertex in enumerate(vertices):
            if i and vertex.group() != vertices[i - 1].group():
                outp_aln.write('\n')
            outp_aln.write('%-20s %s\n' % (vertex.name, seqs[vertex.name]))

    for vertex in vertices:
        seq = seqs[vertex.name].upper().replace('-', '')
        nt_string.write_fasta_entry(vertex.name, seq, outp_fa)


def describe_graph(graph):
    from collections import Counter

    print('Total %d vertices' % len(graph.vertices))
    groups = { vertex.group() for vertex in graph.vertices }
    print('Total %d groups' % len(groups))

    components = list(graph.components())
    dominant_counter = Counter()
    for component in components:
        dominant = sum(vertex.allele() == 1 for vertex in component)
        dominant_counter[dominant] += 1

    print('%d orphans' % dominant_counter[0])
    print('%d single-centre components' % dominant_counter[1])
    print('%d multi-centre components' %
            (sum(dominant_counter.values()) - dominant_counter[1] - dominant_counter[0]))
    print(dominant_counter)



def main():
    import argparse
    from extra._version import __version__
    parser = argparse.ArgumentParser(description='Create coding space using Hamming graph',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -i File [-a File] -o File [-c File] [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', help='File containing multiple alignment in format <name> <seq>',
                         required=True, metavar='File', type=argparse.FileType())
    io_args.add_argument('-o', '--output', help='Output fasta file',
                         required=True, metavar='File', type=argparse.FileType('w'))
    io_args.add_argument('-a', '--alignment', help='Optional: output multiple alignment',
                         metavar='File', type=argparse.FileType('w'))
    io_args.add_argument('-c', '--coverage', help='Optional input file containing segment coverage',
                         metavar='File', type=argparse.FileType())
    io_args.add_argument('-v', '--verbose', help='Describe graph', action='store_true')

    gr_args = parser.add_argument_group('Graph arguments')
    gr_args.add_argument('-t', '--tau', help='Hamming distance between sequences (default: %(default)s)',
                         type=int, metavar='Int', default=4)
    gr_args.add_argument('-1', '--first-allele', help='Remove neighbours of the first allele (*01)',
                         action='store_true')
    gr_args.add_argument('-r', '--resolution',
                         help='Component resolution method:\n'
                              '  edges    - remove vertices with the most edges\n'
                              '  nedges   - remove vertices neighboring vertices with the most edges\n'
                              '  coverage - remove low covered vertices (required -c File)\n'
                              '  alleles  - remove vertices with the biggest allele numbers\n'
                              '(default: %(default)s)',
                         choices=['edges', 'nedges', 'coverage', 'alleles'], default='edges')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)
    
    args = parser.parse_args()
    
    graph, seqs = read_input(args.input, args.tau)
    if args.verbose:
        describe_graph(graph)
    if args.first_allele:
        graph.remove_first_allele_neighbours()

    selector = prepare_selector(args, graph.vertices)
    vertices = list(decompose_graph(graph, selector))
    write_output(vertices, seqs, args.output, args.alignment)
    print('Removed by selector:', Graph.removed_total)
    


if __name__ == '__main__':
    main()

