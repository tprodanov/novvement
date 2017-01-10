import sys
from operator import itemgetter


class Tee:
    def __init__(self, f, g=sys.stdout):
        self.f = f
        self.g = g

    def write(self, *args, **kwargs):
        self.f.write(*args, **kwargs)
        self.g.write(*args, **kwargs)


def v_alignment_to_mutations(v_alignment, potential_mutations):
    while next(v_alignment).startswith('#'):
        pass
    prev_read = None

    for line in v_alignment:
        gene, read, mut_type, position, nt = line.strip().split()
        position = int(position)

        if gene not in potential_mutations:
            continue
        if prev_read != read:
            if prev_read and read_mutations:
                read_mutations.sort(key=itemgetter(0))
                yield prev_read, prev_gene, read_mutations
            prev_read = read
            prev_gene = gene
            read_mutations = []
            gene_mutations = potential_mutations[gene]

        if (position, nt) in gene_mutations:
            read_mutations.append((position, nt))
    if prev_read and read_mutations:
        read_mutations.sort(key=itemgetter(0))
        yield prev_read, prev_gene, read_mutations
