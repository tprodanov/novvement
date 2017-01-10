

class FastaRead:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def __str__(self):
        return '%s\n%s\n' % (self.name, self.seq)


def read_fasta(f):
    name = next(f).strip()
    current_sequence = []
    for line in f:
        if line[0] == '>':
            yield FastaRead(name, ''.join(current_sequence))
            current_sequence = []
            name = line.strip()
        else:
            current_sequence.append(line.strip().upper())
    yield FastaRead(name, ''.join(current_sequence))


class FastqRead:
    def __init__(self, name, seq, quality):
        self.name = name
        self.seq = seq
        self.quality = quality

    @staticmethod
    def from_file(f):
        name = next(f).strip()
        seq = next(f).strip()
        next(f)
        quality = next(f).strip()
        return FastqRead(name, seq, quality)

    def mean_quality(self):
        quality_sum = 0.0;
        for q in self.quality:
            quality_sum += ord(q) - ord('!')
        return quality_sum / len(self.quality)

    def __str__(self):
        return '%s\n%s\n+\n%s\n' % (self.name, self.seq, self.quality)


def read_fastq(f):
    while True:
        yield FastqRead.from_file(f)


class UnexpectedLetter(Exception):
    pass


complement_nt = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
nt_number = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def reverse_complement(seq):
    return ''.join(complement_nt[nt] for nt in reversed(seq))


def count_kmer(seq, start, k):
    assert(start >= 0 and start + k <= len(seq))
    kmer = 0
    for i in range(start, start + k):
        if seq[i] == 'N':  # TODO: What to do when N
            return -1
        kmer = kmer * 4 + nt_number[seq[i]]
    return kmer


def count_kmer_extend_n(seq, start, k):
    n_positions = []
    kmer = 0
    for i in range(start, start + k):
        if seq[i] == 'N':  # TODO: What to do when N
            if len(n_positions) == 2:
                return []
            n_positions.append(i - start)
        else:
            kmer = kmer * 4 + nt_number[seq[i]]
    if not n_positions:
        return [kmer]
    elif len(n_positions) == 1:
        ix_multiplier = 4 ** n_positions[0]
        return [kmer + ix_multiplier * i for i in range(4)]
    else:
        ix_multiplier0 = 4 ** n_positions[0]
        ix_multiplier1 = 4 ** n_positions[1]
        return [kmer + ix_multiplier0 * i + ix_multiplier1 * j for i in range(4) for j in range(4)]
