

class FastaRead:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq.upper()

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
        self.seq = seq.upper()
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


def read_fastaq(f, default_ext='fq'):
    if f.name.endswith('fq') or f.name.endswith('fastq'):
        filetype = 'fq'
    elif f.name.endswith('fa') or f.name.endswith('fasta'):
        filetype = 'fa'
    else:
        filetype = default_ext

    if filetype == 'fq':
        return read_fastq(f)
    else:
        return read_fasta(f)

