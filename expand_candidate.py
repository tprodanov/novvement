#!/usr/bin/env python3


def load_candidate(f):
    import collections
    candidate = collections.defaultdict(set)

    while next(f).startswith('#'):
        pass
    for line in f:
        segment, position, alt = line.strip().split()
        position = int(position)
        candidate[segment].add((position, alt))

    return candidate


class Combination:
    def __init__(self, combination):
        import collections
        self.combination = combination
        self.coverage = 0
        self.expansion = collections.Counter()

    def inc_coverage(self):
        self.coverage += 1

    def add_error(self, pos_alt):
        self.expansion[pos_alt] += 1

    def candidate(self, ratio):
        for pos_alt, count in self.expansion.most_common():
            if count < self.coverage * ratio:
                return
            yield pos_alt

    def write_long(self, segment, ratio, outp):
        combination_str = ','.join('%s:%s' % item for item in sorted(self.combination))
        if not combination_str:
            combination_str = '*'

        for pos, alt in sorted(self.combination):
            current = '{segment}\t{pos}\t{alt}\t{candidate}\t{combination}\n' \
                            .format(segment=segment, pos=pos, alt=alt,
                                    candidate='old', combination=combination_str)
            for i in range(self.coverage):
                outp.write(current)

        for (pos, alt), count in self.expansion.most_common():
            current = '{segment}\t{pos}\t{alt}\t{candidate}\t{combination}\n' \
                            .format(segment=segment, pos=pos, alt=alt,
                                    candidate='new' if count >= self.coverage * ratio else 'none',
                                    combination=combination_str)
            for i in range(count):
                outp.write(current)


def analyze_read(read_errors, segment_candidate, segment_combinations):
    current_candidate = []
    for pos_alt in read_errors:
        if pos_alt in segment_candidate:
            current_candidate.append(pos_alt)

    current_candidate = tuple(sorted(current_candidate))
    if current_candidate not in segment_combinations:
        segment_combinations[current_candidate] = Combination(current_candidate)
    combination = segment_combinations[current_candidate]
    
    combination.inc_coverage()
    for pos_alt in read_errors:
        if pos_alt not in segment_candidate:
            combination.add_error(pos_alt)


def analyze_all_reads(f, candidate):
    import collections
    combinations = collections.defaultdict(dict)
    
    last_read = None
    last_segment = None
    read_errors = None
    while next(f).startswith('#'):
        pass
    for line in f:
        read, segment, position, alt = line.strip().split('\t')
        position = int(position)
        if read != last_read:
            if read_errors:
                analyze_read(read_errors, candidate[last_segment], combinations[last_segment])
            last_read = read
            last_segment = segment
            read_errors = []
        read_errors.append((position, alt))

    if read_errors:
        analyze_read(read_errors, candidate[last_segment], combinations[last_segment])

    return combinations


def write_candidate(candidate, combinations, coverage, ratio, outp):
    import sys

    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('segment\tposition\talt\n')
    for segment, segment_candidate in candidate.items():
        for combination in combinations[segment].values():
            if combination.coverage >= coverage:
                segment_candidate |= set(combination.candidate(ratio))

        for pos, alt in sorted(segment_candidate):
            outp.write('%s\t%d\t%s\n' % (segment, pos, alt))


def write_long(combinations, coverage, ratio, outp):
    import sys

    outp.write('# %s\n' % ' '.join(sys.argv))
    outp.write('segment\tpos\talt\tcandidate\tcombination\n')

    for segment, segment_combinations in combinations.items():
        for combination in segment_combinations.values():
            if combination.coverage >= coverage:
                combination.write_long(segment, ratio, outp)


def main():
    import argparse
    from extra._version import __version__

    parser = argparse.ArgumentParser(description='Expand candidate polymorphisms',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -f File -c File -o File [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-f', '--filtered', help='Input csv file with filtered errors',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-c', '--candidate', help='Input csv file with candidate polymorphisms',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-o', '--output', help='Output csv file',
                         type=argparse.FileType('w'), required=True, metavar='File')
    io_args.add_argument('--long', action='store_true',
                         help='Write output in a long format')

    exp_args = parser.add_argument_group('Expansion arguments')
    exp_args.add_argument('-t', '--threshold', help='Combination coverage threshold (default: %(default)s)',
                          type=int, metavar='Int', default=50)
    exp_args.add_argument('-r', '--ratio', type=float, metavar='Float', default=0.75,
                          help='Candidate polymorphisms coverage ratio to combination coverage\n'
                          '(default: %(default)s)')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
    candidate = load_candidate(args.candidate)
    combinations = analyze_all_reads(args.filtered, candidate)

    if args.long:
        write_long(combinations, args.threshold, args.ratio, args.output)
    else:
        write_candidate(candidate, combinations, args.threshold, args.ratio, args.output)


if __name__ == '__main__':
    main()

