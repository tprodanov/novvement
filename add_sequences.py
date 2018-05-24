#!/usr/bin/env python3


from extra import nt_string


def change_sequence(segment_seq, consensus):
    if consensus == '*':
        return segment_seq
    consensus = [s.split(':') for s in consensus.split(',')]
    consensus = sorted((int(pos), alt) for pos, alt in consensus)

    new_seq = list(segment_seq)
    for pos, alt in consensus:
        if alt.startswith('d'):
            deletion = int(alt[1:])
            for i in range(deletion):
                new_seq[pos + i - 1] = ''
        elif alt.startswith('i'):
            insert = alt[1:]
            new_seq[pos - 1] = new_seq[pos - 1] + insert
        else:
            new_seq[pos - 1] = alt

    return ''.join(new_seq)


def process_line(segments, line, keep_empty, min_coverage, min_labels):
    subset, coverage, coverage_rate, labels, consensus = line.strip().split('\t')
    segment = subset.split(':')[0]
    coverage = int(coverage)
    coverage_rate = float(coverage_rate)
    labels = int(labels)
    if segment not in segments:
        return
    if not keep_empty and consensus == '*':
        return
    if coverage < min_coverage:
        return
    if labels < min_labels:
        return

    return subset, change_sequence(segments[segment], consensus)


def process_file(inp, segments, keep_empty, min_coverage, min_labels, pos_range, full_outp, cropped_outp):
    l, r = pos_range
    l -= 1
    while next(inp).startswith('#'):
        pass

    for line in inp:
        res = process_line(segments, line, keep_empty, min_coverage, min_labels)
        if res is None:
            continue
        subset, seq = res
        nt_string.write_fasta_entry(subset, seq, full_outp)
        nt_string.write_fasta_entry(subset, seq[l:r], cropped_outp)


def main():
    import argparse
    from extra._version import __version__

    parser = argparse.ArgumentParser(description='Convert read subset consensus to sequences',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -i File -s File -o File [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', type=argparse.FileType(), metavar='File',
                         help='Input csv file with summary over read subsets')
    io_args.add_argument('-s', '--segments', type=argparse.FileType(), metavar='File',
                         help='Input fasta file with germline segments')
    io_args.add_argument('-o', '--output', type=argparse.FileType('w'), metavar=('File', 'File'), nargs=2,
                         help='Two output fasta files, containing full and cropped sequences')

    opt_args = parser.add_argument_group('Optional arguments')
    opt_args.add_argument('-r', '--range', type=int, metavar=('Int', 'Int'), nargs=2, default=[40, 290],
                          help='Positions range for the cropped sequences (default: %(default)s')
    opt_args.add_argument('--keep-empty', action='store_true', help='Keep empty consensus sequences')
    opt_args.add_argument('-c', '--coverage', type=int, metavar='Int', default=0,
                          help='Minimum coverage (default: %(default)s)')
    opt_args.add_argument('-l', '--labels', type=int, metavar='Int', default=0,
                          help='Minimum number of labels (default: %(default)s)')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
 
    segments = {segment.name: segment.seq for segment in nt_string.read_fasta(args.segments)}
    process_file(args.input, segments, keep_empty=args.keep_empty, min_coverage=args.coverage,
                 min_labels=args.labels, pos_range=args.range, full_outp=args.output[0], cropped_outp=args.output[1])


if __name__ == '__main__':
    main()
