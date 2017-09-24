#!/usr/bin/env python3


def detect_error_type(nt1, nt2):
    if nt1 == '-':
        return 'deletion'
    if nt2 == '-':
        return 'insert'
    return 'mismatch' if nt1 != nt2 else 'match'


def alignment_errors(read_seq, segment_seq):
    pos = 1
    current_insert = ''
    current_deletion = 0

    for i, (nt1, nt2) in enumerate(zip(read_seq, segment_seq)):
        error_type = detect_error_type(nt1, nt2)
        if error_type == 'mismatch' or error_type == 'match':
            if current_insert:
                yield 'insert', pos, '-' * len(current_insert), current_insert
                current_insert = ''
            elif current_deletion:
                yield 'deletion', pos - current_deletion, \
                        segment_seq[i - current_deletion : i], '-' * current_deletion
                current_deletion = 0

            if error_type == 'mismatch':
                yield error_type, pos, nt2, nt1

        elif error_type == 'insert':
            current_insert += nt1
        elif error_type == 'deletion':
            current_deletion += 1

        if error_type != 'insert':
            pos += 1

    if current_insert:
        yield 'insert', pos, '-' * len(current_insert), current_insert
    elif current_deletion:
        yield 'deletion', pos - current_deletion, \
            segment_seq[i - current_deletion : i], '-' * current_deletion


def short_output(read_name, segment_name, error_type, pos, ref, alt, outp):
    if error_type == 'insert':
        alt = 'i' + alt
    elif error_type == 'deletion':
        alt = 'd' + str(len(alt))
    outp.write('%s\t%s\t%d\t%s\n' % (read_name, segment_name, pos, alt))


def long_output(read_name, segment_name, error_type, pos, ref, alt, outp):
    outp.write('%s\t%s\t%s\t%d\t%s\t%s\n' % (read_name, segment_name, error_type,
               pos, ref, alt))


def analyze_next(inp, outp, outp_fn):
    import re

    read = next(inp)
    segment = next(inp)

    read_name = re.search(r'READ:([^|\n]*)', read.name).group(1)
    segment_name = re.search(r'GENE:([^|\n]*)', segment.name).group(1)

    for error_type, pos, ref, alt in alignment_errors(read.seq, segment.seq):
        outp_fn(read_name=read_name, segment_name=segment_name, error_type=error_type,
                pos=pos, ref=ref, alt=alt, outp=outp)


def analyze_all(inp, outp, long_fmt=False):
    from extra import nt_string
    import sys

    outp.write('# %s\n' % ' '.join(sys.argv))
    if long_fmt:
        outp.write('read\tsegment\ttype\tposition\tref\talt\n')
        outp_fn = long_output
    else:
        outp.write('read\tsegment\tposition\talt\n')
        outp_fn = short_output

    inp = nt_string.read_fastaq(inp, default_ext='fa')

    try:
        while True:
            analyze_next(inp, outp, outp_fn)
    except StopIteration:
        pass


def main():
    import argparse
    from extra._version import __version__

    parser = argparse.ArgumentParser(description='Convert alignment fasta file to csv file',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -i File -o File [--long]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', help='Input alignment file (fasta or fastq)',
                         type=argparse.FileType(), required=True, metavar='File')
    io_args.add_argument('-o', '--output', help='Output csv file',
                         type=argparse.FileType('w'), required=True, metavar='File')
    io_args.add_argument('--long', help='Long output format', action='store_true')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
    analyze_all(args.input, args.output, args.long)
 

if __name__ == '__main__':
    main()
