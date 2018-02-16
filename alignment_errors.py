#!/usr/bin/env python3


from enum import Enum


class ErrorType(Enum):
    match = 0,
    mismatch = 1,
    insert = 2,
    deletion = 3


class Error:
    def __init__(self, error_type, pos, ref, alt):
        self.error_type = error_type
        self.pos = pos
        self.ref = ref
        self.alt = alt

    def long_str(self):
        assert self.error_type != ErrorType.match
        if self.error_type == ErrorType.mismatch:
            return 'mismatch\t%d\t%s\t%s' % (self.pos, self.ref, self.alt)
        if self.error_type == ErrorType.insert:
            return 'insert\t%d\t%s\t%s' % (self.pos, '-' * self.ref, self.alt)
        if self.error_type == ErrorType.deletion:
            return 'deletion\t%d\t%s\t%s' % (self.pos, self.ref, '-' * self.alt)

    def short_str(self):
        assert self.error_type != ErrorType.match
        alt = self.alt
        if self.error_type == ErrorType.insert:
            alt = 'i%s' % alt
        elif self.error_type == ErrorType.deletion:
            alt = 'd%d' % alt
        return '%d:%s' % (self.pos, alt)



def detect_error_type(nt1, nt2):
    if nt1 == '-':
        return ErrorType.deletion
    if nt2 == '-':
        return ErrorType.insert
    return ErrorType.mismatch if nt1 != nt2 else ErrorType.match


def alignment_errors(read_seq, segment_seq):
    pos = 1
    current_insert = ''
    current_deletion = 0

    for i, (nt1, nt2) in enumerate(zip(read_seq, segment_seq)):
        error_type = detect_error_type(nt1, nt2)
        if error_type == ErrorType.mismatch or error_type == ErrorType.match:
            if current_insert:
                yield Error(ErrorType.insert, pos, len(current_insert), current_insert)
                current_insert = ''
            elif current_deletion:
                yield Error(ErrorType.deletion, pos - current_deletion,
                            segment_seq[i - current_deletion : i], current_deletion)
                current_deletion = 0
            elif error_type == ErrorType.mismatch:
                yield Error(ErrorType.mismatch, pos, nt2, nt1)

        elif error_type == ErrorType.insert:
            current_insert += nt1
        elif error_type == ErrorType.deletion:
            current_deletion += 1

        if error_type != ErrorType.insert:
            pos += 1

    if current_insert:
        yield Error(ErrorType.insert, pos, len(current_insert), current_insert)
    elif current_deletion:
        yield Error(ErrorType.deletion, pos - current_deletion,
                    segment_seq[i - current_deletion : i], current_deletion)


def short_output(read_name, segment_name, errors, outp):
    if errors:
        outp.write('%s\t%s\t%s\n' % (read_name, segment_name, ','.join(map(Error.short_str, errors))))
    else:
        outp.write('%s\t%s\t*\n' % (read_name, segment_name))


def long_output(read_name, segment_name, errors, outp):
    if not errors:
        outp.write('%s\t%s\t0\t*\t*\t*\n' % (read_name, segment_name))
    for error in errors:
        outp.write('%s\t%s\t%s\n' % (read_name, segment_name, error.long_str()))


def analyze_next(inp, outp, outp_fn):
    import re

    read = next(inp)
    segment = next(inp)

    read_name = re.search(r'READ:([^|\n]*)', read.name).group(1)
    segment_name = re.search(r'GENE:([^|\n]*)', segment.name).group(1)
    outp_fn(read_name, segment_name, list(alignment_errors(read.seq, segment.seq)), outp)


def analyze_all(inp, outp, long_fmt=False):
    from extra import nt_string
    import sys

    outp.write('# %s\n' % ' '.join(sys.argv))
    if long_fmt:
        outp.write('read\tsegment\ttype\tposition\tref\talt\n')
        outp_fn = long_output
    else:
        outp.write('read\tsegment\terrors\n')
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
