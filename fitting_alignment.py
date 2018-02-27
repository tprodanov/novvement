#!/usr/bin/env python3


def main():
    import argparse
    import os
    from extra._version import __version__
    import subprocess
    import sys

    parser = argparse.ArgumentParser(description='Fitting alignment',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -i File -o File [-t Int] [-d]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', type=argparse.FileType(), metavar='File',
                         help='Input fasta file with cropped sequences')
    io_args.add_argument('-s', '--segments', metavar='File',
                         help='Input fasta file with germline segments')
    io_args.add_argument('-o', '--output', type=argparse.FileType('w'), metavar='File',
                         help='Output csv file')

    align_args = parser.add_argument_group('Alignment arguments')
    align_args.add_argument('-M', '--mismatch', type=float, metavar='Float', default=-1,
                            help='Mismatch penalty (default: %(default)s)')
    align_args.add_argument('-O', '--open-gap', type=float, metavar='Float', default=-4,
                            help='Penalty for opening a gap (default: %(default)s)')
    align_args.add_argument('-E', '--extend-gap', type=float, metavar='Float', default=-1,
                            help='Penalty for extending a gap (default: %(default)s)')
    
    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)
    other.add_argument('--verbose', action='store_true',
                       help='Show command line')

    args = parser.parse_args()

    args.output.write('# %s\n' % ' '.join(sys.argv))
    args.output.write('subset\tsegment\tscore\tmismatches\tindels\tconsensus\n')
    args.output.flush()
    current_dir = os.path.dirname(__file__)
    
    if args.verbose:
        print('cat %s | %s/bin/fitting_alignment %s %.1f %.1f %.1f > %s' % (args.input.name, current_dir,
                    args.segments, args.mismatch, args.open_gap, args.extend_gap, args.output.name))

    returncode = subprocess.run(['%s/bin/fitting_alignment' % current_dir,
                                 args.segments,
                                 str(args.mismatch),
                                 str(args.open_gap),
                                 str(args.extend_gap)],
                                stdin=args.input,
                                stdout=args.output).returncode
    if returncode != 0:
        sys.stderr.write('Error in constructing fitting alignment\n')
        exit(1)


if __name__ == '__main__':
    main()
