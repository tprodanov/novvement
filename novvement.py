#!/usr/bin/env python3


def create_parser():
    import argparse
    from extra._version import __version__
    parser = argparse.ArgumentParser(description='Extract novel segments from naive Rep-seq data',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -i File -s File -o Dir [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', metavar='File', required=True,
                         help='File containing datasets. Format:\n'
                         '<name> <alignment fasta> <cdr3 fasta>\n'
                         'datasets woth different properties (significantly different coverage,\n'
                         'different PCR sthrength, different level of memory cells) should be\n'
                         'separated by lines "! Group name". There must be at least one group.')
    io_args.add_argument('-s', '--segments', metavar='File', required=True,
                         help='File with germline segments in fasta format')
    io_args.add_argument('-o', '--output', metavar='Dir', required=True,
                         help='Output directory')

    rt_args = parser.add_argument_group('Runtime arguments')
    rt_args.add_argument('-t', '--threads', type=int, metavar='Int', default=4,
                         help='Number of threads (default: %(default)s)')

    filt_args = parser.add_argument_group('Filtering errors')
    filt_args.add_argument('--range', help='Positions range (default: %(default)s)',
                           type=int, metavar=('Int', 'Int'), nargs=2, default=[40, 290])
    filt_args.add_argument('--segm-cov', '--segment-coverage', dest='segment_coverage',
                           type=int, metavar='Int', default=100,
                           help='Minimal segment coverage (default: %(default)s)')
    filt_args.add_argument('--keep-n', action='store_true',
                           help='Do not filter out mismatches with N as a substitute')

    div_args = parser.add_argument_group('Dividing read sets')
    div_args.add_argument('--significance', type=float, metavar='Float', default=.001,
                          help='Division significance (default: %(default)s')
    div_args.add_argument('--pairs', type=int, metavar='Int', default=5,
                          help='Maximum number of evaluated pairs of errors (default: %(default)s)')
    div_args.add_argument('--subset-coverage', type=int, metavar='Int', default=50,
                          help='Minimum subset size (default: %(default)s)')
    div_args.add_argument('--single-peak', type=float, metavar=('Float', 'Float'), default=[.2, .3], nargs=2,
                          help='Two rates for single peak detection: noise rate and\n'
                               'true peak rate (default: %(default)s)')
            
    hl_args = parser.add_argument_group('Merging labels using Hamming graph')
    hl_args.add_argument('--labels-tau', metavar='Int', type=int, default=3,
                         help='Hamming distance between labels (cdr3s) (default: %(default)s)')

    fit_args = parser.add_argument_group('Fitting alignment arguments')
    fit_args.add_argument('--mismatch', type=float, metavar='Float', default=-1,
                          help='Mismatch penalty (default: %(default)s)')
    fit_args.add_argument('--open-gap', type=float, metavar='Float', default=-4,
                          help='Penalty for opening a gap (default: %(default)s)')
    fit_args.add_argument('--extend-gap', type=float, metavar='Float', default=-1,
                          help='Penalty for extending a gap (default: %(default)s)')

    cl_args = parser.add_argument_group('Classifying putative novel segments')
    cl_args.add_argument('--threshold', metavar='Float', type=float, default=0.05,
                         help='Maximal p-value (default: %(default)s)')

    dr_args = parser.add_argument_group('Drawing arguments')
    dr_args.add_argument('--skip-plots', action='store_true',
                         help='Do not draw plots')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)
    return parser


def main():
    import os
    import snakemake

    parser = create_parser()
    args = parser.parse_args()

    args.script_dir = os.path.abspath(os.path.dirname(__file__))

    snakefile = os.path.join(args.script_dir, 'Snakefile')
    snakemake.snakemake(snakefile, cores=args.threads, config=vars(args))


if __name__ == '__main__':
    main()
