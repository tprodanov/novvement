#!/usr/bin/env python3


def create_parser():
    import argparse
    from extra._version import __version__
    parser = argparse.ArgumentParser(description='Extract novel segments from naive Rep-seq data',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -i File -s File -o Dir [args]') # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', help='File containing datasets. Format:\n'
                         '<name> <alignment fasta> <cdr3 fasta>',
                         required=True, metavar='File')
    io_args.add_argument('-s', '--segments', help='File with germline segments in fasta format',
                         required=True, metavar='File')
    io_args.add_argument('-S', '--clipped-segments',
                         help='(Optional) Already clipped germline segments in fasta format',
                         metavar='File')
    io_args.add_argument('-o', '--output', help='Output directory',
                         required=True, metavar='Dir')

    rt_args = parser.add_argument_group('Runtime arguments')
    rt_args.add_argument('-t', '--threads', help='Number of threads (default: %(default)s)',
                         type=int, metavar='Int', default=4)

    filt_args = parser.add_argument_group('Filtering errors')
    filt_args.add_argument('--range', help='Positions range (default: %(default)s)',
                           type=int, metavar=('Int', 'Int'), nargs=2, default=[40, 290])
    filt_args.add_argument('--segm-cov', '--segment-coverage', dest='segment_coverage',
                           type=int, metavar='Int', default=100,
                           help='Minimal segment coverage (default: %(default)s)')
    filt_args.add_argument('--keep-n', action='store_true',
                           help='Do not filter out mismatches with N as a substitute')

    det_args = parser.add_argument_group('Candidate polymorphisms detection')
    det_args.add_argument('--det-method', help='Candidate detection method (default: %(default)s)',
                          choices=['gap', 'quantile'], metavar='gap|quantile', default='gap')
    det_args.add_argument('--gap-size', metavar='Float', type=float, default=0.5,
                          help='Relative gap size (if --method gap) (default: %(default)s)')
    det_args.add_argument('--multiplier', metavar='Float', type=float, default=4,
                          help='Quantile multiplier (if --det-method quantile):\n'
                          'We select positions with number of errors at least\n'
                          '<--multiplier> * <--quantile quantile of errors per position>\n'
                          '(default: %(default)s)')
    det_args.add_argument('--quantile', metavar='Float', type=float, default=95,
                          help='Quantile used for comparison (if --det-method quantile)\n'
                          '(default: %(default)s)')
    
    sc_args = parser.add_argument_group('Spliting reads and finding a consensus')
    sc_args.add_argument('--spl-cov', '--split-coverage', dest='split_coverage',
                         type=int, metavar='INT', default=100,
                         help='Minimal coverage of a reads subset (default: %(default)s)')
    sc_args.add_argument('--cons-ratio', '--consensus-ratio', dest='consensus_ratio',
                         type=float, metavar='FLOAT', default=0.51,
                         help='Ratio for a new polymorphisms to be chosen.\n'
                              'Should be higher than 0.5 (default: %(default)s)')
            
    hl_args = parser.add_argument_group('Merging labels using Hamming graph')
    hl_args.add_argument('--labels-tau', metavar='Int', type=int, default=3,
                         help='Hamming distance between labels (cdr3s) (default: %(default)s)')

    f2_args = parser.add_argument_group('Filtering novel segments')
    f2_args.add_argument('--min-dist', metavar='Int', type=int, default=2,
                         help='Minimal Hamming distance to clipped segments (default: %(default)s)')
    
    hs_args = parser.add_argument_group('Merging sequences using Hamming graph')
    hs_args.add_argument('--seqs-tau', metavar='Int', type=int, default=1,
                         help='Hamming distance between sequences (default: %(default)s)')
    hs_args.add_argument('--seqs-same-lengths', action='store_true',
                         help='Forbid edges between sequences with different sizes')

    sign_args = parser.add_argument_group('Significance arguments')
    sign_args.add_argument('--label-cov', '--label-coverage', dest='label_coverage',
                           type=int, metavar='Int', default=5,
                           help='Minimum label coverage to contribute in significance\n'
                           '(default: %(default)s)')

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
