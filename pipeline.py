#!/usr/bin/python3

import argparse
import os
import sys
import io
import subprocess


def mkdir(path):
    try:
        os.mkdir(path)
    except FileExistsError:
        pass


def log_header(log, args, human_args):
    log.write('# Command: %s\n' % ' '.join(sys.argv))
    args = vars(args)

    for k, v in human_args:
        if not v:
            log.write('#\t%s\n' % k)
        else:
            v = args[v]
            log.write('#\t\t%s: %s\n' % (k, v.name if isinstance(v, io.IOBase) else v))
    log.write('\n')


def make_datasets_file(dir, log, datasets):
    with open(os.path.join(dir, 'datasets.csv'), 'w') as f:
        log.write('\n# Generate <name> <v_alignment.csv> file (%s)\n' % f.name)
        for name, _ in datasets:
            f.write('%s %s\n' % (name, os.path.join(name, 'alignment', 'v_alignment.csv')))


def gene_coverage(args, log, datasets):
    dir = args.output
    sys.stdout.write('Gene coverage\n')
    log.write('\n# Gene coverage\n')

    for name, path in datasets:
        with open(os.path.join(path, 'alignment_info.csv')) as inp, open(os.path.join(dir, name, 'alignment', 'gene_coverage.csv'), 'w') as outp:
            log.write("\tcat %s | cut -f5 | sort | uniq -c | awk '{t = $1; $1 = $2; $2 = t; print}' > %s\n" % (inp.name, outp.name))
            p1 = subprocess.Popen(['cut', '-f5'], stdin=inp, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(['sort'], stdin=p1.stdout, stdout=subprocess.PIPE)
            p3 = subprocess.Popen(['uniq', '-c'], stdin=p2.stdout, stdout=subprocess.PIPE)
            p4 = subprocess.Popen(['awk', '{t = $1; $1 = $2; $2 = t; print}'], stdin=p3.stdout, stdout=outp)
            p1.stdout.close()
            p2.stdout.close()
            p3.stdout.close()
            p4.communicate()


def j_hits(args, log, datasets):
    dir = args.output
    sys.stdout.write('J hit\n')
    log.write('\n# J hit\n')

    for name, path in datasets:
        with open(os.path.join(path, 'alignment_info.csv')) as inp, open(os.path.join(dir, name, 'alignment', 'j_hit.csv'), 'w') as outp:
            log.write("\tcat %s | cut -f1,9 > %s\n" % (inp.name, outp.name))
            p1 = subprocess.Popen(['cut', '-f1,9'], stdin=inp, stdout=outp)
            p1.communicate()


def v_alignment(args, log, datasets):
    dir = args.output
    script_path = os.path.dirname(__file__)

    sys.stdout.write('V alignment -> csv\n')
    log.write('\n# V alignment -> csv\n')
    for name, path in datasets:
        sys.stdout.write('%s ' % name)
        sys.stdout.flush()
        command = [os.path.join(script_path, 'v_alignment_mutations.py'),
                   '-i', os.path.join(path, 'v_alignments.fa'),
                   '-o', os.path.join(dir, name, 'alignment', 'v_alignment.csv')]
        log.write('\t%s\n' % ' '.join(command))
        subprocess.run(command)

    sys.stdout.write('\n')


def potential_mutations(args, log, datasets):
    dir = args.output
    script_path = os.path.dirname(__file__)

    sys.stdout.write('Detecting potential mutations\n')
    log.write('\n# Detecting potential mutations\n')
    for name, path in datasets:
        mkdir(os.path.join(dir, name, 'combinations_a'))
        sys.stdout.write('%s ' % name)
        sys.stdout.flush()
        command = [os.path.join(script_path, 'potential_mutations.py'),
                   '-v', os.path.join(dir, name, 'alignment', 'v_alignment.csv'),
                   '-g', os.path.join(dir, name, 'alignment', 'gene_coverage.csv'),
                   '-o', os.path.join(dir, name, 'combinations_a', 'potential_mutations.csv'),
                   '--range', args.range[0], args.range[1],
                   '--coverage', args.gene_coverage,
                   '--rate', args.mismatch_rate]
        command = [str(x) for x in command]
        log.write('\t%s\n' % ' '.join(command))
        subprocess.run(command)

    sys.stdout.write('\n')


def combinations(args, log, datasets, comb_dir):
    dir = args.output
    script_path = os.path.dirname(__file__)
    
    sys.stdout.write('Combining potential mutations\n')
    log.write('\n# Combining potential mutations\n')

    for name, path in datasets:
        sys.stdout.write('%s ' % name)
        sys.stdout.flush()
        command = [os.path.join(script_path, 'combinations.py'),
                   '-v', os.path.join(dir, name, 'alignment', 'v_alignment.csv'),
                   '-j', os.path.join(dir, name, 'alignment', 'j_hit.csv'),
                   '-m', os.path.join(dir, name, comb_dir, 'potential_mutations.csv')]
        if comb_dir.endswith('b'):
            command += ['-d', os.path.join(dir, 'datasets.csv')]
        command += ['-o', os.path.join(dir, name, comb_dir, 'combinations.csv')]
        if args.human_readable:
            command += ['-h', os.path.join(dir, name, comb_dir, 'human_combinations.txt')]
        command += ['--length', args.length,
                    '--cov-single-j', args.cov_single_j, '--cov-mult-j', args.cov_mult_j]

        command = [str(x) for x in command]
        log.write('\t%s\n' % ' '.join(command))
        subprocess.run(command)

    sys.stdout.write('\n')


def expand_mismatches(args, log, datasets):
    dir = args.output
    script_path = os.path.dirname(__file__)

    sys.stdout.write('Expanding combinations\n')
    log.write('\n# Expanding combinations\n')
    for name, path in datasets:
        mkdir(os.path.join(dir, name, 'combinations_b'))
        sys.stdout.write('%s ' % name)
        sys.stdout.flush()
        command = [os.path.join(script_path, 'expand_potential_mutations.py'),
                   '-v', os.path.join(dir, name, 'alignment', 'v_alignment.csv'),
                   '-m', os.path.join(dir, name, 'combinations_a', 'potential_mutations.csv'),
                   '-c', os.path.join(dir, name, 'combinations_a', 'combinations.csv'),
                   '-o', os.path.join(dir, name, 'combinations_b', 'potential_mutations.csv'),
                   '-e', os.path.join(dir, name, 'combinations_b', 'expanded_mutations.csv'),
                   '--range', args.range[0], args.range[1],
                   '--coverage', args.comb_coverage,
                   '--rate', args.expansion_rate]
        command = [str(x) for x in command]
        log.write('\t%s\n' % ' '.join(command))
        subprocess.run(command)

    sys.stdout.write('\n')


def significance(args, log, datasets):
    dir = args.output
    
    with open(os.path.join(dir, 'combinations_path.csv'), 'w') as f:
        log.write('\n# Generate <name> <combinations.csv> file (%s)\n' % f.name)
        for name, _ in datasets:
            f.write('%s %s\n' % (name, os.path.join(name, 'combinations_b', 'combinations.csv')))
    
    script_path = os.path.dirname(__file__)

    sys.stdout.write('Evaluating combination significance\n')
    log.write('\n# Evaluating combination significance\n')
    
    command = [os.path.join(script_path, 'significance.py'),
               '-i', os.path.join(dir, 'combinations_path.csv'),
               '-o', os.path.join(dir, 'significance.csv'),
               '-m', args.min_significance]
    command = [str(x) for x in command]
    log.write('\t%s\n' % ' '.join(command))
    subprocess.run(command)

    log.write('\n# Significance -> Combinations\n')
    with open(os.path.join(dir, 'significance.csv')) as inp, open(os.path.join(dir, 'combinations.csv'), 'w') as outp:
        log.write('\tcat %s | cut -f2- > %s\n' % (inp.name, outp.name))
        p1 = subprocess.Popen(['cut', '-f2-'], stdin=inp, stdout=outp)
        p1.communicate()


def generate(args, log, datasets):
    dir = args.output
    script_path = os.path.dirname(__file__)

    sys.stdout.write('Generating novel segments\n')
    log.write('\n# Generating novel segments\n')
    
    command = [os.path.join(script_path, 'generate_possible_igv.py'),
               '-c', os.path.join(dir, 'combinations.csv'),
               '-i', args.v_segments,
               '-o', os.path.join(dir, 'segments.fa')]
    command = [str(x) for x in command]
    log.write('\t%s\n' % ' '.join(command))
    subprocess.run(command)

    sys.stdout.write('Segments generated ::::: %s\n' % os.path.join(dir, 'segments.fa'))


def run(args, human_args):
    dir = args.output
    mkdir(args.output)
    if not args.force and os.listdir(dir):
        sys.stderr.write('Output folder is not empty. Try -f/--force\n')
        exit(1)
    
    datasets = [tuple(line.strip().split(' ', 1)) for line in args.input if not line.startswith('#')]

    with open(os.path.join(dir, 'log.txt'), 'w') as log:
        log_header(log, args, human_args)

        make_datasets_file(dir, log, datasets)

        for name, path in datasets:
            mkdir(os.path.join(dir, name))
            mkdir(os.path.join(dir, name, 'alignment'))
        
        gene_coverage(args, log, datasets)
        j_hits(args, log, datasets)
        v_alignment(args, log, datasets)
        potential_mutations(args, log, datasets)
        combinations(args, log, datasets, 'combinations_a')
        expand_mismatches(args, log, datasets)
        combinations(args, log, datasets, 'combinations_b')
        significance(args, log, datasets)
        generate(args, log, datasets)


def main():
    parser = argparse.ArgumentParser(description='Pipeline for detecting novel segments', add_help=False)
    io_args = parser.add_argument_group('inpur/output arguments')
    io_args.add_argument('-i', '--input', required=True, metavar='FILE', type=argparse.FileType(),
                         help='file with lines <name> <path to cdr folder>.'
                              'Should contain: v_alignments.fa, alignment_info.csv')
    io_args.add_argument('-v', '--v-segments', required=True, metavar='FILE', help='file containing V segments',
                         dest='v_segments')
    io_args.add_argument('-o', '--output', help='output directory', required=True, metavar='DIR')
    io_args.add_argument('-f', '--force', help='override files in output directory', action='store_true')
    human_args = [('Input/output', None), ('Input', 'input'), ('V segments', 'v_segments'), ('Output', 'output'), ('Force', 'force')]

    mutation_args = parser.add_argument_group('mutation detection arguments')
    mutation_args.add_argument('--range', help='positions range (default: [60, 280])',
                               metavar=('INT', 'INT'), nargs=2, default=[60, 280], type=int)
    mutation_args.add_argument('--gene-coverage', help='gene coverage threshold (default: 200)',
                               type=int, default=200, metavar='INT', dest='gene_coverage')
    mutation_args.add_argument('--mismatch-rate', help='mismatch rate threshold (default: 0.1)',
                               type=float, default=0.1, metavar='FLOAT', dest='mismatch_rate')
    human_args += [('Mutation detection', None), ('Range', 'range'), ('Gene coverage', 'gene_coverage'), ('Mismatch rate', 'mismatch_rate')]

    comb_args = parser.add_argument_group('combination detection arguments')
    comb_args.add_argument('--length', help='min combination length (default: 2)',
                           type=int, default=2, metavar='INT')
    comb_args.add_argument('--cov-single-j', help='min coverage with the single J hit (default: 15)',
                           type=int, default=15, metavar='INT', dest='cov_single_j')
    comb_args.add_argument('--cov-mult-j', help='min coverage with multiple J hits (default: 5)',
                           type=int, default=5, metavar='INT', dest='cov_mult_j')
    comb_args.add_argument('--human-readable', help='include human readable info', action='store_true',
                            dest='human_readable')
    human_args += [('Combination detection', None), ('Length', 'length'), ('Coverage w/ single J hit', 'cov_single_j'),
                   ('Coverage w/ multiple J hits', 'cov_mult_j'), ('Human readable', 'human_readable')]


    exp_args = parser.add_argument_group('combination expansion arguments')
    exp_args.add_argument('--comb-coverage', help='combination coverage threshold (default: 50)',
                          type=int, default=50, metavar='INT', dest='comb_coverage')
    exp_args.add_argument('--expansion-rate', help='mismatch coverage threshold (ratio to combination coverage) (default: 0.9)',
                          type=float, default=0.9, metavar='FLOAT', dest='expansion_rate')
    human_args += [('Combination expansion', None), ('Combination coverage', 'comb_coverage'), ('Expansion rate', 'expansion_rate')]

    gen_args = parser.add_argument_group('segments generation arguments')
    gen_args.add_argument('--min-significance', help='min significance (default: 20)',
                          type=int, metavar='INT', default=20, dest='min_significance')
    human_args += [('Segments generation', None), ('Min significance', 'min_significance')]

    other = parser.add_argument_group('other arguments')
    other.add_argument('-h', '--help', action='help', help='show this help message and exit')

    args = parser.parse_args()

    try:
        run(args, human_args)
    except KeyboardInterrupt:
        sys.stderr.write('\nKeyboard interrupt\n')
        exit(1)


if __name__ == '__main__':
    main()

