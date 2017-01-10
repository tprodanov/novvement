#!/usr/bin/python3

import argparse
import os
import sys
import io
import subprocess
from _version import __version__


def rc_fail():
    sys.stderr.write('\nNon-empty return code\n')
    exit(1)


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
    log.flush()


def make_datasets_file(dir, log, datasets):
    with open(os.path.join(dir, 'datasets.csv'), 'w') as f:
        log.write('\n# generate <name> <v_alignment.csv> file (%s)\n' % f.name)
        for name, _ in datasets:
            f.write('%s %s\n' % (name, os.path.join(name, 'alignment', 'v_alignment.csv')))


def segment_coverage(args, log, datasets):
    dir = args.output
    sys.stdout.write('Segment coverage\n')
    log.write('\n# Segment coverage\n')

    for name, path in datasets:
        with open(os.path.join(path, 'alignment_info.csv')) as inp, open(os.path.join(dir, name, 'alignment', 'segment_coverage.csv'), 'w') as outp:
            log.write("\tcat %s | cut -f5 | sort | uniq -c | awk '{t = $1; $1 = $2; $2 = t; print}' > %s\n" % (inp.name, outp.name))
            p1 = subprocess.Popen(['cut', '-f5'], stdin=inp, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(['sort'], stdin=p1.stdout, stdout=subprocess.PIPE)
            p3 = subprocess.Popen(['uniq', '-c'], stdin=p2.stdout, stdout=subprocess.PIPE)
            p4 = subprocess.Popen(['awk', '{t = $1; $1 = $2; $2 = t; print}'], stdin=p3.stdout, stdout=outp)
            p1.stdout.close()
            p2.stdout.close()
            p3.stdout.close()
            p4.communicate()
            if p4.returncode:
                rc_fail()


def j_hits(args, log, datasets):
    dir = args.output
    sys.stdout.write('J hit\n')
    log.write('\n# J hit\n')

    for name, path in datasets:
        with open(os.path.join(path, 'alignment_info.csv')) as inp, open(os.path.join(dir, name, 'alignment', 'j_hit.csv'), 'w') as outp:
            log.write("\tcat %s | cut -f1,9 > %s\n" % (inp.name, outp.name))
            p1 = subprocess.Popen(['cut', '-f1,9'], stdin=inp, stdout=outp)
            p1.communicate()
            if p1.returncode:
                rc_fail()


def v_alignment(args, log, datasets):
    dir = args.output
    script_path = os.path.dirname(__file__)

    sys.stdout.write('V alignment -> csv\n')
    log.write('\n# V alignment -> csv\n')
    for name, path in datasets:
        sys.stdout.write('%s ' % name)
        sys.stdout.flush()
        command = [os.path.join(script_path, 'v_alignment_mismatches.py'),
                   '-i', os.path.join(path, 'v_alignments.fa'),
                   '-o', os.path.join(dir, name, 'alignment', 'v_alignment.csv')]
        log.write('\t%s\n' % ' '.join(command))
        if subprocess.run(command).returncode:
            rc_fail()

    sys.stdout.write('\n')


def potential_mismatches(args, log, datasets):
    dir = args.output
    script_path = os.path.dirname(__file__)

    sys.stdout.write('Detecting potential mismatches\n')
    log.write('\n# Detecting potential mismatches\n')
    for name, path in datasets:
        mkdir(os.path.join(dir, name, 'combinations_a'))
        sys.stdout.write('%s ' % name)
        sys.stdout.flush()
        command = [os.path.join(script_path, 'potential_mismatches.py'),
                   '-v', os.path.join(dir, name, 'alignment', 'v_alignment.csv'),
                   '-s', os.path.join(dir, name, 'alignment', 'segment_coverage.csv'),
                   '-o', os.path.join(dir, name, 'combinations_a', 'potential_mismatches.csv'),
                   '--range', args.range[0], args.range[1],
                   '--coverage', args.segment_coverage,
                   '--rate', args.mismatch_rate]
        command = [str(x) for x in command]
        log.write('\t%s\n' % ' '.join(command))
        if subprocess.run(command).returncode:
            rc_fail()

    sys.stdout.write('\n')


def combinations(args, log, datasets, comb_dir):
    dir = args.output
    script_path = os.path.dirname(__file__)

    sys.stdout.write('Combining potential mismatches\n')
    log.write('\n# Combining potential mismatches\n')

    for name, path in datasets:
        sys.stdout.write('%s ' % name)
        sys.stdout.flush()
        command = [os.path.join(script_path, 'combinations.py'),
                   '-v', os.path.join(dir, name, 'alignment', 'v_alignment.csv'),
                   '-j', os.path.join(dir, name, 'alignment', 'j_hit.csv'),
                   '-m', os.path.join(dir, name, comb_dir, 'potential_mismatches.csv')]
        if comb_dir.endswith('b'):
            command += ['-d', os.path.join(dir, 'datasets.csv')]
        command += ['-o', os.path.join(dir, name, comb_dir, 'combinations.csv')]
        if args.human_readable:
            command += ['-h', os.path.join(dir, name, comb_dir, 'human_combinations.txt')]
        command += ['--length', args.length,
                    '--cov-single-j', args.cov_single_j, '--cov-mult-j', args.cov_mult_j]

        command = [str(x) for x in command]
        log.write('\t%s\n' % ' '.join(command))
        if subprocess.run(command).returncode:
            rc_fail()

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
        command = [os.path.join(script_path, 'expand_potential_mismatches.py'),
                   '-v', os.path.join(dir, name, 'alignment', 'v_alignment.csv'),
                   '-m', os.path.join(dir, name, 'combinations_a', 'potential_mismatches.csv'),
                   '-c', os.path.join(dir, name, 'combinations_a', 'combinations.csv'),
                   '-o', os.path.join(dir, name, 'combinations_b', 'potential_mismatches.csv'),
                   '-e', os.path.join(dir, name, 'combinations_b', 'expanded_mismatches.csv'),
                   '--range', args.range[0], args.range[1],
                   '--coverage', args.comb_coverage,
                   '--rate', args.expansion_rate]
        command = [str(x) for x in command]
        log.write('\t%s\n' % ' '.join(command))
        if subprocess.run(command).returncode:
            rc_fail()

    sys.stdout.write('\n')


def significance(args, log, datasets):
    dir = args.output

    with open(os.path.join(dir, 'combinations_path.csv'), 'w') as f:
        log.write('\n# generate <name> <combinations.csv> file (%s)\n' % f.name)
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
        if p1.returncode:
            rc_fail()


def filter_combinations(args, log, datasets):
    dir = args.output
    script_path = os.path.dirname(__file__)

    sys.stdout.write('Filter combinations\n')
    log.write('\n# Filter combinations\n')

    command = [os.path.join(script_path, 'similarity_filter.py'),
               '-c', os.path.join(dir, 'combinations.csv'),
               '-v', args.v_segments,
               '-o', os.path.join(dir, 'f_combinations.csv'),
               '-l', os.path.join(dir, 'filter.log'),
               '--range', args.range[0], args.range[1],
               '--source-dist', args.source_dist,
               '--target-dist', args.target_dist,
               '--target-mf', args.target_mf]
    command = [str(x) for x in command]
    log.write('\t%s\n' % ' '.join(command))
    if subprocess.run(command).returncode:
        rc_fail()


def generate(args, log, datasets):
    dir = args.output
    script_path = os.path.dirname(__file__)

    sys.stdout.write('Generating novel segments\n')
    log.write('\n# Generating novel segments\n')

    command = [os.path.join(script_path, 'combinations_to_segments.py'),
               '-c', os.path.join(dir, 'f_combinations.csv'),
               '-v', args.v_segments,
               '-o', os.path.join(dir, 'segments.fa')]
    command = [str(x) for x in command]
    log.write('\t%s\n' % ' '.join(command))
    if subprocess.run(command).returncode:
        rc_fail()

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

        segment_coverage(args, log, datasets)
        j_hits(args, log, datasets)
        v_alignment(args, log, datasets)
        potential_mismatches(args, log, datasets)
        combinations(args, log, datasets, 'combinations_a')
        expand_mismatches(args, log, datasets)
        combinations(args, log, datasets, 'combinations_b')
        significance(args, log, datasets)
        filter_combinations(args, log, datasets)
        generate(args, log, datasets)


def main():
    parser = argparse.ArgumentParser(description='Pipeline for detecting novel segments.\n'
                                                 'Created by Timofey Prodanov (timofey.prodanov@gmail.com)',
                                     usage='%(prog)s -i DATASETS -v V_SEGMENTS -o OUT_DIR',
                                     add_help=False, formatter_class=argparse.RawTextHelpFormatter)
    io_args = parser.add_argument_group('Inpur/output arguments')
    io_args.add_argument('-i', '--input', required=True, metavar='FILE', type=argparse.FileType(),
                         help='file with lines <name> <path to cdr folder>.\n'
                              'Should contain: v_alignments.fa, alignment_info.csv')
    io_args.add_argument('-v', '--v-segments', required=True, metavar='FILE', help='file containing V segments',
                         dest='v_segments')
    io_args.add_argument('-o', '--output', help='output directory', required=True, metavar='DIR')
    io_args.add_argument('-f', '--force', help='override files in output directory', action='store_true')
    human_args = [('Input/output', None), ('Input', 'input'), ('V segments', 'v_segments'), ('Output', 'output'), ('Force', 'force')]

    mismatch_args = parser.add_argument_group('Mismatch detection arguments')
    mismatch_args.add_argument('--range', help='positions range (default: [60, 290])',
                               metavar=('INT', 'INT'), nargs=2, default=[60, 290], type=int)
    mismatch_args.add_argument('--segment-coverage', help='segment coverage threshold (default: 200)',
                               type=int, default=200, metavar='INT', dest='segment_coverage')
    mismatch_args.add_argument('--mismatch-rate', help='mismatch rate threshold (default: 0.1)',
                               type=float, default=0.1, metavar='FLOAT', dest='mismatch_rate')
    human_args += [('mismatch detection', None), ('Range', 'range'), ('Segment coverage', 'segment_coverage'), ('Mismatch rate', 'mismatch_rate')]

    comb_args = parser.add_argument_group('Combination detection arguments')
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

    exp_args = parser.add_argument_group('Combination expansion arguments')
    exp_args.add_argument('--comb-coverage', help='combination coverage threshold (default: 50)',
                          type=int, default=50, metavar='INT', dest='comb_coverage')
    exp_args.add_argument('--expansion-rate', help='mismatch coverage threshold (ratio\nto combination coverage) (default: 0.9)',
                          type=float, default=0.9, metavar='FLOAT', dest='expansion_rate')
    human_args += [('Combination expansion', None), ('Combination coverage', 'comb_coverage'), ('Expansion rate', 'expansion_rate')]

    filter_args = parser.add_argument_group('Segments filtering arguments')
    filter_args.add_argument('--min-significance', help='min significance (default: 20)',
                             type=int, metavar='INT', default=20, dest='min_significance')
    filter_args.add_argument('--source-dist', help='minimum required distance to\n'
                                                   'any source segment (default: --length)',
                             metavar='INT', dest='source_dist', type=int)
    filter_args.add_argument('--target-dist', help='minimum reliable distance to\n'
                                                   'any target combination (default: 3)',
                             metavar='INT', dest='target_dist', type=int, default=3)
    filter_args.add_argument('--target-mf', help='significance multiplication factor\n'
                                                 'to filter out combination with unreliable\n'
                                                 'target distance (default: 4)',
                             metavar='FLOAT', dest='target_mf', type=float, default=4)
    human_args += [('Segments filtering', None), ('Min significance', 'min_significance'),
                   ('Distance to source', 'source_dist'), ('Distance to target', 'target_dist'),
                   ('Significance multiplication factor', 'target_mf')]

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='show this help message and exit')
    other.add_argument('--version', action='version', help='show version',
                       version=__version__)

    args = parser.parse_args()
    if not args.source_dist:
        args.source_dist = args.length

    try:
        run(args, human_args)
    except KeyboardInterrupt:
        sys.stderr.write('\nKeyboard interrupt\n')
        exit(1)


if __name__ == '__main__':
    main()
