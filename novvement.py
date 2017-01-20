#!/usr/bin/python3

import argparse
import os
import sys
import io
import subprocess
import time
import datetime

from extra._version import __version__
from extra.utilities import oprint, eprint


def rc_fail():
    eprint('\nNon-empty return code\n', 'red')
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
    oprint('Segment coverage\n', 'green')
    log.write('\n# Segment coverage\n')

    for name, path in datasets:
        with open(os.path.join(path, 'alignment_info.csv')) as inp, \
                open(os.path.join(dir, name, 'alignment', 'segment_coverage.csv'), 'w') as outp:
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
    oprint('J hit\n', 'green')
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

    oprint('V alignment -> csv\n', 'green')
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

    oprint('Detecting potential mismatches\n', 'green')
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


def combinations_dirname(iteration):
    return 'combinations_' + chr(ord('a') + iteration)


def combinations(args, log, datasets, iteration, last):
    dir = args.output
    script_path = os.path.dirname(__file__)

    comb_dir = combinations_dirname(iteration)

    oprint('Combining potential mismatches. ', 'green')
    oprint('Iteration %d\n' % (iteration + 1), 'cyan')
    log.write('\n# Combining potential mismatches. Iteration %d\n' % (iteration + 1))

    for name, path in datasets:
        sys.stdout.write('%s ' % name)
        sys.stdout.flush()
        command = [os.path.join(script_path, 'combinations.py'),
                   '-v', os.path.join(dir, name, 'alignment', 'v_alignment.csv'),
                   '-j', os.path.join(dir, name, 'alignment', 'j_hit.csv'),
                   '-m', os.path.join(dir, name, comb_dir, 'potential_mismatches.csv')]
        if last:
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


def expand_mismatches(args, log, datasets, iteration):
    dir = args.output
    script_path = os.path.dirname(__file__)

    oprint('Expanding combinations. ', 'green')
    oprint('Iteration %d\n' % (iteration + 1), 'cyan')
    log.write('\n# Expanding combinations. Iteration %d\n' % (iteration + 1))

    prev_dir = combinations_dirname(iteration)
    next_dir = combinations_dirname(iteration + 1)

    for name, path in datasets:
        mkdir(os.path.join(dir, name, next_dir))
        sys.stdout.write('%s ' % name)
        sys.stdout.flush()
        command = [os.path.join(script_path, 'expand_potential_mismatches.py'),
                   '-v', os.path.join(dir, name, 'alignment', 'v_alignment.csv'),
                   '-m', os.path.join(dir, name, prev_dir, 'potential_mismatches.csv'),
                   '-c', os.path.join(dir, name, prev_dir, 'combinations.csv'),
                   '-o', os.path.join(dir, name, next_dir, 'potential_mismatches.csv'),
                   '-e', os.path.join(dir, name, next_dir, 'expanded_mismatches.csv'),
                   '--range', args.range[0], args.range[1],
                   '--coverage', args.comb_coverage,
                   '--rate', args.expansion_rate]
        command = [str(x) for x in command]
        log.write('\t%s\n' % ' '.join(command))
        if subprocess.run(command).returncode:
            rc_fail()

    sys.stdout.write('\n')


def significance(args, log, datasets, last_dir):
    dir = args.output

    with open(os.path.join(dir, 'combinations_path.csv'), 'w') as f:
        log.write('\n# generate <name> <combinations.csv> file (%s)\n' % f.name)
        for name, _ in datasets:
            f.write('%s %s\n' % (name, os.path.join(name, last_dir, 'combinations.csv')))

    script_path = os.path.dirname(__file__)

    oprint('Evaluating combination significance\n', 'green')
    log.write('\n# Evaluating combination significance\n')

    command = [os.path.join(script_path, 'significance.py'),
               '-i', os.path.join(dir, 'combinations_path.csv'),
               '-o', os.path.join(dir, 'significance.csv'),
               '-m', args.min_significance]
    command = [str(x) for x in command]
    log.write('\t%s\n' % ' '.join(command))
    subprocess.run(command)

    log.write('\n# Significance -> Combinations\n')
    with open(os.path.join(dir, 'significance.csv')) as inp, open(os.path.join(dir, 'pre_combinations.csv'), 'w') as outp:
        log.write('\tcat %s | cut -f2- > %s\n' % (inp.name, outp.name))
        p1 = subprocess.Popen(['cut', '-f2-'], stdin=inp, stdout=outp)
        p1.communicate()
        if p1.returncode:
            rc_fail()


def filter_combinations(args, log, datasets):
    dir = args.output
    script_path = os.path.dirname(__file__)

    oprint('Filter combinations\n', 'green')
    log.write('\n# Filter combinations\n')

    command = [os.path.join(script_path, 'similarity_filter.py'),
               '-c', os.path.join(dir, 'pre_combinations.csv'),
               '-v', args.v_segments,
               '-o', os.path.join(dir, 'combinations.csv'),
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

    oprint('Generating novel segments\n', 'green')
    log.write('\n# Generating novel segments\n')

    command = [os.path.join(script_path, 'combinations_to_segments.py'),
               '-i', os.path.join(dir, 'combinations.csv'),
               '-v', args.v_segments,
               '-o', os.path.join(dir, 'segments.fa')]
    command = [str(x) for x in command]
    log.write('\t%s\n' % ' '.join(command))
    if subprocess.run(command).returncode:
        rc_fail()

    oprint('Segments generated', 'yellow', 'bold')
    sys.stdout.write(' ::::: ')
    oprint('%s\n' % os.path.join(dir, 'segments.fa'), 'cyan')


def run(args, human_args):
    start = time.perf_counter()
    dir = args.output
    mkdir(args.output)
    if not args.force and os.listdir(dir):
        eprint('Output folder is not empty. Try -f/--force\n', 'red')
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

        for i in range(args.expansion_cycles):
            combinations(args, log, datasets, i, False)
            expand_mismatches(args, log, datasets, i)
        combinations(args, log, datasets, args.expansion_cycles, True)

        significance(args, log, datasets, combinations_dirname(args.expansion_cycles))
        filter_combinations(args, log, datasets)
        generate(args, log, datasets)

        seconds = time.perf_counter() - start
        time_str = str(datetime.timedelta(seconds=round(seconds)))
        oprint('Execution time: %s\n' % time_str, 'magenta')
        log.write('\n# Execution time: %s\n' % t)


def main():
    parser = argparse.ArgumentParser(description='Pipeline for detecting novel segments.\n'
                                                 'Created by Timofey Prodanov (timofey.prodanov@gmail.com)',
                                     usage='%(prog)s -i File -v File -o Dir',
                                     add_help=False, formatter_class=argparse.RawTextHelpFormatter)
    io_args = parser.add_argument_group('Inpur/output arguments')
    io_args.add_argument('-i', '--input', required=True, metavar='File', type=argparse.FileType(),
                         help='File with lines <name> <path to cdr folder>.\n'
                              'Should contain: v_alignments.fa, alignment_info.csv')
    io_args.add_argument('-v', '--v-segments', required=True, metavar='File',
                         help='Fasta file containing V segments', dest='v_segments')
    io_args.add_argument('-o', '--output', help='Output directory', required=True, metavar='Dir')
    io_args.add_argument('-f', '--force', help='Override files in output directory', action='store_true')
    human_args = [('Input/output', None), ('Input', 'input'), ('V segments', 'v_segments'), ('Output', 'output'),
                  ('Force', 'force')]

    mismatch_args = parser.add_argument_group('Mismatch detection arguments')
    mismatch_args.add_argument('--range', help='Positions range (default: [60, 290])',
                               metavar=('Int', 'Int'), nargs=2, default=[60, 290], type=int)
    mismatch_args.add_argument('--segment-coverage', help='Segment coverage threshold (default: 200)',
                               type=int, default=200, metavar='Int', dest='segment_coverage')
    mismatch_args.add_argument('--mismatch-rate', help='Mismatch rate threshold (default: 0.1)',
                               type=float, default=0.1, metavar='Float', dest='mismatch_rate')
    human_args += [('mismatch detection', None), ('Range', 'range'), ('Segment coverage', 'segment_coverage'),
                   ('Mismatch rate', 'mismatch_rate')]

    comb_args = parser.add_argument_group('Combination detection arguments')
    comb_args.add_argument('--length', help='Min combination length (default: 2)',
                           type=int, default=2, metavar='Int')
    comb_args.add_argument('--cov-single-j', help='Min coverage with the single J hit (default: 15)',
                           type=int, default=15, metavar='Int', dest='cov_single_j')
    comb_args.add_argument('--cov-mult-j', help='Min coverage with multiple J hits (default: 5)',
                           type=int, default=5, metavar='Int', dest='cov_mult_j')
    comb_args.add_argument('--human-readable', help='Include human readable info', action='store_true',
                            dest='human_readable')
    human_args += [('Combination detection', None), ('Length', 'length'), ('Coverage w/ single J hit', 'cov_single_j'),
                   ('Coverage w/ multiple J hits', 'cov_mult_j'), ('Human readable', 'human_readable')]

    exp_args = parser.add_argument_group('Combination expansion arguments')
    exp_args.add_argument('--expansion-cycles', help='Number of combination expansion\n'
                                                     'cycles (default: 1)',
                          type=int, default=1, metavar='Int', dest='expansion_cycles')                    
    exp_args.add_argument('--comb-coverage', help='Combination coverage threshold (default: 50)',
                          type=int, default=50, metavar='Int', dest='comb_coverage')
    exp_args.add_argument('--expansion-rate', help='Mismatch coverage threshold (ratio\n'
                                                   'to combination coverage) (default: 0.9)',
                          type=float, default=0.9, metavar='Float', dest='expansion_rate')
    human_args += [('Combination expansion', None), ('Expansion cycles', 'expansion_cycles'),
                   ('Combination coverage', 'comb_coverage'), ('Expansion rate', 'expansion_rate')]

    filter_args = parser.add_argument_group('Segments filtering arguments')
    filter_args.add_argument('--min-significance', help='Min significance (default: 20)',
                             type=int, metavar='Int', default=20, dest='min_significance')
    filter_args.add_argument('--source-dist', help='Minimum required distance to\n'
                                                   'any source segment (default: --length)',
                             metavar='Int', dest='source_dist', type=int)
    filter_args.add_argument('--target-dist', help='Minimum reliable distance to\n'
                                                   'any target combination (default: 3)',
                             metavar='Int', dest='target_dist', type=int, default=3)
    filter_args.add_argument('--target-mf', help='Significance multiplication factor\n'
                                                 'to filter out combination with unreliable\n'
                                                 'target distance (default: 4)',
                             metavar='Float', dest='target_mf', type=float, default=4)
    human_args += [('Segments filtering', None), ('Min significance', 'min_significance'),
                   ('Distance to source', 'source_dist'), ('Distance to target', 'target_dist'),
                   ('Significance multiplication factor', 'target_mf')]

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other.add_argument('--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()
    if not args.source_dist:
        args.source_dist = args.length

    try:
        run(args, human_args)
    except KeyboardInterrupt:
        eprint('\nKeyboard interrupt\n', 'red')
        exit(1)


if __name__ == '__main__':
    main()
