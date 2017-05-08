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

    for line in __version__.splitlines():
        log.write('# %s\n' % line)
    log.write('\n')
    log.flush()


def make_datasets_file(dir, log, datasets):
    with open(os.path.join(dir, 'data', 'datasets.csv'), 'w') as f:
        log.write('\n# generate <name> <v_alignment.csv> file (%s)\n' % f.name)
        for name, _ in datasets:
            f.write('%s %s\n' % (name, os.path.join('data', name, 'alignment', 'v_alignment.csv')))


def segment_coverage(args, log, datasets):
    dir = args.output
    oprint('Segment coverage\n', 'green')
    log.write('\n# Segment coverage\n')

    for name, path in datasets:
        with open(os.path.join(path, 'alignment_info.csv')) as inp, \
                open(os.path.join(dir, 'data', name, 'alignment', 'segment_coverage.csv'), 'w') as outp:
            log.write("\tcat %s | cut -f%d | sort | uniq -c | awk '{t = $1; $1 = $2; $2 = t; print}' > %s\n" 
                      % (inp.name, args.v_col, outp.name))
            p1 = subprocess.Popen(['cut', '-f', str(args.v_col)], stdin=inp, stdout=subprocess.PIPE)
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
        with open(os.path.join(path, 'alignment_info.csv')) as inp, \
                open(os.path.join(dir, 'data', name, 'alignment', 'j_hit.csv'), 'w') as outp:
            log.write("\tcat %s | cut -f1,%d > %s\n" % (inp.name, args.j_col, outp.name))
            p1 = subprocess.Popen(['cut', '-f1,%d' % args.j_col], stdin=inp, stdout=outp)
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
                   '-o', os.path.join(dir, 'data', name, 'alignment', 'v_alignment.csv'),
                   '--mismatches-only']
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
        mkdir(os.path.join(dir, 'data', name, 'combinations_a'))
        sys.stdout.write('%s ' % name)
        sys.stdout.flush()
        command = [os.path.join(script_path, 'potential_mismatches.py'),
                   '-v', os.path.join(dir, 'data', name, 'alignment', 'v_alignment.csv'),
                   '-s', os.path.join(dir, 'data', name, 'alignment', 'segment_coverage.csv'),
                   '-o', os.path.join(dir, 'data', name, 'combinations_a', 'potential_mismatches.csv'),
                   '--range', args.range[0], args.range[1],
                   '--coverage', args.segment_coverage,
                   '--rate', args.mismatch_rate]
        command = [str(x) for x in command]
        log.write('\t%s\n' % ' '.join(command))
        if subprocess.run(command).returncode:
            rc_fail()

    sys.stdout.write('\n')


def combinations_dirname(i):
    from string import ascii_lowercase
    n = len(ascii_lowercase)

    if i < n:
        suffix = ascii_lowercase[i]
    else:
        suffix = ascii_lowercase[i // n - 1] + ascii_lowercase[i % n]
    return 'combinations_' + suffix


def combinations(args, log, datasets, iteration, detection_coverage):
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
                   '-v', os.path.join(dir, 'data', name, 'alignment', 'v_alignment.csv'),
                   '-j', os.path.join(dir, 'data', name, 'alignment', 'j_hit.csv'),
                   '-m', os.path.join(dir, 'data', name, comb_dir, 'potential_mismatches.csv'),
                   '-o', os.path.join(dir, 'data', name, comb_dir, 'combinations.csv'),
                   '--length', args.length,
                    '--cov-single-j', detection_coverage, '--cov-mult-j', detection_coverage]

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
        mkdir(os.path.join(dir, 'data', name, next_dir))
        sys.stdout.write('%s ' % name)
        sys.stdout.flush()
        command = [os.path.join(script_path, 'expand_potential_mismatches.py'),
                   '-v', os.path.join(dir, 'data', name, 'alignment', 'v_alignment.csv'),
                   '-m', os.path.join(dir, 'data', name, prev_dir, 'potential_mismatches.csv'),
                   '-c', os.path.join(dir, 'data', name, prev_dir, 'combinations.csv'),
                   '-o', os.path.join(dir, 'data', name, next_dir, 'potential_mismatches.csv'),
                   '-e', os.path.join(dir, 'data', name, next_dir, 'expanded_mismatches.csv'),
                   '--range', args.range[0], args.range[1],
                   '--coverage', args.expansion_coverage,
                   '--rate', args.expansion_rate]
        command = [str(x) for x in command]
        log.write('\t%s\n' % ' '.join(command))
        if subprocess.run(command).returncode:
            rc_fail()

    sys.stdout.write('\n')


def group_allelic(args, log, datasets, last_dir):
    dir = args.output
    script_path = os.path.dirname(__file__)

    oprint('Grouping allelic variations\n', 'green')

    mkdir(os.path.join(dir, 'group'))
    filename = os.path.join(dir, 'group', 'paths.txt')
    log.write('\n# Writing paths to combinations -> %s\n' % filename)
    with open(filename, 'w') as f:
        for name, _ in datasets:
            f.write('%s\t%s\n' % (name,
                os.path.join('..', 'data', name, last_dir, 'combinations.csv')))

    log.write('\n# Grouping allelic variations\n')
    command = [os.path.join(script_path, 'group_allelic.py'),
               '-d', filename,
               '-v', args.v_segments,
               '-o', os.path.join(dir, 'group', 'short.txt'),
               '-f', os.path.join(dir, 'group', 'full.txt'),
               '-c', os.path.join(dir, 'group', 'components.txt'),
               '--range', args.range[0], args.range[1],
               '--source-dist', args.length,
               '--j-coverage', args.j_coverage,
               '--dat-coverage', args.dat_coverage,
               '--hamming', args.hamming,
               '--shared', args.shared]
    command = [str(x) for x in command]
    log.write('\t%s\n' % ' '.join(command))
    if subprocess.run(command).returncode:
        rc_fail()

    sys.stdout.write('\n')


def generate(args, log, datasets):
    dir = args.output

    script_path = os.path.dirname(__file__)

    oprint('Generating segments\n', 'green')
    inp_name = os.path.join(dir, 'group', 'components.txt')
    out_name = os.path.join(dir, 'combinations.csv')
    log.write('\n# Filtered combinations at %s\n' % out_name)

    with open(inp_name) as inp, open(out_name, 'w') as outp:
        outp.write('segment\tcombination\tsignificance\n')
        outp.flush()

        grep_command = ('grep', '-Pv', '^([\\s#]|$)')
        cut_command = ('cut', '-f2,5')
        sed_command = ('sed', 's/\\./\\t/')
        awk_command = ('awk', '$1 > 2 {OFS="\\t"; t=$1; $1=$2; $2=$3; $3=t; print}')

        log.write('\tcat %s | %s | %s | %s | %s > %s\n'
                  % (inp_name, "%s %s '%s'" % grep_command, ' '.join(cut_command),
                     "%s '%s'" % sed_command, "%s '%s'" % awk_command, out_name))
        
        p1 = subprocess.Popen(grep_command, stdin=inp, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(cut_command, stdin=p1.stdout, stdout=subprocess.PIPE)
        p3 = subprocess.Popen(sed_command, stdin=p2.stdout, stdout=subprocess.PIPE)
        p4 = subprocess.Popen(awk_command, stdin=p3.stdout, stdout=outp)
        p1.stdout.close()
        p2.stdout.close()
        p3.stdout.close()
        p4.communicate()
        if p4.returncode:
            rc_fail()

    log.write('\n# Generating novel segments\n')

    command = [os.path.join(script_path, 'make_segments.py'),
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
    if not args.force and os.listdir(dir) and args.start == 'A':
        eprint('Output folder is not empty. Try -f/--force\n', 'red')
        exit(1)

    datasets = [tuple(line.strip().split(' ', 1)) for line in args.input if not line.startswith('#')]
    mkdir(os.path.join(dir, 'data'))

    with open(os.path.join(dir, 'log.txt'), 'w') as log:
        log_header(log, args, human_args)
        
        if args.start == 'A':
            make_datasets_file(dir, log, datasets)
            for name, path in datasets:
                mkdir(os.path.join(dir, 'data', name))
                mkdir(os.path.join(dir, 'data', name, 'alignment'))

            segment_coverage(args, log, datasets)
            j_hits(args, log, datasets)
            v_alignment(args, log, datasets)
        
        potential_mismatches(args, log, datasets)

        for i in range(args.expansion_cycles):
            combinations(args, log, datasets, i, args.detection_coverage)
            expand_mismatches(args, log, datasets, i)
        combinations(args, log, datasets, args.expansion_cycles, 0)

        group_allelic(args, log, datasets, combinations_dirname(args.expansion_cycles))
        generate(args, log, datasets)

        seconds = time.perf_counter() - start
        t = str(datetime.timedelta(seconds=round(seconds)))
        oprint('Execution time: %s\n' % t, 'magenta')
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
                         help='Fasta file containing V segments')
    io_args.add_argument('-o', '--output', help='Output directory', required=True, metavar='Dir')
    io_args.add_argument('-f', '--force', help='Override files in output directory', action='store_true')
    io_args.add_argument('--start', help='Start from a specific point in a pipeline:\n'
                                         'A: start from the beginning (default)\n'
                                         'B: start from detecting potential mismatches\n',
                         choices=['A', 'B'], default='A')
    human_args = [('Input/output', None), ('Input', 'input'), ('V segments', 'v_segments'), ('Output', 'output'),
                  ('Force', 'force'), ('Start', 'start')]

    input_fmt = parser.add_argument_group('Input format')
    input_fmt.add_argument('--v-col', help='Number of a V-hit column in \n'
                                           '<alignment_info.csv> (default: 5)',
                           type=int, metavar='Int', default=5)
    input_fmt.add_argument('--j-col', help='Number of a J-hit column in \n'
                                           '<alignment_info.csv> (default: 9)',
                           type=int, metavar='Int', default=9)
    human_args += [('Input format', None), ('V column', 'v_col'), ('J column', 'j_col')]

    mismatch_args = parser.add_argument_group('Mismatch detection')
    mismatch_args.add_argument('--range', help='Positions range (default: [60, 290])',
                               metavar=('Int', 'Int'), nargs=2, default=[60, 290], type=int)
    mismatch_args.add_argument('--segment-coverage', help='Segment coverage threshold (default: 200)',
                               type=int, default=200, metavar='Int')
    mismatch_args.add_argument('--mismatch-rate', help='Mismatch rate threshold (default: 0.1)',
                               type=float, default=0.1, metavar='Float')
    human_args += [('Mismatch detection', None), ('Range', 'range'), ('Segment coverage', 'segment_coverage'),
                   ('Mismatch rate', 'mismatch_rate')]

    comb_args = parser.add_argument_group('Combination detection and expansion')
    comb_args.add_argument('--length', help='Min difference from existing segment (default: 4)',
                           type=int, default=4, metavar='Int')
    comb_args.add_argument('--detection-coverage', help='Min coverage to detect a combination (default: 15)',
                           type=int, default=15, metavar='Int')
    comb_args.add_argument('--expansion-cycles', help='Number of combination expansion\n'
                                                     'cycles (default: 1)',
                          type=int, default=1, metavar='Int')
    comb_args.add_argument('--expansion-coverage',
                          help='Min combination coverage during expansion (default: 50)',
                          type=int, default=50, metavar='Int')
    comb_args.add_argument('--expansion-rate', help='Mismatch coverage threshold (ratio\n'
                                                   'to combination coverage) (default: 0.9)',
                          type=float, default=0.9, metavar='Float')
    human_args += [('Combination detection and expansion', None), ('Length', 'length'),
                   ('Detection coverage', 'detection_coverage'), ('Expansion cycles', 'expansion_cycles'),
                   ('Expansion coverage', 'expansion_coverage'), ('Expansion rate', 'expansion_rate')]

    filter_args = parser.add_argument_group('Segments filtering and grouping')
    filter_args.add_argument('--min-significance', help='Min significance (default: 20)',
                             type=int, metavar='Float', default=20, dest='min_significance')
    filter_args.add_argument('--j-coverage', metavar='Int', type=int, default=50,
                             help='J segment covers novel segment if they appear\n'
                                  'together at least <--j-coverage> times (default: 50)')
    filter_args.add_argument('--dat-coverage', metavar='Int', type=int, default=50,
                             help='Novel segment is covered by a dataset if they appear\n'
                                  'together at least <--dat-coverage> times (default: 50)')
    filter_args.add_argument('--hamming', metavar='Int', type=int, default=2,
                             help='If distance between two novel segments is at most\n'
                                  '<--hamming> - these segments will be placed in\n'
                                  'one group (default: 2)')
    filter_args.add_argument('--shared', metavar='Float', type=float, default=0.9,
                             help='Two novel segments will be placed in one group\n'
                                  'if they share <--shared> of their\n'
                                  'polymorphisms (default: 0.9)')
    human_args += [('Segments filtering and grouping', None), ('Min significance', 'min_significance'),
                   ('J coverage', 'j_coverage'), ('Dataset coverage', 'dat_coverage'),
                   ('Hamming distance', 'hamming'), ('Shared polymorphisms', 'shared')]

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other.add_argument('--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    try:
        run(args, human_args)
    except KeyboardInterrupt:
        eprint('\nKeyboard interrupt\n', 'red')
        exit(1)


if __name__ == '__main__':
    main()
