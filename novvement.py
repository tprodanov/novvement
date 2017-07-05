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
            if isinstance(v, list) or isinstance(v, tuple):
                v = ', '.join(str(x) for x in v)
            elif isinstance(v, io.IOBase):
                v = v.name
            log.write('#\t\t%s: %s\n' % (k, v))

    log.write('#\n# Start: %s\n#\n' % datetime.datetime.now().strftime('%d-%m-%Y %H:%M:%S')) 

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
                   '-i', os.path.join(path, args.v_align_name),
                   '-o', os.path.join(dir, 'data', name, 'alignment', 'v_alignment.csv'),
                   '--mismatches-only']
        log.write('\t%s\n' % ' '.join(command))
        if subprocess.run(command).returncode:
            rc_fail()

    sys.stdout.write('\n')


def potential_mismatches(args, log, datasets):
    dir = args.output
    script_path = os.path.dirname(__file__)

    oprint('Detecting potential polymorphisms\n', 'green')
    log.write('\n# Detecting potential polymorphisms\n')
    for name, path in datasets:
        sys.stdout.write('%s ' % name)
        sys.stdout.flush()
        command = [os.path.join(script_path, 'potential_polymorphisms.py'),
                   '-v', os.path.join(dir, 'data', name, 'alignment', 'v_alignment.csv'),
                   '-s', os.path.join(dir, 'data', name, 'alignment', 'segment_coverage.csv'),
                   '-o', os.path.join(dir, 'data', name, 'alignment', 'potential.csv'),
                   '--range', args.range[0], args.range[1],
                   '--coverage', args.segment_coverage,
                   '--bad', args.mismatch_rate[0],
                   '--good', args.mismatch_rate[1],
                   '--bad-cov', args.bad_coverage,
                   '--mult', args.mult,
                   '--nth', args.nth,
                   '--radius', args.radius]
        command = [str(x) for x in command]
        log.write('\t%s\n' % ' '.join(command))
        if subprocess.run(command).returncode:
            rc_fail()

    sys.stdout.write('\n')


def working_dirname(i):
    from string import ascii_lowercase

    if i == -1:
        return 'alignment'

    n = len(ascii_lowercase)

    if i < n:
        suffix = ascii_lowercase[i]
    else:
        suffix = ascii_lowercase[i // n - 1] + ascii_lowercase[i % n]
    return 'combinations_' + suffix


def combinations(args, log, datasets, iteration, last):
    dir = args.output
    script_path = os.path.dirname(__file__)

    prev_dir = working_dirname(iteration - 1)
    comb_dir = working_dirname(iteration)

    oprint('Combining potential polymorphisms. ', 'green')
    oprint('Iteration %d\n' % (iteration + 1), 'cyan')
    log.write('\n# Combining potential polymorphisms. Iteration %d\n' % (iteration + 1))

    for name, path in datasets:
        mkdir(os.path.join(dir, 'data', name, comb_dir))
        sys.stdout.write('%s ' % name)
        sys.stdout.flush()
        command = [os.path.join(script_path, 'combinations.py'),
                   '-v', os.path.join(dir, 'data', name, 'alignment', 'v_alignment.csv'),
                   '-j', os.path.join(dir, 'data', name, 'alignment', 'j_hit.csv'),
                   '-m', os.path.join(dir, 'data', name, prev_dir, 'potential.csv')]
        if last:
            command += ['-c', os.path.join(dir, 'data', name, comb_dir, 'combinations.csv'),
                        '--detection-cov', 0,
                        '--length', args.length]
        else:
            command += ['-p', os.path.join(dir, 'data', name, comb_dir, 'potential.csv'),
                        '-e', os.path.join(dir, 'data', name, comb_dir, 'expanded.csv'),
                        '--range', args.range[0], args.range[1],
                        '--expansion-cov', args.expansion_coverage,
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
               '-v', *args.v_segments,
               '-o', os.path.join(dir, 'group', 'short.txt'),
               '-f', os.path.join(dir, 'group', 'full.txt'),
               '-c', os.path.join(dir, 'group', 'components.txt'),
               '--range', args.range[0], args.range[1],
               '--source-dist', args.length,
               '--j-coverage', args.j_coverage,
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
        outp.write('segment\tsignificance\tlength\tcombination\n')
        outp.flush()

        grep_command = ('grep', '-Pv', '^([\\s#]|$)')
        cut_command = ('cut', '-f2,5')
        sed_command = ('sed', 's/(\\|)/\\t/g')
        awk_command = ('awk', '$1 >= %d {OFS="\\t"; t=$1; $1=$2; $2=t;  print}' % args.min_significance)

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
               '-v', *args.v_segments,
               '-o', os.path.join(dir, 'segments.fa')]
    command = [str(x) for x in command]
    log.write('\t%s\n' % ' '.join(command))
    if subprocess.run(command).returncode:
        rc_fail()

    oprint('Segments generated', 'yellow', 'bold')
    sys.stdout.write(' ::::: ')
    oprint('%s\n' % os.path.join(dir, 'segments.fa'), 'cyan')


def run(args, human_args):
    if args.range[0] <= 0:
        eprint('Invalid left margin: --range %d %d\n' % tuple(args.range), 'red')
        exit(1)

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
        for i in range(args.cycles):
            combinations(args, log, datasets, i, i == args.cycles - 1)

        group_allelic(args, log, datasets, working_dirname(args.cycles - 1))
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
    io_args.add_argument('-v', '--v-segments', required=True, metavar='File', nargs='+',
                         help='Fasta file containing V segments (can be several files)')
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
    input_fmt.add_argument('--v-align-name', metavar='Str', default='v_alignments.fa',
                           help='Name of a V alignment file (default: v_alignments.fa)')
    human_args += [('Input format', None), ('V column', 'v_col'), ('J column', 'j_col'),
                   ('V alignment file name', 'v_align_name')]

    mismatch_args = parser.add_argument_group('Polymorphisms detection')
    mismatch_args.add_argument('--range', help='Positions range (default: 40, 290)',
                               metavar=('Int', 'Int'), nargs=2, default=[40, 290], type=int)
    mismatch_args.add_argument('--segment-coverage', help='Segment coverage threshold (default: 1)',
                               type=int, default=1, metavar='Int')
    mismatch_args.add_argument('--mismatch-rate', nargs=2, metavar='Float', default=[0.1, 0.2], type=float,
                               help='Two values: <bad-rate> <good-rate>\n'
                                    'mismatches with rate lower than <bad-rate> are discarded\n'
                                    'mismatches with rate higher than <good-rate> are selected\n'
                                    '(default: 0.1 0.2)')
    mismatch_args.add_argument('--bad-coverage', metavar='Int', type=int, default=30,
                               help='Polymorphisms with coverage less than <--bad-coverage>\n'
                                    'are discarded (default: 30)')
    mismatch_args.add_argument('--mult', metavar='Float', type=float, default=7,
                               help='Polymorphisms with rate between <bad-rate> <good-rate>\n'
                                    'should appear <--mult> times more frequent than the <--nth>\n'
                                    'frequent polymorphism in the neighborhood (default: 7)')
    mismatch_args.add_argument('--nth', metavar='Int', type=int, default=2,
                               help='Nth frequent polymorphism in the neighborhood for comparison\n'
                                    '0 is the most frequent (will never trigger)\n'
                                    '<--radius> : median\n'
                                    '2 * <--radius> : the least frequent\n'
                                    '(default: 2)')
    mismatch_args.add_argument('--radius', metavar='Int', type=int, default=5,
                              help='Neighborhood radius (default: 5)')
    human_args += [('Polymorphisms detection', None), ('Range', 'range'),
                   ('Segment coverage', 'segment_coverage'), ('Mismatch rate', 'mismatch_rate'),
                   ('Bad coverage', 'bad_coverage'), ('Multiplication', 'mult'),
                   ('Nth neighbor', 'nth'), ('Neighborhood radius', 'radius')]

    comb_args = parser.add_argument_group('Combination detection and expansion')
    comb_args.add_argument('--length', help='Min difference from existing segment (default: 4)',
                           type=int, default=4, metavar='Int')
    comb_args.add_argument('--cycles', type=int, default=4, metavar='Int',
                           help='Number of combination detection and expansion cycles (default: 4)')
    comb_args.add_argument('--expansion-coverage',
                           help='Min combination coverage during expansion (default: 50)',
                           type=int, default=50, metavar='Int')
    comb_args.add_argument('--expansion-rate', type=float, default=0.75, metavar='Float',
                           help='Polymorphism should appear together with combination more often\n'
                           'than <--expansion-rate> * <combination coverage> (default: 0.75)')
    human_args += [('Combination detection and expansion', None), ('Length', 'length'), ('Cycles', 'cycles'),
                   ('Expansion coverage', 'expansion_coverage'), ('Expansion rate', 'expansion_rate')]

    filter_args = parser.add_argument_group('Segments filtering and grouping')
    filter_args.add_argument('--min-significance', help='Min significance (default: 3)',
                             type=int, metavar='Float', default=3, dest='min_significance')
    filter_args.add_argument('--j-coverage', metavar='Int', type=int, default=15,
                             help='J segment covers novel segment if they appear\n'
                                  'together at least <--j-coverage> times (default: 15)')
    filter_args.add_argument('--hamming', metavar='Int', type=int, default=0,
                             help='If distance between two novel segments is at most\n'
                                  '<--hamming> - these segments will be placed in\n'
                                  'one group (default: 0)')
    filter_args.add_argument('--shared', metavar='Float', type=float, default=0.75,
                             help='Two novel segments will be placed in one group\n'
                                  'if they share <--shared> of their\n'
                                  'polymorphisms (default: 0.75)')
    human_args += [('Segments filtering and grouping', None), ('Min significance', 'min_significance'),
                   ('J coverage', 'j_coverage'),
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
