#!/usr/bin/python3

import argparse
import sys
import os
import subprocess
import io
from operator import itemgetter

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


def load_datasets(f, log):
    datasets = []
    i = 0
    for line in f:
        if line.startswith('#'):
            continue
        line = line.strip()
        if not line:
            continue

        s = line.split(' ', 1)
        if len(s) == 1:
            datasets.append((s[0], []))
        elif not datasets:
            sys.stderr.write('Datasets file should start with line <name>, that indicate '
                             'the beginning of a new individual.\n')
            exit(1)
        else:
            datasets[-1][1].append(tuple(s))
            i += 1

    log.write('# Loaded %d datasets within %d individuals from %s\n\n'
              % (i, len(datasets), f.name))
    return datasets


def validate(args, datasets, log):
    dir = args.output
    script_path = os.path.dirname(__file__)
    datasets_path = os.path.dirname(args.datasets.name)

    sys.stdout.write('Validation\n')
    log.write('# Validation\n')
    for individual, ind_datasets in datasets:
        sys.stdout.write('%s: ' % individual)
        for name, path in ind_datasets:
            sys.stdout.write('%s ' % name)
            sys.stdout.flush()

            command = [os.path.join(script_path, 'memory_validation.py'),
                       '-v', os.path.join(datasets_path, path),
                       '-c', args.combinations.name,
                       '-n', '%s-%s' % (individual, name),
                       '-d', os.path.join(dir, individual, name),
                       '--range', args.range[0], args.range[1],
                       '--neigh', args.neigh,
                       '--neigh-overhead', args.neigh_overhead,
                       '--neigh-good', args.neigh_good,
                       '--neigh-bad', args.neigh_bad]
            command = list(map(str, command))
            log.write('\t%s\n' % ' '.join(command))
            if subprocess.run(command).returncode:
                rc_fail()
        sys.stdout.write('\n')
    log.write('\n')


def combined_validate(args, datasets, log):
    dir = args.output
    script_path = os.path.dirname(__file__)
    datasets_path = os.path.dirname(args.datasets.name)

    sys.stdout.write('Combining individual datasets\n')
    log.write('# Combining individual datasets\n')
    for individual, ind_datasets in datasets:
        sys.stdout.write('%s ' % individual)
        sys.stdout.flush()

        cat_command = ['cat'] + list(map(itemgetter(1), ind_datasets))
        command = [os.path.join(script_path, 'memory_validation.py'),
                   '-v', '-',
                   '-c', args.combinations.name,
                   '-n', individual,
                   '-d', os.path.join(dir, individual),
                   '--range', args.range[0], args.range[1],
                   '--neigh', args.neigh,
                   '--neigh-overhead', args.neigh_overhead,
                   '--neigh-good', args.neigh_good,
                   '--neigh-bad', args.neigh_bad]
        command = list(map(str, command))
        log.write('\t%s | %s\n' % (' '.join(cat_command), ' ' .join(command)))

        p1 = subprocess.Popen(cat_command, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(command, stdin=p1.stdout)
        p1.stdout.close()
        p2.communicate()
        if p2.returncode:
            rc_fail()
    sys.stdout.write('\n')
    log.write('\n')


def tail_command(path, log, outp):
    command = ['tail', '-n+5', path]
    log.write('\t%s >> %s\n' % (' '.join(command), outp.name))
    p = subprocess.Popen(command, stdout=outp)
    p.communicate()


def combine_segments_summary(args, datasets, log):
    dir = args.output
    log.write('# Combine all segments_summary.csv\n')
    target = 'segments_summary.csv'
    with open(os.path.join(dir, target), 'w') as outp:
        outp.write('# %s\n' % ' '.join(sys.argv))
        outp.flush()
        command1 = ['tail', '-n+2', os.path.join(dir, datasets[0][0], target)]
        command2 = ['head', '-n3']
        log.write('\t%s | %s > %s\n' % (' '.join(command1), ' '.join(command2), outp.name))
        p1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(command2, stdin=p1.stdout, stdout=outp)
        p1.stdout.close()
        p2.communicate()
        if p2.returncode:
            rc_fail()

        for individual, ind_datasets in datasets:
            tail_command(os.path.join(dir, individual, target), log, outp)
            for name, _ in ind_datasets:
                tail_command(os.path.join(dir, individual, name, target), log, outp)


def run(args, human_args):
    dir = args.output
    mkdir(dir)
    if not args.force and os.listdir(dir):
        sys.stderr.write('Output folder is not empty. Try -f/--force\n')
        exit(1)

    with open(os.path.join(dir, 'log.txt'), 'w') as log:
        log_header(log, args, human_args)

        datasets = load_datasets(args.datasets, log)
        for individual, ind_datasets in datasets:
            mkdir(os.path.join(dir, individual))
            for name, path in ind_datasets:
                mkdir(os.path.join(dir, individual, name))

        validate(args, datasets, log)
        combined_validate(args, datasets, log)
        combine_segments_summary(args, datasets, log)
        sys.stdout.write(':::: Success ::::\n')


def main():
    parser = argparse.ArgumentParser(description='Use memory Rep-seq to validate novel segments',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='%(prog)s -c File -d File -o Dir [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-c', '--combinations', help='Combinations csv file', metavar='File',
                         type=argparse.FileType(), required=True)
    io_args.add_argument('-d', '--datasets',
                         help='File that specifies paths to datasets from\n'
                              'different individuals. Format: line <name>\n'
                              'indicates beginning of a single individual.\n'
                              'Following lines <name> <path to v_alignment fasta>\n'
                              'specify its datasets',
                         metavar='File', type=argparse.FileType(), required=True)
    io_args.add_argument('-o', '--output', help='Output directory', metavar='Dir', required=True)
    io_args.add_argument('-f', '--force', help='Override files in output directory', action='store_true')
    human_args = [('Input/output arguments', None), ('Combinations', 'combinations'), ('Datasets', 'datasets'),
                  ('Output directory', 'output'), ('Force', 'force')]

    val_args = parser.add_argument_group('Validating arguments')
    val_args.add_argument('--range', help='Positions range (default: [60, 290])',
                          metavar=('Int', 'Int'), nargs=2, default=[60, 290], type=int)
    val_args.add_argument('--neigh', help='Position neighborhood size (default: 3)',
                          metavar='Int', default=3, type=int)
    val_args.add_argument('--neigh-overhead', help='Acceptable mismatch rate overhead\n'
                                                   'in the position neighborhood (default: 1.1)',
                          metavar='Float', default=1.1, type=float, dest='neigh_overhead')
    val_args.add_argument('--neigh-good', help='Position is treated good, if its mismatch rate is\n'
                                               'less than --neigh-good (default: 0.2)',
                          metavar='Float', default=0.2, type=float, dest='neigh_good')
    val_args.add_argument('--neigh-bad', help='Position is treated good, if its mismatch rate is\n'
                                              'more than --neigh-bad (default: 0.8)',
                          metavar='Float', default=0.8, type=float, dest='neigh_bad')
    human_args += [('Validation arguments', None), ('Range', 'range'), ('neighborhood size', 'neigh'),
                   ('Neighborhood overhead', 'neigh_overhead'), ('Neighborhood: good', 'neigh_good'),
                   ('Neighborhood: bad', 'neigh_bad')]

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other.add_argument('--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    try:
        run(args, human_args)
    except KeyboardInterrupt:
        sys.stderr.write('\nKeyboard Interrupt\n')
        exit(1)


if __name__ == '__main__':
    main()
