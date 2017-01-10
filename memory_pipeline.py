#!/usr/bin/python3

import argparse
import sys
import os
import subprocess
from operator import itemgetter


def mkdir(path):
    try:
        os.mkdir(path)
    except FileExistsError:
        pass


def log_header(log, args):
    log.write('# Command: %s\n' % ' '.join(sys.argv))
    log.write('#\tCombinations: %s\n' % args.combinations.name)
    log.write('#\tDatasets: %s\n' % args.datasets.name)
    log.write('#\tOutput: %s\n' % args.output)
    log.write('#\tForce: %s\n' % args.force)
    log.write('\n\n')
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
                       '-d', os.path.join(dir, individual, name)]
            log.write('\t%s\n' % ' '.join(command))
            subprocess.run(command)
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
                   '-d', os.path.join(dir, individual)]
        log.write('\t%s | %s\n' % (' '.join(cat_command), ' ' .join(command)))

        p1 = subprocess.Popen(cat_command, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(command, stdin=p1.stdout)
        p1.stdout.close()
        p2.communicate()
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

        for individual, ind_datasets in datasets:
            tail_command(os.path.join(dir, individual, target), log, outp)
            for name, _ in ind_datasets:
                tail_command(os.path.join(dir, individual, name, target), log, outp)


def run(args):
    dir = args.output
    mkdir(dir)
    if not args.force and os.listdir(dir):
        sys.stderr.write('Output folder is not empty. Try -f/--force\n')
        exit(1)

    with open(os.path.join(dir, 'log.txt'), 'w') as log:
        log_header(log, args)

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
    parser = argparse.ArgumentParser(description='use memory Rep-seq to validate novel segments')
    parser.add_argument('-c', '--combinations', help='combinations csv file', metavar='FILE',
                        type=argparse.FileType(), required=True)
    parser.add_argument('-d', '--datasets', help='datasets file, that specifies paths to datasets from different individuals. '
                                                 'Format: line <name> indicates beginning of a single individual. '
                                                 'Following lines <name> <path to v_alignment fasta> specify this individual datasets. ',
                        metavar='FILE', type=argparse.FileType(), required=True)
    parser.add_argument('-o', '--output', help='output directory', metavar='DIR', required=True)
    parser.add_argument('-f', '--force', help='override files in output directory', action='store_true')

    args = parser.parse_args()
    try:
        run(args)
    except KeyboardInterrupt:
        sys.stderr.write('\nKeyboard Interrupt\n')
        exit(1)


if __name__ == '__main__':
    main()

