#!/usr/bin/env python3

import sklearn.linear_model
import math
import numpy as np
import re
import random
try:
    import matplotlib
    matplotlib.use('Agg')
except ImportError:
    pass


class Segment:
    def __init__(self, dataset, line):
        line = line.strip().split()
        self.dataset = dataset
        self.segment = line[1]
        self.coverage = int(line[2])
        self.labels = int(line[3])
        self.mismatches = int(line[5])
        self.indels = int(line[6])
        self.consensus = line[7]
        self.full_seq = line[8]
        self.cropped_seq = line[9]

    def known(self):
        assert (self.mismatches + self.indels == 0) == (self.consensus == '*')
        return self.consensus == '*'

    def fitness(self):
        return self.coverage + self.labels

    def to_str(self, prob):
        return '{dataset}\t{segment}\t{prob:.4f}\t{coverage}\t{labels}\t{mismatches}\t{indels}\t{consensus}\t' \
               '{full}\t{cropped}'.format(dataset=self.dataset, segment=self.segment, prob=prob,
                    coverage=self.coverage, labels=self.labels, mismatches=self.mismatches, indels=self.indels,
                    consensus=self.consensus, full=self.full_seq, cropped=self.cropped_seq)


class Group:
    coverage_fun = np.log
    coverage_revfun = np.exp
    labels_fun = np.log
    labels_revfun = np.exp

    def __init__(self, name):
        self.name = name
        self.writtable_name = re.sub(r'\s', '_', self.name)
        self.segments = []
        self.probabilities = None
        self.lr = sklearn.linear_model.LogisticRegression(penalty='l2', C=100)

    def add_segment(self, segment):
        self.segments.append(segment)

    def classify(self):
        X = [[Group.coverage_fun(segm.coverage), Group.labels_fun(segm.labels)] for segm in self.segments]
        y = [int(segm.known()) for segm in self.segments]
        self.lr.fit(X, y)
        self.probabilities = self.lr.predict_proba(X)

    def get_segm_prob(self):
        assert len(self.probabilities) == len(self.segments)
        for segm, prob in zip(self.segments, self.probabilities):
            yield segm, prob[1]

    def add_isoclines(self, prob, isoclines, transformed):
        import matplotlib.pyplot as plt
        if transformed: # If the plot is transformed - isoclines are not
            x_fun = x_revfun = y_revfun = y_fun = lambda x: x
        else:
            x_fun = Group.coverage_fun
            y_fun = Group.labels_fun
            x_revfun = Group.coverage_revfun
            y_revfun = Group.labels_revfun

        axes = plt.axis()
        iso_x = x_revfun(np.linspace(*x_fun([max(axes[0], .001), axes[1]]), 100))
        for p in isoclines:
            interc = self.lr.intercept_[0] - math.log(p / (1 - p))
            coeffs = self.lr.coef_[0]

            iso_y = y_revfun(-(interc + coeffs[0] * x_fun(iso_x)) / coeffs[1])
            plt.plot(iso_x, iso_y, alpha=1 if p == prob else .4, linewidth=1, c='k') #, label=str(p))
            y_text = axes[2] + .1 * random.random() * (axes[3] - axes[2])
            x_text = x_revfun(-(interc + coeffs[1] * y_fun(y_text)) / coeffs[0])
            plt.text(x_text, y_text, str(p), fontsize=8, alpha=.5)

    def draw(self, outp_name, prob, isoclines, transformed=False):
        import matplotlib.pyplot as plt
        import matplotlib.colors as plt_col
        from operator import attrgetter
        from itertools import compress

        if transformed:
            x_fun = Group.coverage_fun
            x_revfun = Group.coverage_revfun
            y_fun = Group.labels_fun
            y_revfun = Group.labels_revfun
        else:
            x_fun = x_revfun = y_fun = y_revfun = lambda x: x

        x = list(map(x_fun, map(attrgetter('coverage'), self.segments)))
        y = list(map(y_fun, map(attrgetter('labels'), self.segments)))
        known = list(map(Segment.known, self.segments))

        colors = [plt_col.to_hex('xkcd:deep sky blue'), plt_col.to_hex('xkcd:tomato')]
        plt.clf()
        plt.scatter(list(compress(x, known)), list(compress(y, known)),
                    c=colors[0], s=5, alpha=.5, label='Known segments')
        plt.scatter(list(compress(x, np.logical_not(known))), list(compress(y, np.logical_not(known))),
                    c=colors[1], s=5, alpha=.5, label='Putative novel segments')
        axes = plt.axis()
        self.add_isoclines(prob, isoclines, transformed)
       
        plt.xlim([axes[0], axes[1]])
        plt.ylim([axes[2], axes[3]])
        plt.legend(loc=4)
        plt.title(self.name)

        fmt_x = '%s(%%s)' % x_fun.__name__ if transformed else '%s'
        fmt_y = '%s(%%s)' % y_fun.__name__ if transformed else '%s'
        plt.xlabel(fmt_x % 'Coverage')
        plt.ylabel(fmt_y % 'CDR3s')
        plt.savefig(outp_name, dpi=500)


def write_outp(groups, thr_prob, outp_success, outp_fail):
    import sys

    res = []
    for group in groups:
        for segm, prob in group.get_segm_prob():
            res.append((prob, group.writtable_name, segm))
    res.sort(reverse=True, key=lambda x: (x[0], x[2].fitness()))

    for outp in [outp_success, outp_fail]:
        outp.write('# %s\n' % ' '.join(sys.argv))
        outp.write('group\tdataset\tsegment\tprob\tcoverage\tlabels\tmismatches\t'
                   'indels\tconsensus\tfull_seq\tcropped_seq\n')

    for prob, group_name, segm in res:
        outp = outp_success if prob >= thr_prob else outp_fail
        outp.write('%s\t%s\n' % (group_name, segm.to_str(prob)))


def load_input(inputs, datasets):
    groups = []
    inputs = iter(inputs)
    for line in datasets:
        if line.startswith('#'):
            continue
        if line.startswith('!'):
            name = line.lstrip('!').strip()
            groups.append(Group(name))
            continue

        name = line.strip().split()[0]
        inp = next(inputs)
        while next(inp).startswith('#'):
            pass

        for line2 in inp:
            groups[-1].add_segment(Segment(name, line2))
    return groups


def draw_plots(groups, out_dir, prob, isoclines):
    import os

    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    n = len(groups)
    if n == 0:
        return
    fmt = '%s/group%%0%dd%%s.png' % (out_dir, 1 + int(math.log10(n)))
    
    isoclines = sorted(set(isoclines) | {prob}) if isoclines is not None else []
    for i, group in enumerate(groups):
        group.draw(fmt % (i + 1, ''), prob, isoclines, False)
        group.draw(fmt % (i + 1, '_transformed'), prob, isoclines, True)


def main():
    import argparse
    from extra._version import __version__

    np.seterr(all='ignore')

    parser = argparse.ArgumentParser(description='Combine and classify segments from multiple summary files',
                                     formatter_class=argparse.RawTextHelpFormatter, add_help=False,
                                     usage='-i File+ -d File -o File File [-p Dir] [args]')
    io_args = parser.add_argument_group('Input/output arguments')
    io_args.add_argument('-i', '--input', nargs='+', metavar='File', type=argparse.FileType(),
                         help='Multiple input csv file with summaries over found segments')
    io_args.add_argument('-d', '--datasets', metavar='File', type=argparse.FileType(),
                         help='Csv file where each line starts with the name of a dataset,\n'
                         'order is the same as in --input File+.\n'
                         'Datasets with different properties (significantly different coverage,\n'
                         'different PCR strength) are separated by lines "! Group name"')
    io_args.add_argument('-o', '--output', metavar=('File', 'File'), nargs=2, type=argparse.FileType('w'),
                         help='Two output csv files, with found segments and filtered out segments')
    io_args.add_argument('-p', '--plots', metavar='Dir', required=False,
                         help='Output directory with plots over each group (not required)')

    cl_args = parser.add_argument_group('Classification arguments')
    cl_args.add_argument('--prob', metavar='Float', type=float, default=.8,
                         help='Minimal logistic regression probability (default: %(default)s)')

    dr_args = parser.add_argument_group('Drawing arguments')
    dr_args.add_argument('--isoclines', nargs='*', type=float, metavar='Float',
                         help='Draw isoclines with the provided probabilities, \n'
                         'if --isoclines is provided, --prob isocline is always included')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    groups = load_input(args.input, args.datasets)
    for group in groups:
        group.classify()
    write_outp(groups, args.prob, *args.output)
    if args.plots:
        draw_plots(groups, args.plots, args.prob, args.isoclines)


if __name__ == '__main__':
    main()
