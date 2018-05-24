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

import sklearn.linear_model
import scipy.stats


class Segment:
    def __init__(self, dataset, line):
        line = line.strip().split()
        self.dataset = dataset
        self.segment = line[1]
        self.coverage = int(line[2])
        self.coverage_rate = float(line[3])
        self.labels = int(line[4])
        self.mismatches = int(line[6])
        self.indels = int(line[7])
        self.consensus = line[8]
        self.full_seq = line[9]
        self.cropped_seq = line[10]

    def known(self):
        assert (self.mismatches + self.indels == 0) == (self.consensus == '*')
        return self.consensus == '*'

    def fitness(self):
        return self.coverage_rate, self.labels

    def to_str(self, pvalue):
        return '{dataset}\t{segment}\t{pvalue}\t{coverage}\t{coverage_rate}\t{labels}\t{mismatches}\t' \
               '{indels}\t{consensus}\t{full}\t{cropped}'.format(dataset=self.dataset, segment=self.segment,
                       pvalue=pvalue, coverage=self.coverage, coverage_rate=self.coverage_rate,
                       labels=self.labels, mismatches=self.mismatches, indels=self.indels,
                       consensus=self.consensus, full=self.full_seq, cropped=self.cropped_seq)


class Group:
    def __init__(self, name):
        self.name = name
        self.writtable_name = re.sub(r'\s', '_', self.name)
        self.segments = []
        self.probabilities = None
        self.lr = sklearn.linear_model.LogisticRegression(penalty='l2', C=100)

    def add_segment(self, segment):
        self.segments.append(segment)

    def regression_coefficients(self):
        x = np.log(np.matrix([segm.coverage_rate for segm in self.segments]).T)
        y = np.log([segm.labels for segm in self.segments])
        lr = sklearn.linear_model.LinearRegression()
        lr.fit(x, y)
        return lr.coef_[0], lr.intercept_

    def project(self, coef, interc):
        x = np.log([segm.coverage_rate for segm in self.segments])
        y = np.log([segm.labels for segm in self.segments])
        return (y + x / coef - interc) / (coef + 1 / coef)

    def classify(self):
        self.coef, self.interc = self.regression_coefficients()
        self.projections = self.project(self.coef, self.interc)
        unknown_projections = np.array([x for x, segm in zip(self.projections, self.segments) if not segm.known()])

        self.mu = np.median(unknown_projections)
        self.sigma = np.std(unknown_projections)
        self.pvalue = 1 - scipy.stats.norm.cdf(self.projections, self.mu, self.sigma)

    def get_segm_pvalue(self):
        assert len(self.pvalue) == len(self.segments)
        return zip(self.segments, self.pvalue)

    def draw(self, outp_name, pvalue):
        import matplotlib.pyplot as plt
        import matplotlib.colors as plt_col
        from itertools import compress

        x = np.log([segm.coverage_rate for segm in self.segments])
        y = np.log([segm.labels for segm in self.segments])
        known = list(map(Segment.known, self.segments))

        colors = [plt_col.to_hex('xkcd:deep sky blue'), plt_col.to_hex('xkcd:tomato')]
        plt.clf()
        plt.scatter(list(compress(x, known)), list(compress(y, known)),
                    c=colors[0], s=5, alpha=.5, label='Known segments')
        plt.scatter(list(compress(x, np.logical_not(known))), list(compress(y, np.logical_not(known))),
                    c=colors[1], s=5, alpha=.5, label='Putative novel segments')
        axes = plt.axis()

        x = np.array([axes[0], axes[1]])
        y = self.interc + x * self.coef
        plt.plot(x, y, alpha=.5, c='black')

        x0 = scipy.stats.norm.ppf(1 - pvalue, self.mu, self.sigma)
        y0 = self.interc + x0 * self.coef
        y = np.array([axes[2], axes[3]])
        x = -self.coef * y + self.coef * y0 + x0 
        plt.plot(x, y, alpha=.5, c='black')

#        plt.xlim([axes[0], axes[1]])
#        plt.ylim([axes[2], axes[3]])
        plt.legend(loc=4)
        plt.title(self.name)

        plt.xlabel('log(Coverage rate)')
        plt.ylabel('log(#CDR3s)')
        plt.savefig(outp_name, dpi=500)

    def draw_hist(self, outp_name, pvalue):
        import matplotlib.pyplot as plt
        import matplotlib.colors as plt_col
        from itertools import compress

        known = list(map(Segment.known, self.segments))
        projections = (self.projections - self.mu) / self.sigma

        colors = [plt_col.to_hex('xkcd:deep sky blue'), plt_col.to_hex('xkcd:tomato')]
        plt.clf()
        plt.hist(list(compress(projections, known)),
                 color=colors[0], alpha=.5, label='Known segments', bins=50)
        count, bins, _ = plt.hist(list(compress(projections, np.logical_not(known))),
                                  color=colors[1], alpha=.5, label='Putative novel segments', bins=50)

        normal = scipy.stats.norm.pdf(bins)
        plt.plot(bins, normal / sum(normal) * sum(count), linewidth=2)
        plt.plot([scipy.stats.norm.ppf(1 - pvalue)] * 2, [0, max(count)],
                 c='black', alpha=.5)

        plt.legend(loc=1)
        plt.title(self.name)

        plt.xlabel('Segment projection')
        plt.ylabel('Count')
        plt.savefig(outp_name, dpi=500)

    def print_coefficients(self, pvalue):
        print('%s\tintercept\tcoverage' % self.name)
        print('\t%.5f\t%.5f' % (self.interc, self.coef))

        x0 = scipy.stats.norm.ppf(1 - pvalue, self.mu, self.sigma)
        y0 = self.interc + x0 * self.coef
        print('\t%.5f\t%.5f' % (x0 / self.coef + y0, -1 / self.coef))


def write_outp(groups, thr_pvalue, outp_success, outp_fail):
    import sys

    res = []
    for group in groups:
        for segm, pvalue in group.get_segm_pvalue():
            res.append((group.writtable_name, segm, pvalue))
    res.sort(reverse=True, key=lambda x: (x[2], x[1].fitness()))

    for outp in [outp_success, outp_fail]:
        outp.write('# %s\n' % ' '.join(sys.argv))
        outp.write('group\tdataset\tsegment\tpvalue\tcoverage\tcoverage_rate\t'
                   'labels\tmismatches\tindels\tconsensus\tfull_seq\tcropped_seq\n')

    for group_name, segm, pvalue in res:
        outp = outp_success if pvalue >= thr_pvalue else outp_fail
        outp.write('%s\t%s\n' % (group_name, segm.to_str(pvalue)))


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


def draw_plots(groups, out_dir, pvalue):
    import os

    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    n = len(groups)
    if n == 0:
        return
    fmt = '%s/group%%0%dd%%s.png' % (out_dir, 1 + int(math.log10(n)))
    
    for i, group in enumerate(groups):
        group.draw(fmt % (i + 1, ''), pvalue)
        group.draw_hist(fmt % (i + 1, '_hist'), pvalue)


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
    cl_args.add_argument('-t', '--threshold', metavar='Float', type=float, default=0.05,
                         help='Maximal p-value (default: %(default)s)')
    cl_args.add_argument('--verbose', action='store_true', help='Print classification coefficients')

    other = parser.add_argument_group('Other arguments')
    other.add_argument('-h', '--help', action='help', help='Show this message and exit')
    other.add_argument('-V', '--version', action='version', help='Show version', version=__version__)

    args = parser.parse_args()

    groups = load_input(args.input, args.datasets)
    for group in groups:
        group.classify()
    write_outp(groups, args.threshold, *args.output)
    if args.plots:
        draw_plots(groups, args.plots, args.threshold)

    if args.verbose:
        for group in groups:
            group.print_coefficients(args.threshold)


if __name__ == '__main__':
    main()
