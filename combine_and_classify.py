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
import sklearn.metrics


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

    def to_str(self, dist_to_break):
        return '{dataset}\t{segment}\t{dist_to_break:.1f}\t{coverage}\t{coverage_rate}\t{labels}\t{mismatches}\t' \
               '{indels}\t{consensus}\t{full}\t{cropped}'.format(dataset=self.dataset, segment=self.segment,
                       dist_to_break=dist_to_break, coverage=self.coverage, coverage_rate=self.coverage_rate,
                       labels=self.labels, mismatches=self.mismatches, indels=self.indels,
                       consensus=self.consensus, full=self.full_seq, cropped=self.cropped_seq)


def count_elements_to_left(projections, coef, interc):
    x_min = min(projections)
    x_max = max(projections)

    total = len(projections)

    res = []
    for i in range(101):
        res.append(1 - sum(projections <= x_min + (x_max - x_min) * i / 100) / total)
    return res


def lr_squared_error(curve, i=0, j=None):
    if j is None:
        j = len(curve)
    x = np.matrix(list(range(i, j))).T
    y = np.array(curve[i:j])

    lr = sklearn.linear_model.LinearRegression().fit(x, y)
    pred = lr.predict(x)
    return (j - i) * sklearn.metrics.mean_squared_error(pred, y)


def break_curve(curve):
    m = float('inf')
    best_i = 0
    for i in range(1, len(curve) - 1):
        curr = lr_squared_error(curve, 0, i) + lr_squared_error(curve, i)
        if curr < m:
            m = curr
            best_i = i
    return best_i, m


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
        projections = self.project(self.coef, self.interc)
        unknown_projections = np.array([x for x, segm in zip(projections, self.segments) if not segm.known()])
        curve = count_elements_to_left(unknown_projections, self.coef, self.interc)
        break_point, _ = break_curve(curve)
        
        if break_point <= 5:
            self.interc = None
            self.coef = None
            break_point = 100

        x_min = min(unknown_projections)
        x_max = max(unknown_projections)
        self.break_point_proj = break_point / 100 * (x_max - x_min) + x_min
        self.step = (x_max - x_min) / 100
        self.dist_to_break = (projections - self.break_point_proj) / self.step

    def get_segm_dist_to_break(self):
        assert len(self.dist_to_break) == len(self.segments)
        return zip(self.segments, self.dist_to_break)

    def add_isoclines(self, dist_to_break, isoclines):
        import matplotlib.pyplot as plt
        axes = plt.axis()

        x = np.array([axes[0], axes[1]])
        y = self.interc + x * self.coef
        plt.plot(x, y, alpha=.5, c='black')

        for dist in isoclines:
            x0 = dist * self.step + self.break_point_proj
            y0 = self.interc + x0 * self.coef
            y = np.array([axes[2], axes[3]])
            x = -self.coef * y + self.coef * y0 + x0
            plt.plot(x, y, alpha=.5, c='green' if dist == dist_to_break else 'black')

    def draw(self, outp_name, prob, isoclines):
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
        self.add_isoclines(prob, isoclines)
       
        plt.xlim([axes[0], axes[1]])
        plt.ylim([axes[2], axes[3]])
        plt.legend(loc=4)
        plt.title(self.name)

        plt.xlabel('log(Coverage rate)')
        plt.ylabel('log(#CDR3s)')
        plt.savefig(outp_name, dpi=500)

    def print_coefficients(self):
        print('%s  intercept   coverage' % self.name)
        print(' ' * len(self.name), end='  ')
        print('%9.4f' % self.interc, end='  ')
        print('%9.4f' % self.coef)


def write_outp(groups, thr_dist, outp_success, outp_fail):
    import sys

    res = []
    for group in groups:
        for segm, dist_to_break in group.get_segm_dist_to_break():
            res.append((group.writtable_name, segm, dist_to_break))
    res.sort(reverse=True, key=lambda x: (x[2], x[1].fitness()))

    for outp in [outp_success, outp_fail]:
        outp.write('# %s\n' % ' '.join(sys.argv))
        outp.write('group\tdataset\tsegment\tdist_to_break\tcoverage\tcoverage_rate\t'
                   'labels\tmismatches\tindels\tconsensus\tfull_seq\tcropped_seq\n')

    for group_name, segm, dist_to_break in res:
        outp = outp_success if dist_to_break >= thr_dist else outp_fail
        outp.write('%s\t%s\n' % (group_name, segm.to_str(dist_to_break)))


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


def draw_plots(groups, out_dir, thr_dist):
    import os

    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    n = len(groups)
    if n == 0:
        return
    fmt = '%s/group%%0%dd.png' % (out_dir, 1 + int(math.log10(n)))
    
    isoclines = [0, thr_dist]
    for i, group in enumerate(groups):
        group.draw(fmt % (i + 1), thr_dist, isoclines)


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
    cl_args.add_argument('-t', '--threshold', metavar='Float', type=float, default=10,
                         help='Minimal distance to the break point (default: %(default)s)')
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
            group.print_coefficients()


if __name__ == '__main__':
    main()
