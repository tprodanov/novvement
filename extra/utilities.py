import sys
from operator import itemgetter
from termcolor import colored

class Tee:
    def __init__(self, f, g=sys.stdout):
        self.f = f
        self.g = g

    def write(self, *args, **kwargs):
        self.f.write(*args, **kwargs)
        self.g.write(*args, **kwargs)


def v_alignment_to_mismatches(v_alignment, potential_mismatches):
    while next(v_alignment).startswith('#'):
        pass
    prev_read = None

    for line in v_alignment:
        segment, read, mut_type, position, nt = line.strip().split()
        position = int(position)

        if segment not in potential_mismatches:
            continue
        if prev_read != read:
            if prev_read and read_mismatches:
                read_mismatches.sort(key=itemgetter(0))
                yield prev_read, prev_segment, read_mismatches
            prev_read = read
            prev_segment = segment
            read_mismatches = []
            segment_mismatches = potential_mismatches[segment]

        if (position, nt) in segment_mismatches:
            read_mismatches.append((position, nt))
    if prev_read and read_mismatches:
        read_mismatches.sort(key=itemgetter(0))
        yield prev_read, prev_segment, read_mismatches


def print_colored(msg, stdout, col, attrs):
    (sys.stdout if stdout else sys.stderr).write(colored(msg, col, attrs=attrs))


def oprint(msg, col, *attrs):
    print_colored(msg, True, col, attrs)


def eprint(msg, col, *attrs):
    print_colored(msg, False, col, attrs)


