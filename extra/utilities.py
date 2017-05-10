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


def print_colored(msg, stdout, col, attrs):
    stream = sys.stdout if stdout else sys.stderr
    stream.write(colored(msg, col, attrs=attrs))
    stream.flush()


def oprint(msg, col, *attrs):
    print_colored(msg, True, col, attrs)


def eprint(msg, col, *attrs):
    print_colored(msg, False, col, attrs)


